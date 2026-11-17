from __future__ import annotations

from dataclasses import dataclass, replace
import os
from pathlib import Path
import re
from typing import TYPE_CHECKING

from .executor import execute_manifest_with_resume
from .executor import RESUME_ACTION_USE_EXISTING
from .layout import (
    layout_from_analyze_manifest,
    layout_from_prepare_manifest,
    layout_from_run_manifest,
    resolve_project_output_root,
)
from .models import (
    PROJECT_MANIFEST_TYPE,
    AnalyzeWorkflowManifest,
    PrepareWorkflowManifest,
    ProjectStageExecutionRecord,
    ProjectWorkflowManifest,
    ProjectWorkflowResult,
    RunWorkflowManifest,
    StageFileRecord,
)
from .resume import STAGE_STATUS_COMPLETED

if TYPE_CHECKING:
    from .. import EVBSystem


PathLike = os.PathLike[str] | str
StageOutputManifest = PrepareWorkflowManifest | RunWorkflowManifest | AnalyzeWorkflowManifest
PROJECT_RESTART_POLICY_CASCADE = "cascade"
PROJECT_RESTART_POLICY_STAGE_LOCAL = "stage_local"
PROJECT_ACTION_SKIPPED = "skipped"
VALID_PROJECT_RESTART_POLICIES = frozenset(
    {
        PROJECT_RESTART_POLICY_CASCADE,
        PROJECT_RESTART_POLICY_STAGE_LOCAL,
    }
)
VALID_PROJECT_RECOVERY_START_STAGES = frozenset({"prepare", "run", "analyze"})


@dataclass(frozen=True)
class ProjectStagePlan:
    """Resolved stage input plan for one step of a project workflow."""

    stage_name: str
    label: str
    input_manifest_label: str
    input_manifest_path: Path
    manifest: StageOutputManifest


def execute_project_manifest(
    path: PathLike,
    *,
    force_rerun: bool = False,
    restart_policy: str = PROJECT_RESTART_POLICY_CASCADE,
    start_stage: str | None = None,
    analyze_label: str | None = None,
    evb_system_factory: type["EVBSystem"] | object = None,
) -> ProjectWorkflowManifest:
    """Load and execute a chained project workflow manifest from disk."""

    manifest = ProjectWorkflowManifest.read_json(path)
    _, output_manifest = run_project_workflow(
        manifest,
        force_rerun=force_rerun,
        restart_policy=restart_policy,
        start_stage=start_stage,
        analyze_label=analyze_label,
        evb_system_factory=evb_system_factory,
    )
    return output_manifest


def recover_project_manifest(
    path: PathLike,
    *,
    start_stage: str | None = None,
    analyze_label: str | None = None,
    evb_system_factory: type["EVBSystem"] | object = None,
) -> tuple[ProjectWorkflowManifest, str]:
    """Recover a project workflow by selecting a restart policy from current status."""

    _validate_recovery_selector_conflicts(start_stage=start_stage, analyze_label=analyze_label)

    if analyze_label is not None:
        normalized_analyze_label = _normalize_recovery_analyze_label(analyze_label)
        output_manifest = execute_project_manifest(
            path,
            analyze_label=normalized_analyze_label,
            evb_system_factory=evb_system_factory,
        )
        return output_manifest, f"analyze_label:{normalized_analyze_label}"

    if start_stage is not None:
        normalized_start_stage = _normalize_recovery_start_stage(start_stage)
        output_manifest = execute_project_manifest(
            path,
            start_stage=normalized_start_stage,
            evb_system_factory=evb_system_factory,
        )
        return output_manifest, f"from_{normalized_start_stage}"

    restart_policy = suggest_project_restart_policy(path)
    output_manifest = execute_project_manifest(
        path,
        restart_policy=restart_policy,
        evb_system_factory=evb_system_factory,
    )
    return output_manifest, restart_policy


def run_project_workflow(
    manifest: ProjectWorkflowManifest,
    *,
    force_rerun: bool = False,
    restart_policy: str = PROJECT_RESTART_POLICY_CASCADE,
    start_stage: str | None = None,
    analyze_label: str | None = None,
    evb_system_factory: type["EVBSystem"] | object = None,
) -> tuple[tuple[StageOutputManifest, ...], ProjectWorkflowManifest]:
    """Execute a chained prepare/run/analyze workflow using stage-local resume semantics."""

    project_root = _resolve_project_workflow_root(manifest)
    stage_plans = resolve_project_stage_manifests(manifest)
    normalized_restart_policy = _normalize_restart_policy(restart_policy)
    normalized_start_stage = _normalize_recovery_start_stage(start_stage)
    normalized_analyze_label = _normalize_recovery_analyze_label(analyze_label)

    _validate_recovery_selector_conflicts(
        start_stage=normalized_start_stage,
        analyze_label=normalized_analyze_label,
    )
    if normalized_analyze_label is not None:
        _validate_target_analyze_label(stage_plans, normalized_analyze_label)

    stage_outputs: list[StageOutputManifest] = []
    stage_results: list[ProjectStageExecutionRecord] = []
    resolved_input_manifests: list[StageFileRecord] = []

    for stage_plan in stage_plans:
        stage_plan.input_manifest_path.parent.mkdir(parents=True, exist_ok=True)
        stage_plan.manifest.write_json(stage_plan.input_manifest_path)
        resolved_input_manifests.append(
            StageFileRecord(stage_plan.input_manifest_label, stage_plan.input_manifest_path)
        )

    for stage_plan in stage_plans:
        if _should_skip_stage_for_analyze_label(
            stage_plan,
            analyze_label=normalized_analyze_label,
        ):
            stage_results.append(
                ProjectStageExecutionRecord(
                    stage_name=stage_plan.stage_name,
                    label=stage_plan.label,
                    input_manifest_path=stage_plan.input_manifest_path,
                    output_manifest_path=_project_stage_output_path(manifest, stage_plan),
                    action=PROJECT_ACTION_SKIPPED,
                )
            )
            continue

        stage_force_rerun = force_rerun or _should_force_rerun_from_stage(
            stage_plan,
            start_stage=normalized_start_stage,
        ) or _should_force_rerun_for_analyze_label(
            stage_plan,
            analyze_label=normalized_analyze_label,
        ) or (
            normalized_start_stage is None
            and normalized_analyze_label is None
            and normalized_restart_policy == PROJECT_RESTART_POLICY_CASCADE
            and _should_force_rerun(stage_results)
        )
        output_manifest, action = execute_manifest_with_resume(
            stage_plan.input_manifest_path,
            force_rerun=stage_force_rerun,
            evb_system_factory=evb_system_factory,
        )
        stage_outputs.append(output_manifest)
        stage_results.append(
            ProjectStageExecutionRecord(
                stage_name=stage_plan.stage_name,
                label=stage_plan.label,
                input_manifest_path=stage_plan.input_manifest_path,
                output_manifest_path=_manifest_output_path(output_manifest),
                action=action,
            )
        )

    result = ProjectWorkflowResult(
        stage_results=tuple(stage_results),
        resolved_input_manifests=tuple(resolved_input_manifests),
        manifest_path=project_root / "manifests" / "project_manifest.json",
    )
    output_manifest = ProjectWorkflowManifest(
        project=manifest.project,
        spec=manifest.spec,
        result=result,
    )
    output_manifest.write_json(result.manifest_path)
    return tuple(stage_outputs), output_manifest


def resolve_project_stage_manifests(manifest: ProjectWorkflowManifest) -> tuple[ProjectStagePlan, ...]:
    """Resolve the chained stage input manifests for a project workflow."""

    project_root = _resolve_project_workflow_root(manifest)
    manifests_dir = project_root / "manifests"
    stage_plans: list[ProjectStagePlan] = []

    prepare_output_manifest_path: Path | None = None
    run_output_manifest_path: Path | None = None

    if manifest.spec.prepare_manifest_path is not None:
        prepare_template = PrepareWorkflowManifest.read_json(manifest.spec.prepare_manifest_path)
        prepare_manifest = replace(prepare_template, project=manifest.project, result=None)
        stage_plans.append(
            ProjectStagePlan(
                stage_name="prepare",
                label="prepare",
                input_manifest_label="prepare_input_manifest",
                input_manifest_path=manifests_dir / "prepare_input_manifest.json",
                manifest=prepare_manifest,
            )
        )
        prepare_output_manifest_path = _default_output_manifest_path(prepare_manifest)

    if manifest.spec.run_manifest_path is not None:
        run_template = RunWorkflowManifest.read_json(manifest.spec.run_manifest_path)
        run_spec = run_template.spec
        if prepare_output_manifest_path is not None:
            run_spec = replace(run_spec, prepared_manifest_path=prepare_output_manifest_path)
        run_manifest = replace(run_template, project=manifest.project, spec=run_spec, result=None)
        stage_plans.append(
            ProjectStagePlan(
                stage_name="run",
                label=str(run_manifest.spec.workflow_kind),
                input_manifest_label="run_input_manifest",
                input_manifest_path=manifests_dir / "run_input_manifest.json",
                manifest=run_manifest,
            )
        )
        run_output_manifest_path = _default_output_manifest_path(run_manifest)

    for index, analyze_manifest_file in enumerate(manifest.spec.analyze_manifest_files):
        analyze_template = AnalyzeWorkflowManifest.read_json(analyze_manifest_file.path)
        analyze_spec = analyze_template.spec
        if run_output_manifest_path is not None:
            analyze_spec = replace(analyze_spec, run_manifest_path=run_output_manifest_path)
        analyze_manifest = replace(analyze_template, project=manifest.project, spec=analyze_spec, result=None)
        stage_plans.append(
            ProjectStagePlan(
                stage_name="analyze",
                label=analyze_manifest_file.label,
                input_manifest_label=f"analyze_input_manifest_{analyze_manifest_file.label}",
                input_manifest_path=manifests_dir / _resolved_analyze_manifest_name(index, analyze_manifest_file.label),
                manifest=analyze_manifest,
            )
        )

    return tuple(stage_plans)


def suggest_project_restart_policy(target: PathLike | "WorkflowStatusReport") -> str:
    """Choose a project restart policy from the current stage-status view."""

    report = _coerce_status_report(target)
    if report.manifest_type != PROJECT_MANIFEST_TYPE:
        raise ValueError(
            f"Project restart policy requires manifest_type={PROJECT_MANIFEST_TYPE!r}, "
            f"received {report.manifest_type!r}."
        )

    first_incomplete = next(
        (
            summary
            for summary in report.stage_statuses
            if summary.status != STAGE_STATUS_COMPLETED
        ),
        None,
    )
    if first_incomplete is None:
        return PROJECT_RESTART_POLICY_STAGE_LOCAL
    if first_incomplete.stage_name == "analyze":
        return PROJECT_RESTART_POLICY_STAGE_LOCAL
    return PROJECT_RESTART_POLICY_CASCADE


def _resolve_project_workflow_root(manifest: ProjectWorkflowManifest) -> Path:
    candidate = manifest.spec.output_directory.expanduser()
    if candidate.is_absolute():
        return candidate
    return resolve_project_output_root(manifest.project) / candidate


def _manifest_output_path(manifest: StageOutputManifest) -> Path:
    if manifest.result is None or manifest.result.manifest_path is None:
        raise ValueError(f"{type(manifest).__name__} did not produce an output manifest path.")
    return manifest.result.manifest_path


def _default_output_manifest_path(manifest: StageOutputManifest) -> Path:
    if isinstance(manifest, PrepareWorkflowManifest):
        return layout_from_prepare_manifest(manifest).manifest_path()
    if isinstance(manifest, RunWorkflowManifest):
        return layout_from_run_manifest(manifest).manifest_path()
    if isinstance(manifest, AnalyzeWorkflowManifest):
        return layout_from_analyze_manifest(manifest).manifest_path()
    raise TypeError(f"Unsupported workflow manifest type: {type(manifest).__name__}.")


def _project_stage_output_path(
    manifest: ProjectWorkflowManifest,
    stage_plan: ProjectStagePlan,
) -> Path:
    if manifest.result is not None:
        for stage_result in manifest.result.stage_results:
            if (
                stage_result.stage_name == stage_plan.stage_name
                and stage_result.label == stage_plan.label
            ):
                return stage_result.output_manifest_path
    return _default_output_manifest_path(stage_plan.manifest)


def _resolved_analyze_manifest_name(index: int, label: str) -> str:
    return f"analyze_{index:02d}_{_slug(label)}_input_manifest.json"


def _slug(value: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9_-]+", "_", value.strip())
    normalized = normalized.strip("_-")
    return normalized or "analyze"


def _normalize_restart_policy(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized not in VALID_PROJECT_RESTART_POLICIES:
        raise ValueError(
            f"Unsupported project restart policy {value!r}. "
            f"Expected one of {tuple(sorted(VALID_PROJECT_RESTART_POLICIES))!r}."
        )
    return normalized


def _normalize_recovery_start_stage(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip().lower()
    if normalized not in VALID_PROJECT_RECOVERY_START_STAGES:
        raise ValueError(
            f"Unsupported project recovery start stage {value!r}. "
            f"Expected one of {tuple(sorted(VALID_PROJECT_RECOVERY_START_STAGES))!r}."
        )
    return normalized


def _normalize_recovery_analyze_label(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        raise ValueError("Analyze recovery label must not be empty.")
    return normalized


def _validate_recovery_selector_conflicts(
    *,
    start_stage: str | None,
    analyze_label: str | None,
) -> None:
    if start_stage is not None and analyze_label is not None:
        raise ValueError("Project recovery accepts only one of start_stage or analyze_label.")


def _validate_target_analyze_label(
    stage_plans: tuple[ProjectStagePlan, ...],
    analyze_label: str,
) -> None:
    matching_labels = [
        stage_plan
        for stage_plan in stage_plans
        if stage_plan.stage_name == "analyze" and stage_plan.label == analyze_label
    ]
    if len(matching_labels) == 1:
        return

    available_labels = sorted(
        {
            stage_plan.label
            for stage_plan in stage_plans
            if stage_plan.stage_name == "analyze"
        }
    )
    if not matching_labels:
        raise ValueError(
            f"Unknown analyze recovery label {analyze_label!r}. "
            f"Expected one of {tuple(available_labels)!r}."
        )
    raise ValueError(
        f"Analyze recovery label {analyze_label!r} is ambiguous across project stages."
    )


def _should_force_rerun(stage_results: list[ProjectStageExecutionRecord]) -> bool:
    return any(record.action != RESUME_ACTION_USE_EXISTING for record in stage_results)


def _should_force_rerun_from_stage(
    stage_plan: ProjectStagePlan,
    *,
    start_stage: str | None,
) -> bool:
    if start_stage is None:
        return False
    return _stage_order(stage_plan.stage_name) >= _stage_order(start_stage)


def _should_force_rerun_for_analyze_label(
    stage_plan: ProjectStagePlan,
    *,
    analyze_label: str | None,
) -> bool:
    if analyze_label is None:
        return False
    return stage_plan.stage_name == "analyze" and stage_plan.label == analyze_label


def _should_skip_stage_for_analyze_label(
    stage_plan: ProjectStagePlan,
    *,
    analyze_label: str | None,
) -> bool:
    if analyze_label is None:
        return False
    return stage_plan.stage_name == "analyze" and stage_plan.label != analyze_label


def _stage_order(stage_name: str) -> int:
    if stage_name == "prepare":
        return 0
    if stage_name == "run":
        return 1
    if stage_name == "analyze":
        return 2
    raise ValueError(f"Unsupported project stage name {stage_name!r}.")


def _coerce_status_report(target: PathLike | "WorkflowStatusReport") -> "WorkflowStatusReport":
    if hasattr(target, "stage_statuses") and hasattr(target, "manifest_type"):
        return target
    from .status import collect_workflow_status

    return collect_workflow_status(target)


__all__ = [
    "PROJECT_ACTION_SKIPPED",
    "PROJECT_RESTART_POLICY_CASCADE",
    "PROJECT_RESTART_POLICY_STAGE_LOCAL",
    "VALID_PROJECT_RECOVERY_START_STAGES",
    "VALID_PROJECT_RESTART_POLICIES",
    "ProjectStagePlan",
    "execute_project_manifest",
    "recover_project_manifest",
    "resolve_project_stage_manifests",
    "run_project_workflow",
    "suggest_project_restart_policy",
]
