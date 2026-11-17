from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
from typing import Any

from .layout import layout_from_analyze_manifest, layout_from_prepare_manifest, layout_from_run_manifest
from .models import (
    ANALYZE_MANIFEST_TYPE,
    PREPARE_MANIFEST_TYPE,
    PROJECT_MANIFEST_TYPE,
    RUN_MANIFEST_TYPE,
    AnalyzeWorkflowManifest,
    PrepareWorkflowManifest,
    ProjectWorkflowManifest,
    RunWorkflowManifest,
)
from .project import resolve_project_stage_manifests
from .resume import STAGE_STATUS_COMPLETED, read_stage_resume_state, stage_resume_state_path


PathLike = os.PathLike[str] | str
WorkflowManifest = PrepareWorkflowManifest | RunWorkflowManifest | AnalyzeWorkflowManifest | ProjectWorkflowManifest
StageWorkflowManifest = PrepareWorkflowManifest | RunWorkflowManifest | AnalyzeWorkflowManifest


class _SerializableRecord:
    def to_dict(self) -> dict[str, object]:
        return _jsonify(asdict(self))

    def to_json(self, *, indent: int = 2, sort_keys: bool = True) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=sort_keys)


@dataclass(frozen=True)
class WorkflowStageStatusSummary(_SerializableRecord):
    """Compact status view for one workflow stage."""

    stage_name: str
    label: str
    manifest_type: str
    input_manifest_path: Path
    status_path: Path
    status_exists: bool
    status: str | None = None
    attempt: int | None = None
    active_substage: str | None = None
    detail: str | None = None
    error_message: str | None = None
    stage_output_manifest_path: Path | None = None
    last_action: str | None = None


@dataclass(frozen=True)
class WorkflowSuggestedCommand(_SerializableRecord):
    """Suggested next command derived from the current workflow status."""

    kind: str
    summary: str
    command: str
    target_stage_name: str | None = None
    target_label: str | None = None


@dataclass(frozen=True)
class WorkflowStatusReport(_SerializableRecord):
    """Serializable workflow-status summary for one stage or one project chain."""

    source_manifest_path: Path
    manifest_type: str
    workflow_name: str
    project_root: Path | None = None
    stage_statuses: tuple[WorkflowStageStatusSummary, ...] = ()
    suggested_commands: tuple[WorkflowSuggestedCommand, ...] = ()


def collect_workflow_status(path: PathLike) -> WorkflowStatusReport:
    """Collect a compact workflow-status report for a stage or project manifest."""

    source_path = Path(path).expanduser()
    manifest = _load_workflow_manifest(source_path)
    if isinstance(manifest, ProjectWorkflowManifest):
        return _collect_project_workflow_status(source_path, manifest)
    stage_summary = _stage_status_summary(
        manifest=manifest,
        label=_default_stage_label(manifest),
        input_manifest_path=source_path,
    )
    return WorkflowStatusReport(
        source_manifest_path=source_path,
        manifest_type=_manifest_type(manifest),
        workflow_name=manifest.project.workflow_name,
        project_root=manifest.project.project_root,
        stage_statuses=(stage_summary,),
        suggested_commands=_stage_suggested_commands(source_path, stage_summary),
    )


def _collect_project_workflow_status(
    source_path: Path,
    manifest: ProjectWorkflowManifest,
) -> WorkflowStatusReport:
    stage_action_by_key: dict[tuple[str, str], str] = {}
    input_manifest_by_label: dict[str, Path] = {}
    if manifest.result is not None:
        stage_action_by_key = {
            (record.stage_name, record.label): record.action
            for record in manifest.result.stage_results
        }
        input_manifest_by_label = {
            record.label: record.path
            for record in manifest.result.resolved_input_manifests
        }

    stage_summaries = []
    for stage_plan in resolve_project_stage_manifests(manifest):
        input_manifest_path = input_manifest_by_label.get(
            stage_plan.input_manifest_label,
            stage_plan.input_manifest_path,
        )
        stage_summaries.append(
            _stage_status_summary(
                manifest=stage_plan.manifest,
                label=stage_plan.label,
                input_manifest_path=input_manifest_path,
                last_action=stage_action_by_key.get((stage_plan.stage_name, stage_plan.label)),
            )
        )

    return WorkflowStatusReport(
        source_manifest_path=source_path,
        manifest_type=PROJECT_MANIFEST_TYPE,
        workflow_name=manifest.project.workflow_name,
        project_root=manifest.project.project_root,
        stage_statuses=tuple(stage_summaries),
        suggested_commands=_project_suggested_commands(source_path, stage_summaries),
    )


def _stage_status_summary(
    *,
    manifest: StageWorkflowManifest,
    label: str,
    input_manifest_path: Path,
    last_action: str | None = None,
) -> WorkflowStageStatusSummary:
    status_path = stage_resume_state_path(manifest)
    if not status_path.exists():
        return WorkflowStageStatusSummary(
            stage_name=_stage_name(manifest),
            label=label,
            manifest_type=_manifest_type(manifest),
            input_manifest_path=input_manifest_path,
            status_path=status_path,
            status_exists=False,
            stage_output_manifest_path=_default_output_manifest_path(manifest),
            last_action=last_action,
        )

    state = read_stage_resume_state(status_path)
    return WorkflowStageStatusSummary(
        stage_name=state.stage_name,
        label=label,
        manifest_type=_manifest_type(manifest),
        input_manifest_path=input_manifest_path,
        status_path=status_path,
        status_exists=True,
        status=state.status,
        attempt=state.attempt,
        active_substage=state.active_substage,
        detail=state.detail,
        error_message=state.error_message,
        stage_output_manifest_path=state.stage_manifest_path or _default_output_manifest_path(manifest),
        last_action=last_action,
    )


def _load_workflow_manifest(path: Path) -> WorkflowManifest:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError(f"Workflow manifest at {path} must decode to a JSON object.")
    manifest_type = str(payload.get("manifest_type", "")).strip().lower()
    if manifest_type == PREPARE_MANIFEST_TYPE:
        return PrepareWorkflowManifest.from_dict(payload)
    if manifest_type == RUN_MANIFEST_TYPE:
        return RunWorkflowManifest.from_dict(payload)
    if manifest_type == ANALYZE_MANIFEST_TYPE:
        return AnalyzeWorkflowManifest.from_dict(payload)
    if manifest_type == PROJECT_MANIFEST_TYPE:
        return ProjectWorkflowManifest.from_dict(payload)
    raise ValueError(f"Unsupported workflow manifest_type={manifest_type!r} in {path}.")


def _stage_name(manifest: StageWorkflowManifest) -> str:
    if isinstance(manifest, PrepareWorkflowManifest):
        return "prepare"
    if isinstance(manifest, RunWorkflowManifest):
        return "run"
    if isinstance(manifest, AnalyzeWorkflowManifest):
        return "analyze"
    raise TypeError(f"Unsupported workflow manifest type: {type(manifest).__name__}.")


def _manifest_type(manifest: WorkflowManifest) -> str:
    if isinstance(manifest, PrepareWorkflowManifest):
        return PREPARE_MANIFEST_TYPE
    if isinstance(manifest, RunWorkflowManifest):
        return RUN_MANIFEST_TYPE
    if isinstance(manifest, AnalyzeWorkflowManifest):
        return ANALYZE_MANIFEST_TYPE
    if isinstance(manifest, ProjectWorkflowManifest):
        return PROJECT_MANIFEST_TYPE
    raise TypeError(f"Unsupported workflow manifest type: {type(manifest).__name__}.")


def _default_stage_label(manifest: StageWorkflowManifest) -> str:
    if isinstance(manifest, PrepareWorkflowManifest):
        return "prepare"
    if isinstance(manifest, RunWorkflowManifest):
        return str(manifest.spec.workflow_kind)
    return str(manifest.spec.analysis_kind)


def _default_output_manifest_path(manifest: StageWorkflowManifest) -> Path:
    if isinstance(manifest, PrepareWorkflowManifest):
        return layout_from_prepare_manifest(manifest).manifest_path()
    if isinstance(manifest, RunWorkflowManifest):
        return layout_from_run_manifest(manifest).manifest_path()
    if isinstance(manifest, AnalyzeWorkflowManifest):
        return layout_from_analyze_manifest(manifest).manifest_path()
    raise TypeError(f"Unsupported workflow manifest type: {type(manifest).__name__}.")


def _jsonify(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _jsonify(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_jsonify(item) for item in value]
    if isinstance(value, list):
        return [_jsonify(item) for item in value]
    return value


def _stage_suggested_commands(
    source_path: Path,
    summary: WorkflowStageStatusSummary,
) -> tuple[WorkflowSuggestedCommand, ...]:
    if summary.status == STAGE_STATUS_COMPLETED:
        return ()

    primary_command = _primary_stage_command(summary)
    return (
        WorkflowSuggestedCommand(
            kind="rerun_stage",
            summary=f"Rerun the {summary.stage_name} stage explicitly.",
            command=f"evbopenmm-workflow {primary_command} {source_path}",
            target_stage_name=summary.stage_name,
            target_label=summary.label,
        ),
        WorkflowSuggestedCommand(
            kind="resume_stage",
            summary="Use stage-local restart handling for this manifest.",
            command=f"evbopenmm-workflow resume {source_path}",
            target_stage_name=summary.stage_name,
            target_label=summary.label,
        ),
    )


def _project_suggested_commands(
    source_path: Path,
    stage_summaries: list[WorkflowStageStatusSummary],
) -> tuple[WorkflowSuggestedCommand, ...]:
    incomplete = [summary for summary in stage_summaries if summary.status != STAGE_STATUS_COMPLETED]
    if not incomplete:
        return ()

    suggestions: list[WorkflowSuggestedCommand] = [
        WorkflowSuggestedCommand(
            kind="recover_auto",
            summary="Recover the project using the current status view.",
            command=f"evbopenmm-workflow recover {source_path}",
        )
    ]

    first_incomplete = incomplete[0]
    if first_incomplete.stage_name in {"prepare", "run", "analyze"}:
        suggestions.append(
            WorkflowSuggestedCommand(
                kind="recover_from_stage",
                summary=f"Rerun from the first incomplete stage: {first_incomplete.stage_name}.",
                command=(
                    "evbopenmm-workflow recover --from-stage "
                    f"{first_incomplete.stage_name} {source_path}"
                ),
                target_stage_name=first_incomplete.stage_name,
                target_label=first_incomplete.label,
            )
        )

    if any(summary.stage_name == "analyze" for summary in incomplete):
        for summary in incomplete:
            if summary.stage_name != "analyze":
                continue
            suggestions.append(
                WorkflowSuggestedCommand(
                    kind="recover_analyze_label",
                    summary=f"Rerun only the analyze stage labeled {summary.label}.",
                    command=(
                        "evbopenmm-workflow recover --analyze-label "
                        f"{summary.label} {source_path}"
                    ),
                    target_stage_name=summary.stage_name,
                    target_label=summary.label,
                )
            )

    return tuple(suggestions)


def _primary_stage_command(summary: WorkflowStageStatusSummary) -> str:
    if summary.stage_name == "prepare":
        return "prepare"
    if summary.stage_name == "run":
        return "run"
    if summary.stage_name == "analyze":
        return "analyze"
    raise ValueError(f"Unsupported workflow stage name {summary.stage_name!r}.")


__all__ = [
    "WorkflowSuggestedCommand",
    "WorkflowStageStatusSummary",
    "WorkflowStatusReport",
    "collect_workflow_status",
]
