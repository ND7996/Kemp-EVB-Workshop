from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from .models import (
    AnalyzeWorkflowManifest,
    AnalyzeWorkflowSpec,
    PrepareWorkflowManifest,
    PrepareWorkflowSpec,
    RunWorkflowManifest,
    RunWorkflowSpec,
    WorkflowProject,
)


PREPARE_STAGE_NAME = "prepare"
RUN_STAGE_NAME = "run"
ANALYZE_STAGE_NAME = "analyze"


@dataclass(frozen=True)
class WorkflowStageLayout:
    """Canonical directory contract shared by all workflow stages."""

    stage_name: str
    root: Path
    manifests_dir: Path
    logs_dir: Path
    artifacts_dir: Path
    restart_dir: Path

    def directories(self) -> tuple[Path, ...]:
        return _unique_paths(
            (
                self.root,
                self.manifests_dir,
                self.logs_dir,
                self.artifacts_dir,
                self.restart_dir,
            )
        )

    def manifest_path(self, filename: str | None = None) -> Path:
        resolved_name = filename or f"{self.stage_name}_manifest.json"
        return self.manifests_dir / _artifact_filename(resolved_name, suffix="")

    def log_path(self, label: str = "stage") -> Path:
        return self.logs_dir / _artifact_filename(label, ".log")

    def restart_path(self, label: str = "stage") -> Path:
        return self.restart_dir / _artifact_filename(label, ".json")

    def artifact_path(self, label: str, suffix: str) -> Path:
        return self.artifacts_dir / _artifact_filename(label, suffix)


@dataclass(frozen=True)
class PrepareWorkflowLayout(WorkflowStageLayout):
    """Canonical directory contract for the prepare stage."""

    state_artifacts_dir: Path
    shared_artifacts_dir: Path

    def directories(self) -> tuple[Path, ...]:
        return _unique_paths(super().directories() + (self.state_artifacts_dir, self.shared_artifacts_dir))

    def state_directory(self, state_name: str) -> Path:
        return self.state_artifacts_dir / _slug(state_name)

    def state_artifact_path(self, state_name: str, label: str, suffix: str) -> Path:
        return self.state_directory(state_name) / _artifact_filename(label, suffix)

    def shared_artifact_path(self, label: str, suffix: str) -> Path:
        return self.shared_artifacts_dir / _artifact_filename(label, suffix)


@dataclass(frozen=True)
class RunWorkflowLayout(WorkflowStageLayout):
    """Canonical directory contract for the run stage."""

    coordinate_artifacts_dir: Path
    record_artifacts_dir: Path
    summary_artifacts_dir: Path

    def directories(self) -> tuple[Path, ...]:
        return _unique_paths(
            super().directories()
            + (
                self.coordinate_artifacts_dir,
                self.record_artifacts_dir,
                self.summary_artifacts_dir,
            )
        )

    def coordinate_path(self, label: str, suffix: str = ".dcd") -> Path:
        return self.coordinate_artifacts_dir / _artifact_filename(label, suffix)

    def record_path(self, label: str, suffix: str = ".json") -> Path:
        return self.record_artifacts_dir / _artifact_filename(label, suffix)

    def summary_path(self, label: str, suffix: str = ".json") -> Path:
        return self.summary_artifacts_dir / _artifact_filename(label, suffix)


@dataclass(frozen=True)
class AnalyzeWorkflowLayout(WorkflowStageLayout):
    """Canonical directory contract for the analyze stage."""

    table_artifacts_dir: Path
    plot_artifacts_dir: Path
    report_artifacts_dir: Path

    def directories(self) -> tuple[Path, ...]:
        return _unique_paths(
            super().directories()
            + (
                self.table_artifacts_dir,
                self.plot_artifacts_dir,
                self.report_artifacts_dir,
            )
        )

    def table_path(self, label: str, suffix: str = ".csv") -> Path:
        return self.table_artifacts_dir / _artifact_filename(label, suffix)

    def plot_path(self, label: str, suffix: str = ".png") -> Path:
        return self.plot_artifacts_dir / _artifact_filename(label, suffix)

    def report_path(self, label: str, suffix: str = ".json") -> Path:
        return self.report_artifacts_dir / _artifact_filename(label, suffix)


def resolve_project_output_root(project: WorkflowProject) -> Path:
    root = project.output_root.expanduser()
    if root.is_absolute():
        return root
    if project.project_root is None:
        return root
    return project.project_root.expanduser() / root


def resolve_prepare_layout(project: WorkflowProject, spec: PrepareWorkflowSpec) -> PrepareWorkflowLayout:
    root = _resolve_stage_root(project, spec.output_directory)
    artifacts_dir = root / "artifacts"
    return PrepareWorkflowLayout(
        stage_name=PREPARE_STAGE_NAME,
        root=root,
        manifests_dir=root / "manifests",
        logs_dir=root / "logs",
        artifacts_dir=artifacts_dir,
        restart_dir=root / "restart",
        state_artifacts_dir=artifacts_dir / "states",
        shared_artifacts_dir=artifacts_dir / "shared",
    )


def resolve_run_layout(project: WorkflowProject, spec: RunWorkflowSpec) -> RunWorkflowLayout:
    root = _resolve_stage_root(project, spec.output_directory)
    artifacts_dir = root / "artifacts"
    return RunWorkflowLayout(
        stage_name=RUN_STAGE_NAME,
        root=root,
        manifests_dir=root / "manifests",
        logs_dir=root / "logs",
        artifacts_dir=artifacts_dir,
        restart_dir=root / "restart",
        coordinate_artifacts_dir=artifacts_dir / "coordinates",
        record_artifacts_dir=artifacts_dir / "records",
        summary_artifacts_dir=artifacts_dir / "summaries",
    )


def resolve_analyze_layout(project: WorkflowProject, spec: AnalyzeWorkflowSpec) -> AnalyzeWorkflowLayout:
    root = _resolve_stage_root(project, spec.output_directory)
    artifacts_dir = root / "artifacts"
    return AnalyzeWorkflowLayout(
        stage_name=ANALYZE_STAGE_NAME,
        root=root,
        manifests_dir=root / "manifests",
        logs_dir=root / "logs",
        artifacts_dir=artifacts_dir,
        restart_dir=root / "restart",
        table_artifacts_dir=artifacts_dir / "tables",
        plot_artifacts_dir=artifacts_dir / "plots",
        report_artifacts_dir=artifacts_dir / "reports",
    )


def layout_from_prepare_manifest(manifest: PrepareWorkflowManifest) -> PrepareWorkflowLayout:
    return resolve_prepare_layout(manifest.project, manifest.spec)


def layout_from_run_manifest(manifest: RunWorkflowManifest) -> RunWorkflowLayout:
    return resolve_run_layout(manifest.project, manifest.spec)


def layout_from_analyze_manifest(manifest: AnalyzeWorkflowManifest) -> AnalyzeWorkflowLayout:
    return resolve_analyze_layout(manifest.project, manifest.spec)


def ensure_stage_layout(layout: WorkflowStageLayout) -> WorkflowStageLayout:
    for directory in layout.directories():
        directory.mkdir(parents=True, exist_ok=True)
    return layout


def _resolve_stage_root(project: WorkflowProject, output_directory: Path) -> Path:
    candidate = output_directory.expanduser()
    if candidate.is_absolute():
        return candidate
    return resolve_project_output_root(project) / candidate


def _artifact_filename(label: str, suffix: str) -> str:
    if suffix and not suffix.startswith("."):
        suffix = f".{suffix}"
    if not suffix and "." in label and "/" not in label and "\\" not in label:
        return _slug(label, preserve_dots=True)
    return f"{_slug(label)}{suffix}"


def _slug(value: str, *, preserve_dots: bool = False) -> str:
    pattern = r"[^A-Za-z0-9._-]+" if preserve_dots else r"[^A-Za-z0-9_-]+"
    normalized = re.sub(pattern, "_", value.strip())
    normalized = normalized.strip("._-")
    return normalized or "artifact"


def _unique_paths(paths: tuple[Path, ...]) -> tuple[Path, ...]:
    ordered: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        ordered.append(path)
    return tuple(ordered)


__all__ = [
    "PREPARE_STAGE_NAME",
    "RUN_STAGE_NAME",
    "ANALYZE_STAGE_NAME",
    "WorkflowStageLayout",
    "PrepareWorkflowLayout",
    "RunWorkflowLayout",
    "AnalyzeWorkflowLayout",
    "resolve_project_output_root",
    "resolve_prepare_layout",
    "resolve_run_layout",
    "resolve_analyze_layout",
    "layout_from_prepare_manifest",
    "layout_from_run_manifest",
    "layout_from_analyze_manifest",
    "ensure_stage_layout",
]
