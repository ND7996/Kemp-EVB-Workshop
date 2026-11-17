from __future__ import annotations

import json
import os
from pathlib import Path
from typing import TYPE_CHECKING

from .analyze import analyze_workflow
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
from .prepare import run_prepare_workflow
from .run import run_workflow

if TYPE_CHECKING:
    from collections.abc import Mapping

    from .. import EVBSystem


PathLike = os.PathLike[str] | str
WorkflowManifest = PrepareWorkflowManifest | RunWorkflowManifest | AnalyzeWorkflowManifest | ProjectWorkflowManifest
WorkflowOutputManifest = WorkflowManifest


def load_workflow_manifest(path: PathLike) -> WorkflowManifest:
    """Load a typed workflow manifest from JSON by inspecting its manifest_type."""

    manifest_path = Path(path).expanduser()
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Workflow manifest at {manifest_path} must decode to a JSON object.")

    manifest_type = str(payload.get("manifest_type", "")).strip().lower()
    if manifest_type == PREPARE_MANIFEST_TYPE:
        return PrepareWorkflowManifest.from_dict(payload)
    if manifest_type == RUN_MANIFEST_TYPE:
        return RunWorkflowManifest.from_dict(payload)
    if manifest_type == ANALYZE_MANIFEST_TYPE:
        return AnalyzeWorkflowManifest.from_dict(payload)
    if manifest_type == PROJECT_MANIFEST_TYPE:
        return ProjectWorkflowManifest.from_dict(payload)
    raise ValueError(
        f"Unsupported workflow manifest_type={manifest_type!r} in {manifest_path}. "
        f"Expected one of {PREPARE_MANIFEST_TYPE!r}, {RUN_MANIFEST_TYPE!r}, "
        f"{ANALYZE_MANIFEST_TYPE!r}, or {PROJECT_MANIFEST_TYPE!r}."
    )


def execute_prepare(path: PathLike, *, evb_system_factory: type["EVBSystem"] | object = None) -> PrepareWorkflowManifest:
    """Load and execute a prepare manifest from disk."""

    manifest = _expect_manifest_type(load_workflow_manifest(path), PrepareWorkflowManifest)
    kwargs = {}
    if evb_system_factory is not None:
        kwargs["evb_system_factory"] = evb_system_factory
    _, output_manifest = run_prepare_workflow(manifest, **kwargs)
    return output_manifest


def execute_run(path: PathLike, *, evb_system_factory: type["EVBSystem"] | object = None) -> RunWorkflowManifest:
    """Load and execute a run manifest from disk."""

    manifest = _expect_manifest_type(load_workflow_manifest(path), RunWorkflowManifest)
    kwargs = {}
    if evb_system_factory is not None:
        kwargs["evb_system_factory"] = evb_system_factory
    _, output_manifest = run_workflow(manifest, **kwargs)
    return output_manifest


def execute_analyze(path: PathLike) -> AnalyzeWorkflowManifest:
    """Load and execute an analyze manifest from disk."""

    manifest = _expect_manifest_type(load_workflow_manifest(path), AnalyzeWorkflowManifest)
    _, output_manifest = analyze_workflow(manifest)
    return output_manifest


def execute_project(
    path: PathLike,
    *,
    force_rerun: bool = False,
    restart_policy: str = "cascade",
    start_stage: str | None = None,
    analyze_label: str | None = None,
    evb_system_factory: type["EVBSystem"] | object = None,
) -> ProjectWorkflowManifest:
    """Load and execute a project manifest from disk."""

    from .project import execute_project_manifest

    manifest = _expect_manifest_type(load_workflow_manifest(path), ProjectWorkflowManifest)
    manifest_path = _resolved_manifest_path(path, manifest)
    return execute_project_manifest(
        manifest_path,
        force_rerun=force_rerun,
        restart_policy=restart_policy,
        start_stage=start_stage,
        analyze_label=analyze_label,
        evb_system_factory=evb_system_factory,
    )


def execute_manifest(path: PathLike, *, evb_system_factory: type["EVBSystem"] | object = None) -> WorkflowOutputManifest:
    """Execute any supported workflow manifest from disk and return its output manifest."""

    manifest = load_workflow_manifest(path)
    if isinstance(manifest, PrepareWorkflowManifest):
        return execute_prepare(path, evb_system_factory=evb_system_factory)
    if isinstance(manifest, RunWorkflowManifest):
        return execute_run(path, evb_system_factory=evb_system_factory)
    if isinstance(manifest, AnalyzeWorkflowManifest):
        return execute_analyze(path)
    if isinstance(manifest, ProjectWorkflowManifest):
        return execute_project(path, evb_system_factory=evb_system_factory)
    raise TypeError(f"Unsupported workflow manifest type: {type(manifest).__name__}.")


def _expect_manifest_type(manifest: WorkflowManifest, manifest_type: type[WorkflowManifest]) -> WorkflowManifest:
    if not isinstance(manifest, manifest_type):
        raise TypeError(
            f"Expected {manifest_type.__name__}, received {type(manifest).__name__}."
        )
    return manifest


def _resolved_manifest_path(path: PathLike, manifest: WorkflowManifest) -> Path:
    if isinstance(path, (str, os.PathLike)):
        return Path(path).expanduser()
    raise TypeError(f"Unsupported manifest path type for {type(manifest).__name__}: {type(path).__name__}.")


__all__ = [
    "load_workflow_manifest",
    "execute_prepare",
    "execute_run",
    "execute_analyze",
    "execute_project",
    "execute_manifest",
]
