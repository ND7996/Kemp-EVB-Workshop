from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from .execute import execute_manifest, load_workflow_manifest
from .layout import layout_from_analyze_manifest, layout_from_prepare_manifest, layout_from_run_manifest
from .models import AnalyzeWorkflowManifest, PrepareWorkflowManifest, RunWorkflowManifest
from .resume import (
    STAGE_STATUS_COMPLETED,
    STAGE_STATUS_FAILED,
    STAGE_STATUS_INTERRUPTED,
    STAGE_STATUS_RUNNING,
    StageResumeState,
    read_stage_resume_state,
    stage_resume_state_path,
)

if TYPE_CHECKING:
    from .. import EVBSystem


PathLike = os.PathLike[str] | str
WorkflowManifest = PrepareWorkflowManifest | RunWorkflowManifest | AnalyzeWorkflowManifest
WorkflowOutputManifest = WorkflowManifest
RESUME_ACTION_USE_EXISTING = "use_existing"
RESUME_ACTION_RERUN = "rerun"
RESUME_ACTION_RESUME = "resume"
VALID_RESUME_ACTIONS = frozenset(
    {
        RESUME_ACTION_USE_EXISTING,
        RESUME_ACTION_RERUN,
        RESUME_ACTION_RESUME,
    }
)


def execute_manifest_with_resume(
    path: PathLike,
    *,
    force_rerun: bool = False,
    evb_system_factory: type["EVBSystem"] | object = None,
) -> tuple[WorkflowOutputManifest, str]:
    """Execute a workflow manifest with simple restart-aware stage-status handling."""

    manifest = load_workflow_manifest(path)
    if force_rerun:
        return _execute_and_mark(path, evb_system_factory=evb_system_factory), RESUME_ACTION_RERUN

    resume_state = _read_existing_resume_state(manifest)
    existing_output_manifest = _read_existing_output_manifest(manifest)

    if resume_state is None:
        return _execute_and_mark(path, evb_system_factory=evb_system_factory), RESUME_ACTION_RERUN

    if resume_state.status == STAGE_STATUS_COMPLETED and existing_output_manifest is not None:
        return existing_output_manifest, RESUME_ACTION_USE_EXISTING

    if resume_state.status in {STAGE_STATUS_FAILED, STAGE_STATUS_INTERRUPTED, STAGE_STATUS_RUNNING}:
        return _execute_and_mark(path, evb_system_factory=evb_system_factory), RESUME_ACTION_RESUME

    return _execute_and_mark(path, evb_system_factory=evb_system_factory), RESUME_ACTION_RERUN


def _execute_and_mark(
    path: PathLike,
    *,
    evb_system_factory: type["EVBSystem"] | object = None,
) -> WorkflowOutputManifest:
    return execute_manifest(path, evb_system_factory=evb_system_factory)


def _read_existing_resume_state(manifest: WorkflowManifest) -> StageResumeState | None:
    resume_path = stage_resume_state_path(manifest)
    if not resume_path.exists():
        return None
    return read_stage_resume_state(resume_path)


def _read_existing_output_manifest(manifest: WorkflowManifest) -> WorkflowOutputManifest | None:
    candidate_paths = []
    if getattr(manifest, "result", None) is not None and getattr(manifest.result, "manifest_path", None) is not None:
        candidate_paths.append(manifest.result.manifest_path)
    candidate_paths.append(_default_output_manifest_path(manifest))

    for candidate_path in candidate_paths:
        if candidate_path is None:
            continue
        candidate = Path(candidate_path).expanduser()
        if not candidate.exists():
            continue
        return load_workflow_manifest(candidate)
    return None


def _default_output_manifest_path(manifest: WorkflowManifest) -> Path:
    if isinstance(manifest, PrepareWorkflowManifest):
        return layout_from_prepare_manifest(manifest).manifest_path()
    if isinstance(manifest, RunWorkflowManifest):
        return layout_from_run_manifest(manifest).manifest_path()
    if isinstance(manifest, AnalyzeWorkflowManifest):
        return layout_from_analyze_manifest(manifest).manifest_path()
    raise TypeError(f"Unsupported workflow manifest type: {type(manifest).__name__}.")


__all__ = [
    "RESUME_ACTION_USE_EXISTING",
    "RESUME_ACTION_RERUN",
    "RESUME_ACTION_RESUME",
    "VALID_RESUME_ACTIONS",
    "execute_manifest_with_resume",
]
