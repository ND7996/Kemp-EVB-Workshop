from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
from typing import Any

from .layout import (
    ANALYZE_STAGE_NAME,
    PREPARE_STAGE_NAME,
    RUN_STAGE_NAME,
    WorkflowStageLayout,
    layout_from_analyze_manifest,
    layout_from_prepare_manifest,
    layout_from_run_manifest,
)
from .models import AnalyzeWorkflowManifest, PrepareWorkflowManifest, RunWorkflowManifest, StageFileRecord


PathLike = os.PathLike[str] | str
WORKFLOW_STAGE_RESUME_TYPE = "stage_resume_state"
DEFAULT_STAGE_STATUS_LABEL = "stage_status"
STAGE_STATUS_PENDING = "pending"
STAGE_STATUS_RUNNING = "running"
STAGE_STATUS_INTERRUPTED = "interrupted"
STAGE_STATUS_FAILED = "failed"
STAGE_STATUS_COMPLETED = "completed"
VALID_STAGE_STATUSES = frozenset(
    {
        STAGE_STATUS_PENDING,
        STAGE_STATUS_RUNNING,
        STAGE_STATUS_INTERRUPTED,
        STAGE_STATUS_FAILED,
        STAGE_STATUS_COMPLETED,
    }
)
VALID_STAGE_NAMES = frozenset({PREPARE_STAGE_NAME, RUN_STAGE_NAME, ANALYZE_STAGE_NAME})

WorkflowManifest = PrepareWorkflowManifest | RunWorkflowManifest | AnalyzeWorkflowManifest
LayoutOrManifest = WorkflowStageLayout | WorkflowManifest


class _SerializableRecord:
    def to_dict(self) -> dict[str, object]:
        return _jsonify(asdict(self))

    def to_json(self, *, indent: int = 2, sort_keys: bool = True) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=sort_keys)


@dataclass(frozen=True)
class SubstageStatusRecord(_SerializableRecord):
    """Serializable progress record for one workflow substage."""

    name: str
    status: str
    detail: str | None = None
    restart_metadata_file: Path | None = None
    generated_files: tuple[StageFileRecord, ...] = ()

    def __post_init__(self) -> None:
        _validate_stage_status(self.status, field_name="status")

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "SubstageStatusRecord":
        return cls(
            name=_require_str(payload, "name"),
            status=_require_str(payload, "status"),
            detail=_optional_str(payload.get("detail"), field_name="detail"),
            restart_metadata_file=_optional_path(
                payload.get("restart_metadata_file"),
                field_name="restart_metadata_file",
            ),
            generated_files=_record_tuple(payload.get("generated_files"), StageFileRecord, field_name="generated_files"),
        )


@dataclass(frozen=True)
class StageResumeState(_SerializableRecord):
    """Serializable restart/resume snapshot for one workflow stage."""

    stage_name: str
    stage_root: Path
    status: str
    attempt: int = 1
    stage_manifest_path: Path | None = None
    active_substage: str | None = None
    substages: tuple[SubstageStatusRecord, ...] = ()
    restart_metadata_file: Path | None = None
    generated_files: tuple[StageFileRecord, ...] = ()
    detail: str | None = None
    error_message: str | None = None
    record_type: str = WORKFLOW_STAGE_RESUME_TYPE

    def __post_init__(self) -> None:
        if self.stage_name not in VALID_STAGE_NAMES:
            raise ValueError(f"Unsupported stage name {self.stage_name!r}.")
        _validate_stage_status(self.status, field_name="status")
        if int(self.attempt) < 1:
            raise ValueError("'attempt' must be >= 1.")
        if self.record_type != WORKFLOW_STAGE_RESUME_TYPE:
            raise ValueError(
                f"Expected record_type={WORKFLOW_STAGE_RESUME_TYPE!r}, received {self.record_type!r}."
            )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "StageResumeState":
        record_type = _require_str(payload, "record_type")
        if record_type != WORKFLOW_STAGE_RESUME_TYPE:
            raise ValueError(
                f"Expected record_type={WORKFLOW_STAGE_RESUME_TYPE!r}, received {record_type!r}."
            )
        return cls(
            stage_name=_require_str(payload, "stage_name"),
            stage_root=_optional_path(payload.get("stage_root"), field_name="stage_root", required=True),
            status=_require_str(payload, "status"),
            attempt=int(payload.get("attempt", 1)),
            stage_manifest_path=_optional_path(payload.get("stage_manifest_path"), field_name="stage_manifest_path"),
            active_substage=_optional_str(payload.get("active_substage"), field_name="active_substage"),
            substages=_record_tuple(payload.get("substages"), SubstageStatusRecord, field_name="substages"),
            restart_metadata_file=_optional_path(
                payload.get("restart_metadata_file"),
                field_name="restart_metadata_file",
            ),
            generated_files=_record_tuple(payload.get("generated_files"), StageFileRecord, field_name="generated_files"),
            detail=_optional_str(payload.get("detail"), field_name="detail"),
            error_message=_optional_str(payload.get("error_message"), field_name="error_message"),
            record_type=record_type,
        )

    def write_json(self, path: PathLike, *, indent: int = 2, sort_keys: bool = True) -> Path:
        output_path = Path(path).expanduser()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(self.to_json(indent=indent, sort_keys=sort_keys) + "\n", encoding="utf-8")
        return output_path

    @classmethod
    def read_json(cls, path: PathLike) -> "StageResumeState":
        payload = json.loads(Path(path).expanduser().read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping):
            raise ValueError("Stage resume JSON must decode to a mapping.")
        return cls.from_dict(payload)


def stage_resume_state_from_manifest(
    manifest: WorkflowManifest,
    *,
    status: str,
    attempt: int = 1,
    active_substage: str | None = None,
    substages: Sequence[SubstageStatusRecord] = (),
    restart_metadata_file: Path | None = None,
    generated_files: Sequence[StageFileRecord] = (),
    detail: str | None = None,
    error_message: str | None = None,
) -> StageResumeState:
    """Build a stage-resume record anchored to an existing workflow manifest/layout."""

    layout = _layout_from_manifest(manifest)
    manifest_path = _manifest_path_for_resume(manifest, layout)
    return StageResumeState(
        stage_name=layout.stage_name,
        stage_root=layout.root,
        status=status,
        attempt=attempt,
        stage_manifest_path=manifest_path,
        active_substage=active_substage,
        substages=tuple(substages),
        restart_metadata_file=restart_metadata_file,
        generated_files=tuple(generated_files),
        detail=detail,
        error_message=error_message,
    )


def stage_resume_state_path(target: LayoutOrManifest | PathLike, *, label: str = DEFAULT_STAGE_STATUS_LABEL) -> Path:
    """Resolve the canonical path for a stage-resume record."""

    if isinstance(target, (str, os.PathLike)):
        return Path(target).expanduser()
    layout = _coerce_layout(target)
    return layout.restart_path(label)


def write_stage_resume_state(
    target: LayoutOrManifest | PathLike,
    state: StageResumeState,
    *,
    label: str = DEFAULT_STAGE_STATUS_LABEL,
    indent: int = 2,
    sort_keys: bool = True,
) -> Path:
    """Write a stage-resume record to the canonical restart path or an explicit path."""

    return state.write_json(stage_resume_state_path(target, label=label), indent=indent, sort_keys=sort_keys)


def read_stage_resume_state(
    target: LayoutOrManifest | PathLike,
    *,
    label: str = DEFAULT_STAGE_STATUS_LABEL,
) -> StageResumeState:
    """Read a stage-resume record from the canonical restart path or an explicit path."""

    return StageResumeState.read_json(stage_resume_state_path(target, label=label))


def _coerce_layout(target: LayoutOrManifest) -> WorkflowStageLayout:
    if isinstance(target, WorkflowStageLayout):
        return target
    return _layout_from_manifest(target)


def _layout_from_manifest(manifest: WorkflowManifest) -> WorkflowStageLayout:
    if isinstance(manifest, PrepareWorkflowManifest):
        return layout_from_prepare_manifest(manifest)
    if isinstance(manifest, RunWorkflowManifest):
        return layout_from_run_manifest(manifest)
    if isinstance(manifest, AnalyzeWorkflowManifest):
        return layout_from_analyze_manifest(manifest)
    raise TypeError(f"Unsupported workflow manifest type: {type(manifest).__name__}.")


def _manifest_path_for_resume(manifest: WorkflowManifest, layout: WorkflowStageLayout) -> Path:
    if isinstance(manifest, PrepareWorkflowManifest) and manifest.result is not None and manifest.result.manifest_path is not None:
        return manifest.result.manifest_path
    if isinstance(manifest, RunWorkflowManifest) and manifest.result is not None and manifest.result.manifest_path is not None:
        return manifest.result.manifest_path
    if isinstance(manifest, AnalyzeWorkflowManifest) and manifest.result is not None and manifest.result.manifest_path is not None:
        return manifest.result.manifest_path
    return layout.manifest_path()


def _validate_stage_status(value: object, *, field_name: str) -> None:
    normalized = str(value)
    if normalized not in VALID_STAGE_STATUSES:
        raise ValueError(
            f"Unsupported {field_name!r} value {normalized!r}. "
            f"Expected one of {tuple(sorted(VALID_STAGE_STATUSES))!r}."
        )


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


def _require_str(payload: Mapping[str, object], field_name: str) -> str:
    value = payload.get(field_name)
    if not isinstance(value, str):
        raise ValueError(f"'{field_name}' must be a string.")
    return value


def _optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"'{field_name}' must be a string or null.")
    return value


def _optional_path(value: object, *, field_name: str, required: bool = False) -> Path | None:
    if value is None:
        if required:
            raise ValueError(f"'{field_name}' is required.")
        return None
    if not isinstance(value, (str, os.PathLike)):
        raise ValueError(f"'{field_name}' must be a filesystem path.")
    return Path(value).expanduser()


def _record_tuple(value: object, record_type, *, field_name: str) -> tuple:
    if value is None:
        return tuple()
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise ValueError(f"'{field_name}' must be a sequence.")
    records: list[object] = []
    for item in value:
        if not isinstance(item, Mapping):
            raise ValueError(f"Entries in '{field_name}' must be mappings.")
        records.append(record_type.from_dict(item))
    return tuple(records)


__all__ = [
    "WORKFLOW_STAGE_RESUME_TYPE",
    "DEFAULT_STAGE_STATUS_LABEL",
    "STAGE_STATUS_PENDING",
    "STAGE_STATUS_RUNNING",
    "STAGE_STATUS_INTERRUPTED",
    "STAGE_STATUS_FAILED",
    "STAGE_STATUS_COMPLETED",
    "VALID_STAGE_STATUSES",
    "SubstageStatusRecord",
    "StageResumeState",
    "stage_resume_state_from_manifest",
    "stage_resume_state_path",
    "write_stage_resume_state",
    "read_stage_resume_state",
]
