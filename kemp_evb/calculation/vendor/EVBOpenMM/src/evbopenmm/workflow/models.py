from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any, Mapping
import os


PathLike = os.PathLike[str] | str
WORKFLOW_MANIFEST_VERSION = "evb_workflow_v1"
PREPARE_MANIFEST_TYPE = "prepare"
RUN_MANIFEST_TYPE = "run"
ANALYZE_MANIFEST_TYPE = "analyze"
PROJECT_MANIFEST_TYPE = "project"


class _SerializableRecord:
    """Shared JSON helpers for workflow records."""

    def to_dict(self) -> dict[str, object]:
        return _jsonify(asdict(self))

    def to_json(self, *, indent: int = 2, sort_keys: bool = True) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=sort_keys)


class _WorkflowManifestRecord(_SerializableRecord):
    """Adds file-based JSON helpers for top-level workflow manifests."""

    def write_json(self, path: PathLike, *, indent: int = 2, sort_keys: bool = True) -> Path:
        output_path = Path(path).expanduser()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(self.to_json(indent=indent, sort_keys=sort_keys) + "\n", encoding="utf-8")
        return output_path


@dataclass(frozen=True)
class WorkflowProject(_SerializableRecord):
    """Project-level workflow metadata shared by prepare/run/analyze manifests."""

    workflow_name: str
    output_root: Path
    project_root: Path | None = None
    manifest_version: str = WORKFLOW_MANIFEST_VERSION
    description: str | None = None
    tags: tuple[str, ...] = ()

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "WorkflowProject":
        return cls(
            workflow_name=_require_str(payload, "workflow_name"),
            output_root=_optional_path(payload.get("output_root"), field_name="output_root", required=True),
            project_root=_optional_path(payload.get("project_root"), field_name="project_root"),
            manifest_version=_require_str(payload, "manifest_version"),
            description=_optional_str(payload.get("description"), field_name="description"),
            tags=_string_tuple(payload.get("tags"), field_name="tags"),
        )


@dataclass(frozen=True)
class StateInputRecord(_SerializableRecord):
    """Single named state input used by the prepare stage."""

    state_name: str
    pdb_path: Path

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "StateInputRecord":
        return cls(
            state_name=_require_str(payload, "state_name"),
            pdb_path=_optional_path(payload.get("pdb_path"), field_name="pdb_path", required=True),
        )


@dataclass(frozen=True)
class AtomReferenceRecord(_SerializableRecord):
    """Serializable EVB atom reference used by workflow manifests."""

    state_name: str
    chain_id: str
    residue_id: int
    atom_name: str

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "AtomReferenceRecord":
        return cls(
            state_name=_require_str(payload, "state_name"),
            chain_id=_require_str(payload, "chain_id"),
            residue_id=_require_int(payload, "residue_id"),
            atom_name=_require_str(payload, "atom_name"),
        )

    def to_atom_reference(self) -> tuple[str, str, int, str]:
        return (self.state_name, self.chain_id, self.residue_id, self.atom_name)


@dataclass(frozen=True)
class AtomMappingRecord(_SerializableRecord):
    """Serializable atom-identity mapping entry for EVB state preparation."""

    source: AtomReferenceRecord
    target: AtomReferenceRecord

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "AtomMappingRecord":
        return cls(
            source=AtomReferenceRecord.from_dict(_mapping_value(payload, "source")),
            target=AtomReferenceRecord.from_dict(_mapping_value(payload, "target")),
        )

    def to_mapping_pair(self) -> tuple[tuple[str, str, int, str], tuple[str, str, int, str]]:
        return self.source.to_atom_reference(), self.target.to_atom_reference()


@dataclass(frozen=True)
class LigandSdfInputRecord(_SerializableRecord):
    """Serializable SDF/MOL input descriptor for OpenFF parameterization."""

    residue_name: str
    path: Path
    index: int | None = None
    charge_method: str | None = None

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "LigandSdfInputRecord":
        return cls(
            residue_name=_require_str(payload, "residue_name"),
            path=_optional_path(payload.get("path"), field_name="path", required=True),
            index=_optional_int(payload.get("index"), field_name="index"),
            charge_method=_optional_str(payload.get("charge_method"), field_name="charge_method"),
        )

    def to_parameterizer_entry(self) -> dict[str, object]:
        payload: dict[str, object] = {"path": str(self.path)}
        if self.index is not None:
            payload["index"] = int(self.index)
        if self.charge_method is not None:
            payload["charge_method"] = self.charge_method
        return payload


@dataclass(frozen=True)
class StageFileRecord(_SerializableRecord):
    """Named workflow artifact path used across stage specs and results."""

    label: str
    path: Path

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "StageFileRecord":
        return cls(
            label=_require_str(payload, "label"),
            path=_optional_path(payload.get("path"), field_name="path", required=True),
        )


@dataclass(frozen=True)
class PreparedStateArtifact(_SerializableRecord):
    """Prepared artifact paths for one EVB state."""

    state_name: str
    topology_path: Path | None = None
    system_xml_path: Path | None = None
    positions_path: Path | None = None
    metadata_path: Path | None = None

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "PreparedStateArtifact":
        return cls(
            state_name=_require_str(payload, "state_name"),
            topology_path=_optional_path(payload.get("topology_path"), field_name="topology_path"),
            system_xml_path=_optional_path(payload.get("system_xml_path"), field_name="system_xml_path"),
            positions_path=_optional_path(payload.get("positions_path"), field_name="positions_path"),
            metadata_path=_optional_path(payload.get("metadata_path"), field_name="metadata_path"),
        )


@dataclass(frozen=True)
class LambdaScheduleStage(_SerializableRecord):
    """Serializable per-stage lambda schedule entry for runtime workflows."""

    state_lambdas: Mapping[str, float]
    label: str | None = None
    steps: int | None = None
    temperature_kelvin: float | None = None
    timestep_fs: float | None = None

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "LambdaScheduleStage":
        return cls(
            state_lambdas=_string_float_mapping(payload.get("state_lambdas"), field_name="state_lambdas"),
            label=_optional_str(payload.get("label"), field_name="label"),
            steps=_optional_int(payload.get("steps"), field_name="steps"),
            temperature_kelvin=_optional_float(payload.get("temperature_kelvin"), field_name="temperature_kelvin"),
            timestep_fs=_optional_float(payload.get("timestep_fs"), field_name="timestep_fs"),
        )


@dataclass(frozen=True)
class PrepareWorkflowSpec(_SerializableRecord):
    """Serializable inputs for a reproducible EVB prepare stage."""

    state_inputs: tuple[StateInputRecord, ...]
    output_directory: Path
    force_field: str = "amber14"
    atom_mapping: tuple[AtomMappingRecord, ...] = ()
    ligand_charges: Mapping[str, float] = field(default_factory=dict)
    ligand_sdf_inputs: tuple[LigandSdfInputRecord, ...] = ()
    ligand_smiles: Mapping[str, str] = field(default_factory=dict)
    smirnoff_forcefield: str = "openff-2.2.0.offxml"
    charge_method: str = "am1bcc"
    verbose: bool = False
    evb_initial_state: str | None = None
    auto_build_forces: bool = True
    reference_state: str | None = None
    force_field_files: tuple[Path, ...] = ()
    parameterizer: str | None = None
    parameterizer_inputs: tuple[StageFileRecord, ...] = ()
    water_model: str | None = None
    solvent_padding_nm: float | None = None
    box_size_nm: tuple[float, float, float] | None = None
    ionic_strength_molar: float | None = None
    neutralize: bool = False
    positive_ion: str = "Na+"
    negative_ion: str = "Cl-"
    num_added: int | None = None
    solvent_sphere_radius_angstrom: float | None = None
    solvent_residue_names: tuple[str, ...] = ()
    scaas_enabled: bool = False
    scaas_radius_angstrom: float | None = None
    scaas_polarisation: bool = False
    scaas_backend: str = "cpp"

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "PrepareWorkflowSpec":
        box_size_raw = payload.get("box_size_nm")
        box_size_nm = None if box_size_raw is None else _float_tuple(box_size_raw, field_name="box_size_nm")
        if box_size_nm is not None and len(box_size_nm) != 3:
            raise ValueError("'box_size_nm' must contain exactly three components.")
        return cls(
            state_inputs=_record_tuple(payload.get("state_inputs"), StateInputRecord, field_name="state_inputs"),
            output_directory=_optional_path(payload.get("output_directory"), field_name="output_directory", required=True),
            force_field=str(payload.get("force_field", "amber14")),
            atom_mapping=_record_tuple(payload.get("atom_mapping"), AtomMappingRecord, field_name="atom_mapping"),
            ligand_charges=_string_float_mapping(payload.get("ligand_charges"), field_name="ligand_charges"),
            ligand_sdf_inputs=_record_tuple(
                payload.get("ligand_sdf_inputs"),
                LigandSdfInputRecord,
                field_name="ligand_sdf_inputs",
            ),
            ligand_smiles=_string_mapping(payload.get("ligand_smiles"), field_name="ligand_smiles"),
            smirnoff_forcefield=str(payload.get("smirnoff_forcefield", "openff-2.2.0.offxml")),
            charge_method=str(payload.get("charge_method", "am1bcc")),
            verbose=bool(payload.get("verbose", False)),
            evb_initial_state=_optional_str(payload.get("evb_initial_state"), field_name="evb_initial_state"),
            auto_build_forces=bool(payload.get("auto_build_forces", True)),
            reference_state=_optional_str(payload.get("reference_state"), field_name="reference_state"),
            force_field_files=_path_tuple(payload.get("force_field_files"), field_name="force_field_files"),
            parameterizer=_optional_str(payload.get("parameterizer"), field_name="parameterizer"),
            parameterizer_inputs=_record_tuple(
                payload.get("parameterizer_inputs"),
                StageFileRecord,
                field_name="parameterizer_inputs",
            ),
            water_model=_optional_str(payload.get("water_model"), field_name="water_model"),
            solvent_padding_nm=_optional_float(payload.get("solvent_padding_nm"), field_name="solvent_padding_nm"),
            box_size_nm=box_size_nm,
            ionic_strength_molar=_optional_float(
                payload.get("ionic_strength_molar"),
                field_name="ionic_strength_molar",
            ),
            neutralize=bool(payload.get("neutralize", False)),
            positive_ion=str(payload.get("positive_ion", "Na+")),
            negative_ion=str(payload.get("negative_ion", "Cl-")),
            num_added=_optional_int(payload.get("num_added"), field_name="num_added"),
            solvent_sphere_radius_angstrom=_optional_float(
                payload.get("solvent_sphere_radius_angstrom"),
                field_name="solvent_sphere_radius_angstrom",
            ),
            solvent_residue_names=_string_tuple(payload.get("solvent_residue_names"), field_name="solvent_residue_names"),
            scaas_enabled=bool(payload.get("scaas_enabled", False)),
            scaas_radius_angstrom=_optional_float(
                payload.get("scaas_radius_angstrom"),
                field_name="scaas_radius_angstrom",
            ),
            scaas_polarisation=bool(payload.get("scaas_polarisation", False)),
            scaas_backend=str(payload.get("scaas_backend", "cpp")),
        )


@dataclass(frozen=True)
class PrepareWorkflowResult(_SerializableRecord):
    """Serializable outputs emitted by the prepare stage."""

    prepared_states: tuple[PreparedStateArtifact, ...] = ()
    atom_mapping: tuple[AtomMappingRecord, ...] = ()
    generated_files: tuple[StageFileRecord, ...] = ()
    manifest_path: Path | None = None

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "PrepareWorkflowResult":
        return cls(
            prepared_states=_record_tuple(
                payload.get("prepared_states"),
                PreparedStateArtifact,
                field_name="prepared_states",
            ),
            atom_mapping=_record_tuple(payload.get("atom_mapping"), AtomMappingRecord, field_name="atom_mapping"),
            generated_files=_record_tuple(payload.get("generated_files"), StageFileRecord, field_name="generated_files"),
            manifest_path=_optional_path(payload.get("manifest_path"), field_name="manifest_path"),
        )


@dataclass(frozen=True)
class RunWorkflowSpec(_SerializableRecord):
    """Serializable inputs for a reproducible EVB run stage."""

    workflow_kind: str
    prepared_manifest_path: Path
    output_directory: Path
    states: tuple[str, ...] = ()
    platform_name: str | None = None
    platform_properties: Mapping[str, str] = field(default_factory=dict)
    integrator_name: str | None = None
    random_seed: int | None = None
    initialize_velocities: str = "first"
    temperature_kelvin: float | None = None
    pressure_bar: float | None = None
    timestep_fs: float | None = None
    friction_per_ps: float | None = None
    equilibration_steps: int | None = None
    production_steps: int | None = None
    sample_interval: int | None = None
    coordinate_report_interval_ps: float | None = None
    coordinate_csv_interval_steps: int | None = None
    fail_on_instability: bool = True
    lambda_schedule: tuple[LambdaScheduleStage, ...] = ()

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "RunWorkflowSpec":
        return cls(
            workflow_kind=_require_str(payload, "workflow_kind"),
            prepared_manifest_path=_optional_path(
                payload.get("prepared_manifest_path"),
                field_name="prepared_manifest_path",
                required=True,
            ),
            output_directory=_optional_path(payload.get("output_directory"), field_name="output_directory", required=True),
            states=_string_tuple(payload.get("states"), field_name="states"),
            platform_name=_optional_str(payload.get("platform_name"), field_name="platform_name"),
            platform_properties=_string_mapping(payload.get("platform_properties"), field_name="platform_properties"),
            integrator_name=_optional_str(payload.get("integrator_name"), field_name="integrator_name"),
            random_seed=_optional_int(payload.get("random_seed"), field_name="random_seed"),
            initialize_velocities=str(payload.get("initialize_velocities", "first")),
            temperature_kelvin=_optional_float(payload.get("temperature_kelvin"), field_name="temperature_kelvin"),
            pressure_bar=_optional_float(payload.get("pressure_bar"), field_name="pressure_bar"),
            timestep_fs=_optional_float(payload.get("timestep_fs"), field_name="timestep_fs"),
            friction_per_ps=_optional_float(payload.get("friction_per_ps"), field_name="friction_per_ps"),
            equilibration_steps=_optional_int(payload.get("equilibration_steps"), field_name="equilibration_steps"),
            production_steps=_optional_int(payload.get("production_steps"), field_name="production_steps"),
            sample_interval=_optional_int(payload.get("sample_interval"), field_name="sample_interval"),
            coordinate_report_interval_ps=_optional_float(
                payload.get("coordinate_report_interval_ps"),
                field_name="coordinate_report_interval_ps",
            ),
            coordinate_csv_interval_steps=_optional_int(
                payload.get("coordinate_csv_interval_steps"),
                field_name="coordinate_csv_interval_steps",
            ),
            fail_on_instability=bool(payload.get("fail_on_instability", True)),
            lambda_schedule=_record_tuple(payload.get("lambda_schedule"), LambdaScheduleStage, field_name="lambda_schedule"),
        )


@dataclass(frozen=True)
class RunWorkflowResult(_SerializableRecord):
    """Serializable outputs emitted by the run stage."""

    coordinate_files: tuple[StageFileRecord, ...] = ()
    sample_record_files: tuple[StageFileRecord, ...] = ()
    summary_files: tuple[StageFileRecord, ...] = ()
    manifest_path: Path | None = None

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "RunWorkflowResult":
        return cls(
            coordinate_files=_record_tuple(payload.get("coordinate_files"), StageFileRecord, field_name="coordinate_files"),
            sample_record_files=_record_tuple(
                payload.get("sample_record_files"),
                StageFileRecord,
                field_name="sample_record_files",
            ),
            summary_files=_record_tuple(payload.get("summary_files"), StageFileRecord, field_name="summary_files"),
            manifest_path=_optional_path(payload.get("manifest_path"), field_name="manifest_path"),
        )


@dataclass(frozen=True)
class AnalyzeWorkflowSpec(_SerializableRecord):
    """Serializable inputs for a reproducible EVB analyze stage."""

    analysis_kind: str
    run_manifest_path: Path
    output_directory: Path
    input_record_files: tuple[StageFileRecord, ...] = ()
    reference_files: tuple[Path, ...] = ()
    temperature_kelvin: float | None = None
    skip_samples: int = 0
    n_gap_bins: int | None = None
    min_points_per_bin: int | None = None
    minimum_replicas_per_bin: int = 1
    alpha_kcal_mol: float | None = None
    hij_kcal_mol: float | None = None
    target_dg_dagger_kcal_mol: float | None = None
    target_dg0_kcal_mol: float | None = None
    weight_dg_dagger: float = 1.0
    weight_dg0: float = 1.0
    coupling_pattern: tuple[tuple[str, str], ...] = ()
    diabatic_offset_pattern: tuple[str, ...] | None = None
    reaction_coordinate_states: tuple[int, int] | None = None
    coupling_bounds_kcal_mol: tuple[float, float] = (-20.0, 20.0)
    offset_bounds_kcal_mol: tuple[float, float] = (-50.0, 50.0)
    scipy_method: str = "Powell"
    scipy_tolerance: float | None = None

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "AnalyzeWorkflowSpec":
        coupling_pattern_raw = payload.get("coupling_pattern") or ()
        coupling_pattern: list[tuple[str, str]] = []
        for idx, entry in enumerate(coupling_pattern_raw):
            if not isinstance(entry, (list, tuple)) or len(entry) != 2:
                raise ValueError(
                    f"'coupling_pattern'[{idx}] must be a [state_i, state_j] pair."
                )
            coupling_pattern.append((str(entry[0]), str(entry[1])))

        diabatic_offset_pattern_raw = payload.get("diabatic_offset_pattern")
        if diabatic_offset_pattern_raw is None:
            diabatic_offset_pattern: tuple[str, ...] | None = None
        else:
            if not isinstance(diabatic_offset_pattern_raw, (list, tuple)):
                raise ValueError("'diabatic_offset_pattern' must be a list of state names or null.")
            diabatic_offset_pattern = tuple(str(name) for name in diabatic_offset_pattern_raw)

        reaction_coordinate_states_raw = payload.get("reaction_coordinate_states")
        if reaction_coordinate_states_raw is None:
            reaction_coordinate_states: tuple[int, int] | None = None
        else:
            if not isinstance(reaction_coordinate_states_raw, (list, tuple)) or len(reaction_coordinate_states_raw) != 2:
                raise ValueError("'reaction_coordinate_states' must be a pair of integer state indices.")
            reaction_coordinate_states = (
                int(reaction_coordinate_states_raw[0]),
                int(reaction_coordinate_states_raw[1]),
            )

        def _bounds_pair(payload_key: str, default: tuple[float, float]) -> tuple[float, float]:
            raw = payload.get(payload_key)
            if raw is None:
                return default
            parsed = _float_tuple(raw, field_name=payload_key)
            if len(parsed) != 2:
                raise ValueError(f"'{payload_key}' must contain exactly two components.")
            return parsed

        return cls(
            analysis_kind=_require_str(payload, "analysis_kind"),
            run_manifest_path=_optional_path(payload.get("run_manifest_path"), field_name="run_manifest_path", required=True),
            output_directory=_optional_path(payload.get("output_directory"), field_name="output_directory", required=True),
            input_record_files=_record_tuple(
                payload.get("input_record_files"),
                StageFileRecord,
                field_name="input_record_files",
            ),
            reference_files=_path_tuple(payload.get("reference_files"), field_name="reference_files"),
            temperature_kelvin=_optional_float(payload.get("temperature_kelvin"), field_name="temperature_kelvin"),
            skip_samples=int(payload.get("skip_samples", 0)),
            n_gap_bins=_optional_int(payload.get("n_gap_bins"), field_name="n_gap_bins"),
            min_points_per_bin=_optional_int(payload.get("min_points_per_bin"), field_name="min_points_per_bin"),
            minimum_replicas_per_bin=int(payload.get("minimum_replicas_per_bin", 1)),
            alpha_kcal_mol=_optional_float(payload.get("alpha_kcal_mol"), field_name="alpha_kcal_mol"),
            hij_kcal_mol=_optional_float(payload.get("hij_kcal_mol"), field_name="hij_kcal_mol"),
            target_dg_dagger_kcal_mol=_optional_float(
                payload.get("target_dg_dagger_kcal_mol"),
                field_name="target_dg_dagger_kcal_mol",
            ),
            target_dg0_kcal_mol=_optional_float(payload.get("target_dg0_kcal_mol"), field_name="target_dg0_kcal_mol"),
            weight_dg_dagger=float(payload.get("weight_dg_dagger", 1.0)),
            weight_dg0=float(payload.get("weight_dg0", 1.0)),
            coupling_pattern=tuple(coupling_pattern),
            diabatic_offset_pattern=diabatic_offset_pattern,
            reaction_coordinate_states=reaction_coordinate_states,
            coupling_bounds_kcal_mol=_bounds_pair("coupling_bounds_kcal_mol", (-20.0, 20.0)),
            offset_bounds_kcal_mol=_bounds_pair("offset_bounds_kcal_mol", (-50.0, 50.0)),
            scipy_method=str(payload.get("scipy_method", "Powell")),
            scipy_tolerance=_optional_float(
                payload.get("scipy_tolerance"),
                field_name="scipy_tolerance",
            ),
        )


@dataclass(frozen=True)
class AnalyzeWorkflowResult(_SerializableRecord):
    """Serializable outputs emitted by the analyze stage."""

    table_files: tuple[StageFileRecord, ...] = ()
    plot_files: tuple[StageFileRecord, ...] = ()
    report_files: tuple[StageFileRecord, ...] = ()
    manifest_path: Path | None = None

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "AnalyzeWorkflowResult":
        return cls(
            table_files=_record_tuple(payload.get("table_files"), StageFileRecord, field_name="table_files"),
            plot_files=_record_tuple(payload.get("plot_files"), StageFileRecord, field_name="plot_files"),
            report_files=_record_tuple(payload.get("report_files"), StageFileRecord, field_name="report_files"),
            manifest_path=_optional_path(payload.get("manifest_path"), field_name="manifest_path"),
        )


@dataclass(frozen=True)
class ProjectWorkflowSpec(_SerializableRecord):
    """Serializable inputs for a chained prepare/run/analyze project workflow."""

    output_directory: Path
    prepare_manifest_path: Path | None = None
    run_manifest_path: Path | None = None
    analyze_manifest_files: tuple[StageFileRecord, ...] = ()

    def __post_init__(self) -> None:
        if (
            self.prepare_manifest_path is None
            and self.run_manifest_path is None
            and not self.analyze_manifest_files
        ):
            raise ValueError(
                "Project workflow spec must define at least one stage manifest path."
            )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "ProjectWorkflowSpec":
        return cls(
            output_directory=_optional_path(payload.get("output_directory"), field_name="output_directory", required=True),
            prepare_manifest_path=_optional_path(payload.get("prepare_manifest_path"), field_name="prepare_manifest_path"),
            run_manifest_path=_optional_path(payload.get("run_manifest_path"), field_name="run_manifest_path"),
            analyze_manifest_files=_record_tuple(
                payload.get("analyze_manifest_files"),
                StageFileRecord,
                field_name="analyze_manifest_files",
            ),
        )


@dataclass(frozen=True)
class ProjectStageExecutionRecord(_SerializableRecord):
    """Serializable execution outcome for one chained workflow stage."""

    stage_name: str
    label: str
    input_manifest_path: Path
    output_manifest_path: Path
    action: str

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "ProjectStageExecutionRecord":
        return cls(
            stage_name=_require_str(payload, "stage_name"),
            label=_require_str(payload, "label"),
            input_manifest_path=_optional_path(
                payload.get("input_manifest_path"),
                field_name="input_manifest_path",
                required=True,
            ),
            output_manifest_path=_optional_path(
                payload.get("output_manifest_path"),
                field_name="output_manifest_path",
                required=True,
            ),
            action=_require_str(payload, "action"),
        )


@dataclass(frozen=True)
class ProjectWorkflowResult(_SerializableRecord):
    """Serializable outputs emitted by the chained project workflow executor."""

    stage_results: tuple[ProjectStageExecutionRecord, ...] = ()
    resolved_input_manifests: tuple[StageFileRecord, ...] = ()
    manifest_path: Path | None = None

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "ProjectWorkflowResult":
        return cls(
            stage_results=_record_tuple(
                payload.get("stage_results"),
                ProjectStageExecutionRecord,
                field_name="stage_results",
            ),
            resolved_input_manifests=_record_tuple(
                payload.get("resolved_input_manifests"),
                StageFileRecord,
                field_name="resolved_input_manifests",
            ),
            manifest_path=_optional_path(payload.get("manifest_path"), field_name="manifest_path"),
        )


@dataclass(frozen=True)
class PrepareWorkflowManifest(_WorkflowManifestRecord):
    """Top-level prepare manifest."""

    project: WorkflowProject
    spec: PrepareWorkflowSpec
    result: PrepareWorkflowResult | None = None
    manifest_type: str = PREPARE_MANIFEST_TYPE

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "PrepareWorkflowManifest":
        _validate_manifest_type(payload, expected=PREPARE_MANIFEST_TYPE)
        return cls(
            project=WorkflowProject.from_dict(_mapping_value(payload, "project")),
            spec=PrepareWorkflowSpec.from_dict(_mapping_value(payload, "spec")),
            result=_optional_record(payload.get("result"), PrepareWorkflowResult, field_name="result"),
            manifest_type=_require_str(payload, "manifest_type"),
        )

    @classmethod
    def read_json(cls, path: PathLike) -> "PrepareWorkflowManifest":
        return cls.from_dict(_load_json_mapping(path))


@dataclass(frozen=True)
class RunWorkflowManifest(_WorkflowManifestRecord):
    """Top-level run manifest."""

    project: WorkflowProject
    spec: RunWorkflowSpec
    result: RunWorkflowResult | None = None
    manifest_type: str = RUN_MANIFEST_TYPE

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "RunWorkflowManifest":
        _validate_manifest_type(payload, expected=RUN_MANIFEST_TYPE)
        return cls(
            project=WorkflowProject.from_dict(_mapping_value(payload, "project")),
            spec=RunWorkflowSpec.from_dict(_mapping_value(payload, "spec")),
            result=_optional_record(payload.get("result"), RunWorkflowResult, field_name="result"),
            manifest_type=_require_str(payload, "manifest_type"),
        )

    @classmethod
    def read_json(cls, path: PathLike) -> "RunWorkflowManifest":
        return cls.from_dict(_load_json_mapping(path))


@dataclass(frozen=True)
class AnalyzeWorkflowManifest(_WorkflowManifestRecord):
    """Top-level analyze manifest."""

    project: WorkflowProject
    spec: AnalyzeWorkflowSpec
    result: AnalyzeWorkflowResult | None = None
    manifest_type: str = ANALYZE_MANIFEST_TYPE

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "AnalyzeWorkflowManifest":
        _validate_manifest_type(payload, expected=ANALYZE_MANIFEST_TYPE)
        return cls(
            project=WorkflowProject.from_dict(_mapping_value(payload, "project")),
            spec=AnalyzeWorkflowSpec.from_dict(_mapping_value(payload, "spec")),
            result=_optional_record(payload.get("result"), AnalyzeWorkflowResult, field_name="result"),
            manifest_type=_require_str(payload, "manifest_type"),
        )

    @classmethod
    def read_json(cls, path: PathLike) -> "AnalyzeWorkflowManifest":
        return cls.from_dict(_load_json_mapping(path))


@dataclass(frozen=True)
class ProjectWorkflowManifest(_WorkflowManifestRecord):
    """Top-level chained project workflow manifest."""

    project: WorkflowProject
    spec: ProjectWorkflowSpec
    result: ProjectWorkflowResult | None = None
    manifest_type: str = PROJECT_MANIFEST_TYPE

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "ProjectWorkflowManifest":
        _validate_manifest_type(payload, expected=PROJECT_MANIFEST_TYPE)
        return cls(
            project=WorkflowProject.from_dict(_mapping_value(payload, "project")),
            spec=ProjectWorkflowSpec.from_dict(_mapping_value(payload, "spec")),
            result=_optional_record(payload.get("result"), ProjectWorkflowResult, field_name="result"),
            manifest_type=_require_str(payload, "manifest_type"),
        )

    @classmethod
    def read_json(cls, path: PathLike) -> "ProjectWorkflowManifest":
        return cls.from_dict(_load_json_mapping(path))


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


def _load_json_mapping(path: PathLike) -> Mapping[str, object]:
    payload = json.loads(Path(path).expanduser().read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError("Manifest JSON must decode to a mapping.")
    return payload


def _validate_manifest_type(payload: Mapping[str, object], *, expected: str) -> None:
    actual = _require_str(payload, "manifest_type")
    if actual != expected:
        raise ValueError(f"Expected manifest_type={expected!r}, received {actual!r}.")


def _mapping_value(payload: Mapping[str, object], field_name: str) -> Mapping[str, object]:
    value = payload.get(field_name)
    if not isinstance(value, Mapping):
        raise ValueError(f"'{field_name}' must be a mapping.")
    return value


def _require_str(payload: Mapping[str, object], field_name: str) -> str:
    value = payload.get(field_name)
    if not isinstance(value, str) or not value:
        raise ValueError(f"'{field_name}' must be a non-empty string.")
    return value


def _require_int(payload: Mapping[str, object], field_name: str) -> int:
    value = payload.get(field_name)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"'{field_name}' must be an integer.")
    return int(value)


def _optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"'{field_name}' must be a string when provided.")
    return value


def _optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"'{field_name}' must be an integer when provided.")
    return int(value)


def _optional_float(value: object, *, field_name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"'{field_name}' must be numeric when provided.")
    return float(value)


def _optional_path(value: object, *, field_name: str, required: bool = False) -> Path | None:
    if value is None:
        if required:
            raise ValueError(f"'{field_name}' is required.")
        return None
    if not isinstance(value, (str, os.PathLike)):
        raise ValueError(f"'{field_name}' must be path-like.")
    return Path(value).expanduser()


def _string_tuple(value: object, *, field_name: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"'{field_name}' must be a sequence of strings.")
    result: list[str] = []
    for item in value:
        if not isinstance(item, str):
            raise ValueError(f"'{field_name}' must contain only strings.")
        result.append(item)
    return tuple(result)


def _float_tuple(value: object, *, field_name: str) -> tuple[float, ...]:
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"'{field_name}' must be a sequence of numeric values.")
    result: list[float] = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise ValueError(f"'{field_name}' must contain only numeric values.")
        result.append(float(item))
    return tuple(result)


def _path_tuple(value: object, *, field_name: str) -> tuple[Path, ...]:
    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"'{field_name}' must be a sequence of paths.")
    return tuple(_optional_path(item, field_name=field_name, required=True) for item in value)


def _string_mapping(value: object, *, field_name: str) -> dict[str, str]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"'{field_name}' must be a mapping of strings.")
    result: dict[str, str] = {}
    for key, item in value.items():
        if not isinstance(key, str) or not isinstance(item, str):
            raise ValueError(f"'{field_name}' must contain only string keys and values.")
        result[key] = item
    return result


def _string_float_mapping(value: object, *, field_name: str) -> dict[str, float]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"'{field_name}' must be a mapping of state names to floats.")
    result: dict[str, float] = {}
    for key, item in value.items():
        if not isinstance(key, str):
            raise ValueError(f"'{field_name}' must contain only string keys.")
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise ValueError(f"'{field_name}' must contain only numeric values.")
        result[key] = float(item)
    return result


def _record_tuple(value: object, record_type: type, *, field_name: str) -> tuple[Any, ...]:
    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"'{field_name}' must be a sequence.")
    result: list[Any] = []
    for item in value:
        if not isinstance(item, Mapping):
            raise ValueError(f"'{field_name}' must contain only mappings.")
        result.append(record_type.from_dict(item))
    return tuple(result)


def _optional_record(value: object, record_type: type, *, field_name: str) -> Any | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ValueError(f"'{field_name}' must be a mapping when provided.")
    return record_type.from_dict(value)


__all__ = [
    "WORKFLOW_MANIFEST_VERSION",
    "PREPARE_MANIFEST_TYPE",
    "RUN_MANIFEST_TYPE",
    "ANALYZE_MANIFEST_TYPE",
    "PROJECT_MANIFEST_TYPE",
    "WorkflowProject",
    "StateInputRecord",
    "AtomReferenceRecord",
    "AtomMappingRecord",
    "LigandSdfInputRecord",
    "StageFileRecord",
    "PreparedStateArtifact",
    "LambdaScheduleStage",
    "PrepareWorkflowSpec",
    "PrepareWorkflowResult",
    "RunWorkflowSpec",
    "RunWorkflowResult",
    "AnalyzeWorkflowSpec",
    "AnalyzeWorkflowResult",
    "ProjectWorkflowSpec",
    "ProjectStageExecutionRecord",
    "ProjectWorkflowResult",
    "PrepareWorkflowManifest",
    "RunWorkflowManifest",
    "AnalyzeWorkflowManifest",
    "ProjectWorkflowManifest",
]
