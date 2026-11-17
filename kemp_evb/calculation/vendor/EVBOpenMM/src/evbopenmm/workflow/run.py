from __future__ import annotations

import csv
import hashlib
from collections.abc import Callable, Mapping, Sequence
import json
from pathlib import Path
from typing import Any

import openmm

from .. import EVBSystem
from ..core import chain_identifier
from ..input_reader.system_reader import read_pdb_state
from .layout import ensure_stage_layout, layout_from_prepare_manifest, layout_from_run_manifest
from .models import (
    LambdaScheduleStage,
    PrepareWorkflowManifest,
    RunWorkflowManifest,
    RunWorkflowResult,
    RunWorkflowSpec,
    StageFileRecord,
)
from .records import write_endpoint_md_comparison, write_fep_window_records
from .resume import (
    STAGE_STATUS_COMPLETED,
    STAGE_STATUS_FAILED,
    STAGE_STATUS_RUNNING,
    SubstageStatusRecord,
    stage_resume_state_from_manifest,
    write_stage_resume_state,
)


def build_evb_system_from_prepare_manifest(
    manifest: PrepareWorkflowManifest,
    *,
    evb_system_factory: Callable[..., EVBSystem] = EVBSystem,
) -> EVBSystem:
    """Reconstruct a prepared EVBSystem from persisted prepare-manifest artifacts."""

    if manifest.result is None:
        raise ValueError("Prepare manifest does not contain a result block.")
    if not manifest.result.prepared_states:
        raise ValueError("Prepare manifest does not list any prepared states.")

    prepare_layout = layout_from_prepare_manifest(manifest)
    prepared_state_inputs = {
        entry.state_name: entry.topology_path
        for entry in manifest.result.prepared_states
        if entry.topology_path is not None
    }
    if len(prepared_state_inputs) != len(manifest.result.prepared_states):
        raise ValueError("Prepare manifest is missing one or more prepared state topology paths.")

    prepared_state_system_xmls = _prepared_state_system_xml_inputs(manifest.result.prepared_states)
    _validate_prepared_state_metadata(manifest.result.prepared_states)
    _validate_prepare_parameter_bundle(manifest, prepare_layout)
    atom_mapping_records = manifest.result.atom_mapping or manifest.spec.atom_mapping
    atom_mapping_payload = (
        {source: target for source, target in (entry.to_mapping_pair() for entry in atom_mapping_records)}
        if atom_mapping_records
        else None
    )
    prepared_factory = getattr(evb_system_factory, "from_prepared_state_systems", None)
    if prepared_state_system_xmls and callable(prepared_factory):
        evb_system = prepared_factory(
            prepared_state_inputs,
            state_system_xmls=prepared_state_system_xmls,
            force_field=manifest.spec.force_field,
            atom_mapping=atom_mapping_payload,
            verbose=bool(manifest.spec.verbose),
            evb_initial_state=manifest.spec.evb_initial_state or manifest.spec.reference_state,
            auto_build_forces=True,
        )
    else:
        evb_system = evb_system_factory(
            prepared_state_inputs,
            force_field=manifest.spec.force_field,
            atom_mapping=atom_mapping_payload,
            ligand_charges=dict(manifest.spec.ligand_charges) if manifest.spec.ligand_charges else None,
            ligand_sdf_files=(
                {entry.residue_name: entry.to_parameterizer_entry() for entry in manifest.spec.ligand_sdf_inputs}
                if manifest.spec.ligand_sdf_inputs
                else None
            ),
            ligand_smiles=dict(manifest.spec.ligand_smiles) if manifest.spec.ligand_smiles else None,
            smirnoff_forcefield=manifest.spec.smirnoff_forcefield,
            charge_method=manifest.spec.charge_method,
            parameterization_output=prepare_layout.shared_artifacts_dir / "parameterization",
            verbose=bool(manifest.spec.verbose),
            evb_initial_state=manifest.spec.evb_initial_state or manifest.spec.reference_state,
            auto_build_forces=True,
        )
    if manifest.spec.scaas_enabled:
        if manifest.spec.scaas_radius_angstrom is None:
            raise ValueError("'scaas_radius_angstrom' is required when 'scaas_enabled=True'.")
        evb_system.enable_scaas(
            solvent_radius_angstrom=float(manifest.spec.scaas_radius_angstrom),
            polarisation=bool(manifest.spec.scaas_polarisation),
            backend=str(manifest.spec.scaas_backend),
            solvent_residue_names=manifest.spec.solvent_residue_names or None,
        )
    return evb_system


def _prepared_state_system_xml_inputs(prepared_states) -> dict[str, Path] | None:
    resolved: dict[str, Path] = {}
    missing_states: list[str] = []
    for entry in prepared_states:
        if entry.system_xml_path is None:
            missing_states.append(entry.state_name)
            continue
        resolved[entry.state_name] = entry.system_xml_path
    if not resolved:
        return None
    if missing_states:
        raise ValueError(
            "Prepare manifest contains only a partial set of serialized state-system XML artifacts. "
            "Missing -> " + ", ".join(sorted(missing_states))
        )
    return resolved


def _validate_prepared_state_metadata(prepared_states) -> None:
    for entry in prepared_states:
        if entry.topology_path is None:
            continue
        if not entry.topology_path.is_file():
            raise FileNotFoundError(f"Prepared topology PDB not found for '{entry.state_name}': {entry.topology_path}")
        if entry.system_xml_path is not None and not entry.system_xml_path.is_file():
            raise FileNotFoundError(
                f"Serialized state-system XML not found for '{entry.state_name}': {entry.system_xml_path}"
            )
        if entry.metadata_path is None:
            continue
        if not entry.metadata_path.is_file():
            raise FileNotFoundError(f"Prepared metadata JSON not found for '{entry.state_name}': {entry.metadata_path}")

        payload = _read_json_file(entry.metadata_path)
        if payload.get("state_name") is not None and str(payload["state_name"]) != entry.state_name:
            raise ValueError(
                f"Prepared metadata state-name mismatch for '{entry.state_name}': "
                f"metadata reports '{payload['state_name']}'."
            )

        state = read_pdb_state(entry.topology_path)
        if payload.get("n_atoms") is not None and int(payload["n_atoms"]) != state.n_atoms:
            raise ValueError(
                f"Prepared topology atom-count mismatch for '{entry.state_name}': "
                f"metadata reports {payload['n_atoms']}, PDB contains {state.n_atoms} atoms."
            )
        if payload.get("residue_names") is not None:
            metadata_names = tuple(sorted(str(name) for name in payload["residue_names"]))
            topology_names = tuple(sorted(state.residue_names))
            if metadata_names != topology_names:
                raise ValueError(
                    f"Prepared topology residue-name mismatch for '{entry.state_name}': "
                    f"metadata reports {metadata_names}, PDB contains {topology_names}."
                )
        signature = payload.get("topology_signature_sha256")
        if signature is not None and str(signature) != _topology_signature(state.topology):
            raise ValueError(
                f"Prepared topology signature mismatch for '{entry.state_name}'. "
                "The prepared PDB no longer matches the stored metadata."
            )


def _validate_prepare_parameter_bundle(manifest: PrepareWorkflowManifest, prepare_layout) -> None:
    payload = _read_prepare_parameter_bundle(manifest, prepare_layout)
    if payload is None:
        return
    bundle_type = payload.get("bundle_type")
    if bundle_type != "evb_parameter_bundle":
        raise ValueError(f"Unsupported parameter bundle type '{bundle_type}'.")
    bundle_version = payload.get("bundle_version")
    if bundle_version != "evb_parameter_bundle_v1":
        raise ValueError(f"Unsupported parameter bundle version '{bundle_version}'.")

    expected_state_names = tuple(entry.state_name for entry in manifest.result.prepared_states)
    bundle_state_names = tuple(str(name) for name in payload.get("state_names", ()))
    if bundle_state_names and bundle_state_names != expected_state_names:
        raise ValueError(
            "Parameter bundle state names do not match the prepare manifest: "
            f"bundle={bundle_state_names}, manifest={expected_state_names}."
        )

    bundle_force_field = payload.get("force_field")
    if isinstance(bundle_force_field, Mapping):
        bundle_force_field_name = bundle_force_field.get("name")
        if bundle_force_field_name is not None and str(bundle_force_field_name) != manifest.spec.force_field:
            raise ValueError(
                "Parameter bundle force-field name does not match the prepare manifest: "
                f"bundle={bundle_force_field_name}, manifest={manifest.spec.force_field}."
            )
        bundle_create_system_options = bundle_force_field.get("create_system_options")
        if bundle_create_system_options is not None and not isinstance(bundle_create_system_options, Mapping):
            raise ValueError("Parameter bundle force-field create_system_options must be a JSON object.")

    artifacts_payload = payload.get("artifacts")
    bundle_prepared_states = None
    if isinstance(artifacts_payload, Mapping):
        bundle_prepared_states = artifacts_payload.get("prepared_states")
    if isinstance(bundle_prepared_states, Sequence):
        names = tuple(str(entry.get("state_name")) for entry in bundle_prepared_states if isinstance(entry, Mapping))
        if names and names != expected_state_names:
            raise ValueError(
                "Parameter bundle prepared-state records do not match the prepare manifest: "
                f"bundle={names}, manifest={expected_state_names}."
            )

    software_payload = payload.get("software")
    saved_openmm_version = software_payload.get("openmm_version") if isinstance(software_payload, Mapping) else None
    current_openmm_version = getattr(openmm, "__version__", None)
    if saved_openmm_version and current_openmm_version and not _openmm_versions_compatible(
        str(saved_openmm_version),
        str(current_openmm_version),
    ):
        raise ValueError(
            "Parameter bundle OpenMM version is incompatible with the current runtime: "
            f"bundle={saved_openmm_version}, runtime={current_openmm_version}."
        )


def _resolve_parameter_bundle_path(manifest: PrepareWorkflowManifest, prepare_layout) -> Path | None:
    if manifest.result is not None:
        for record in manifest.result.generated_files:
            if record.label == "parameter_bundle":
                return record.path
    candidate = prepare_layout.shared_artifact_path("parameter_bundle", ".json")
    if candidate.is_file():
        return candidate
    return None


def _read_prepare_parameter_bundle(
    manifest: PrepareWorkflowManifest,
    prepare_layout,
) -> Mapping[str, object] | None:
    bundle_path = _resolve_parameter_bundle_path(manifest, prepare_layout)
    if bundle_path is None:
        return None
    if not bundle_path.is_file():
        raise FileNotFoundError(f"Prepare parameter bundle JSON not found: {bundle_path}")
    return _read_json_file(bundle_path)


def _read_json_file(path: Path) -> Mapping[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Failed to parse JSON file: {path}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"JSON payload must be an object: {path}")
    return payload


def _topology_signature(topology) -> str:
    digest = hashlib.sha256()
    for atom in topology.atoms():
        residue = atom.residue
        residue_id = str((residue.id or "")).strip() or str(residue.index + 1)
        fields = (
            chain_identifier(residue.chain) or "",
            residue_id,
            (residue.name or "").strip(),
            (atom.name or "").strip(),
            atom.element.symbol if atom.element is not None else "",
        )
        digest.update("|".join(fields).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _openmm_versions_compatible(saved_version: str, current_version: str) -> bool:
    saved_parts = _parse_version_prefix(saved_version)
    current_parts = _parse_version_prefix(current_version)
    if saved_parts is None or current_parts is None:
        return saved_version == current_version
    return saved_parts[:2] == current_parts[:2]


def _parse_version_prefix(value: str) -> tuple[int, ...] | None:
    prefix = value.strip().split("+", 1)[0].split("-", 1)[0]
    parts: list[int] = []
    for token in prefix.split("."):
        digits = "".join(ch for ch in token if ch.isdigit())
        if not digits:
            break
        parts.append(int(digits))
    if not parts:
        return None
    return tuple(parts)


def run_workflow(
    manifest: RunWorkflowManifest,
    *,
    evb_system_factory: Callable[..., EVBSystem] = EVBSystem,
) -> tuple[EVBSystem, RunWorkflowManifest]:
    """Execute a Phase B run workflow from a persisted prepare manifest."""

    layout = ensure_stage_layout(layout_from_run_manifest(manifest))
    current_substage = "build_evb_system"
    completed_substages: list[SubstageStatusRecord] = []
    write_stage_resume_state(
        layout,
        stage_resume_state_from_manifest(
            manifest,
            status=STAGE_STATUS_RUNNING,
            active_substage=current_substage,
            substages=(SubstageStatusRecord(current_substage, STAGE_STATUS_RUNNING),),
            detail="Run workflow started.",
        ),
    )

    try:
        prepare_manifest = PrepareWorkflowManifest.read_json(manifest.spec.prepared_manifest_path)
        evb_system = build_evb_system_from_prepare_manifest(prepare_manifest, evb_system_factory=evb_system_factory)
        prepare_layout = layout_from_prepare_manifest(prepare_manifest)
        prepared_runtime = _prepared_runtime_payload(
            evb_system=evb_system,
            prepare_manifest=prepare_manifest,
            prepare_layout=prepare_layout,
        )
        completed_substages.append(SubstageStatusRecord(current_substage, STAGE_STATUS_COMPLETED))

        workflow_kind = str(manifest.spec.workflow_kind).strip().lower()
        current_substage = workflow_kind
        if workflow_kind == "equilibrate_schedule":
            _, coordinate_files, sample_record_files, summary_files = _run_equilibration_workflow(
                evb_system=evb_system,
                spec=manifest.spec,
                layout=layout,
                prepared_runtime=prepared_runtime,
            )
        elif workflow_kind == "endpoint_md_comparison":
            coordinate_files, sample_record_files, summary_files = _run_endpoint_comparison_workflow(
                evb_system=evb_system,
                spec=manifest.spec,
                layout=layout,
                prepared_runtime=prepared_runtime,
            )
        elif workflow_kind in {"mapping_fep_umbrella", "mapping_fep_evb_profile"}:
            coordinate_files, sample_record_files, summary_files = _run_mapping_sampling_workflow(
                evb_system=evb_system,
                spec=manifest.spec,
                layout=layout,
                workflow_kind=workflow_kind,
                prepared_runtime=prepared_runtime,
            )
        else:
            raise NotImplementedError(
                "Phase B run v1 currently supports workflow_kind='equilibrate_schedule', "
                "'endpoint_md_comparison', 'mapping_fep_umbrella', and 'mapping_fep_evb_profile'."
            )
        workflow_files = tuple(coordinate_files) + tuple(sample_record_files) + tuple(summary_files)
        completed_substages.append(
            SubstageStatusRecord(
                current_substage,
                STAGE_STATUS_COMPLETED,
                generated_files=workflow_files,
            )
        )

        current_substage = "write_run_outputs"
        result = RunWorkflowResult(
            coordinate_files=coordinate_files,
            sample_record_files=sample_record_files,
            summary_files=summary_files,
            manifest_path=layout.manifest_path(),
        )
        output_manifest = RunWorkflowManifest(
            project=manifest.project,
            spec=manifest.spec,
            result=result,
        )
        output_manifest.write_json(result.manifest_path)
        stage_files = _run_stage_output_files(result)
        completed_substages.append(
            SubstageStatusRecord(
                current_substage,
                STAGE_STATUS_COMPLETED,
                generated_files=stage_files,
            )
        )
        _write_completed_stage_resume_state(
            layout,
            output_manifest,
            substages=completed_substages,
            generated_files=stage_files,
        )
        return evb_system, output_manifest
    except Exception as exc:  # noqa: BLE001 - record failed stage resume state and re-raise
        _write_failed_stage_resume_state(
            layout,
            manifest,
            current_substage=current_substage,
            completed_substages=completed_substages,
            exc=exc,
        )
        raise


def _run_equilibration_workflow(
    *,
    evb_system: EVBSystem,
    spec: RunWorkflowSpec,
    layout,
    prepared_runtime: Mapping[str, object],
) -> tuple[dict[str, object], tuple[StageFileRecord, ...], tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    stage_schedule = _resolved_lambda_schedule(spec, evb_system.state_names)
    timesteps = tuple(
        stage.timestep_fs if stage.timestep_fs is not None else float(spec.timestep_fs or 1.0)
        for stage in stage_schedule
    )
    steps = tuple(
        int(stage.steps if stage.steps is not None else int(spec.equilibration_steps or 100))
        for stage in stage_schedule
    )
    temperatures = tuple(
        stage.temperature_kelvin if stage.temperature_kelvin is not None else float(spec.temperature_kelvin or 300.0)
        for stage in stage_schedule
    )
    lambda_schedule = tuple(dict(stage.state_lambdas) for stage in stage_schedule)

    coordinate_report_file = None
    coordinate_topology_file = None
    coordinate_files: list[StageFileRecord] = []
    if spec.coordinate_report_interval_ps is not None:
        coordinate_report_file = layout.coordinate_path("equilibration")
        coordinate_topology_file = layout.coordinate_path("equilibration_topology", ".pdb")
        coordinate_files.extend(
            (
                StageFileRecord("equilibration_dcd", coordinate_report_file),
                StageFileRecord("equilibration_topology", coordinate_topology_file),
            )
        )

    coordinate_csv_file = None
    sample_record_files: list[StageFileRecord] = []
    if spec.coordinate_csv_interval_steps is not None:
        coordinate_csv_file = layout.record_path("equilibration_coordinates", ".csv")
        sample_record_files.append(StageFileRecord("equilibration_coordinate_csv", coordinate_csv_file))

    summary_csv_path = layout.summary_path("equilibration_summary", ".csv")
    summary_json_path = layout.summary_path("equilibration_result", ".json")
    result_payload = dict(
        evb_system.equilibrate_schedule(
            timesteps=timesteps,
            steps=steps,
            temperatures=temperatures,
            lambda_schedule=lambda_schedule,
            friction=float(spec.friction_per_ps or 1.0),
            initialize_velocities=spec.initialize_velocities,
            random_seed=spec.random_seed,
            platform=_platform_from_spec(spec),
            platform_properties=spec.platform_properties or None,
            coordinate_report_interval_ps=spec.coordinate_report_interval_ps,
            coordinate_report_file=coordinate_report_file,
            coordinate_topology_file=coordinate_topology_file,
            summary_file=summary_csv_path,
            coordinate_csv_file=coordinate_csv_file,
            coordinate_csv_interval_steps=spec.coordinate_csv_interval_steps,
            fail_on_instability=bool(spec.fail_on_instability),
        )
    )
    result_payload["prepared_runtime"] = dict(prepared_runtime)
    _write_json(summary_json_path, result_payload)

    summary_files = (
        StageFileRecord("equilibration_summary_csv", summary_csv_path),
        StageFileRecord("equilibration_result_json", summary_json_path),
    )
    return result_payload, tuple(coordinate_files), tuple(sample_record_files), summary_files


def _run_endpoint_comparison_workflow(
    *,
    evb_system: EVBSystem,
    spec: RunWorkflowSpec,
    layout,
    prepared_runtime: Mapping[str, object],
) -> tuple[tuple[StageFileRecord, ...], tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    result_payload = dict(
        evb_system.run_endpoint_md_comparison(
            states=tuple(spec.states) or None,
            equilibration_steps=int(spec.equilibration_steps or 0),
            production_steps=int(spec.production_steps or 5000),
            sample_interval=int(spec.sample_interval or 10),
            temperature=float(spec.temperature_kelvin or 300.0),
            timestep=float(spec.timestep_fs or 1.0),
            friction=float(spec.friction_per_ps or 1.0),
            random_seed=spec.random_seed,
            platform=_platform_from_spec(spec),
            platform_properties=spec.platform_properties or None,
            initialize_velocities=_initialize_velocities_flag(spec.initialize_velocities),
        )
    )

    record_path = layout.record_path("endpoint_comparison_records", ".json")
    summary_json_path = layout.summary_path("endpoint_comparison_summary", ".json")
    summary_csv_path = layout.summary_path("endpoint_comparison_summary", ".csv")
    write_endpoint_md_comparison(record_path, result_payload)

    summary_payload = _endpoint_summary_payload(result_payload, prepared_runtime=prepared_runtime)
    _write_json(summary_json_path, summary_payload)
    _write_endpoint_summary_csv(summary_csv_path, summary_payload)

    sample_record_files = (
        StageFileRecord("endpoint_comparison_records_json", record_path),
    )
    summary_files = (
        StageFileRecord("endpoint_comparison_summary_json", summary_json_path),
        StageFileRecord("endpoint_comparison_summary_csv", summary_csv_path),
    )
    return tuple(), sample_record_files, summary_files


def _run_mapping_sampling_workflow(
    *,
    evb_system: EVBSystem,
    spec: RunWorkflowSpec,
    layout,
    workflow_kind: str,
    prepared_runtime: Mapping[str, object],
) -> tuple[tuple[StageFileRecord, ...], tuple[StageFileRecord, ...], tuple[StageFileRecord, ...]]:
    state_names = tuple(evb_system.state_names)
    if spec.states and tuple(spec.states) != state_names:
        raise ValueError(
            "Mapping workflows currently require the full EVB state set in canonical order. "
            f"Expected {state_names!r}, received {tuple(spec.states)!r}."
        )

    lambda_schedule = _resolved_mapping_lambda_schedule(spec, state_names)

    coordinate_report_file = None
    coordinate_topology_file = None
    coordinate_files: list[StageFileRecord] = []
    if spec.coordinate_report_interval_ps is not None:
        coordinate_report_file = layout.coordinate_path(f"{workflow_kind}_coordinates")
        coordinate_topology_file = layout.coordinate_path(f"{workflow_kind}_topology", ".pdb")
        coordinate_files.extend(
            (
                StageFileRecord(f"{workflow_kind}_dcd", coordinate_report_file),
                StageFileRecord(f"{workflow_kind}_topology", coordinate_topology_file),
            )
        )

    record_path = layout.record_path(f"{workflow_kind}_window_records", ".json")
    summary_path = layout.summary_path(f"{workflow_kind}_summary", ".json")
    window_records = tuple(
        evb_system.collect_mapping_fep_samples(
            lambda_schedule=lambda_schedule,
            equilibration_steps=int(spec.equilibration_steps or 0),
            production_steps=int(spec.production_steps or 5000),
            sample_interval=int(spec.sample_interval or 10),
            temperature=float(spec.temperature_kelvin or 300.0),
            timestep=float(spec.timestep_fs or 1.0),
            friction=float(spec.friction_per_ps or 1.0),
            random_seed=spec.random_seed,
            platform=_platform_from_spec(spec),
            platform_properties=spec.platform_properties or None,
            coordinate_report_interval_ps=spec.coordinate_report_interval_ps,
            coordinate_report_file=coordinate_report_file,
            coordinate_topology_file=coordinate_topology_file,
        )
    )
    write_fep_window_records(record_path, window_records)
    _write_json(
        summary_path,
        _mapping_summary_payload(
            workflow_kind=workflow_kind,
            state_names=state_names,
            lambda_schedule=lambda_schedule,
            window_records=window_records,
            spec=spec,
            prepared_runtime=prepared_runtime,
        ),
    )

    sample_record_files = (
        StageFileRecord(f"{workflow_kind}_window_records_json", record_path),
    )
    summary_files = (
        StageFileRecord(f"{workflow_kind}_summary_json", summary_path),
    )
    return tuple(coordinate_files), sample_record_files, summary_files


def _resolved_lambda_schedule(spec: RunWorkflowSpec, state_names: tuple[str, ...]):
    if spec.lambda_schedule:
        return spec.lambda_schedule
    default_value = 1.0 / float(len(state_names))
    return (
        LambdaScheduleStage(
            state_lambdas={state_name: default_value for state_name in state_names},
            steps=spec.equilibration_steps or 100,
            temperature_kelvin=spec.temperature_kelvin or 300.0,
            timestep_fs=spec.timestep_fs or 1.0,
        ),
    )


def _resolved_mapping_lambda_schedule(
    spec: RunWorkflowSpec,
    state_names: tuple[str, ...],
) -> tuple[dict[str, float], ...]:
    if not spec.lambda_schedule:
        raise ValueError("Mapping workflows require an explicit lambda_schedule with at least two windows.")

    windows: list[dict[str, float]] = []
    expected_states = set(state_names)
    for index, stage in enumerate(spec.lambda_schedule):
        provided_states = set(stage.state_lambdas)
        missing_states = sorted(expected_states - provided_states)
        extra_states = sorted(provided_states - expected_states)
        if missing_states or extra_states:
            raise ValueError(
                "Mapping lambda schedule must match the EVB state names exactly "
                f"(window_index={index}, missing={tuple(missing_states)!r}, extra={tuple(extra_states)!r})."
            )
        windows.append({state_name: float(stage.state_lambdas[state_name]) for state_name in state_names})

    if len(windows) < 2:
        raise ValueError("Mapping workflows require at least two lambda windows.")
    return tuple(windows)


def _platform_from_spec(spec: RunWorkflowSpec) -> openmm.Platform | None:
    if spec.platform_name is None:
        return None
    return openmm.Platform.getPlatformByName(spec.platform_name)


def _initialize_velocities_flag(value: object) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "always", "first"}:
        return True
    if normalized in {"0", "false", "no", "never", "none"}:
        return False
    raise ValueError(
        "Unsupported 'initialize_velocities' value for endpoint comparison. "
        "Use one of: 'first', 'always', 'never', 'true', or 'false'."
    )


def _mapping_summary_payload(
    *,
    workflow_kind: str,
    state_names: Sequence[str],
    lambda_schedule: Sequence[Mapping[str, float]],
    window_records: Sequence[object],
    spec: RunWorkflowSpec,
    prepared_runtime: Mapping[str, object],
) -> dict[str, object]:
    return {
        "workflow_kind": workflow_kind,
        "state_names": tuple(state_names),
        "lambda_schedule": tuple({state: float(value) for state, value in window.items()} for window in lambda_schedule),
        "window_count": len(window_records),
        "sample_counts": tuple(len(window.samples) for window in window_records),
        "equilibration_steps": int(spec.equilibration_steps or 0),
        "production_steps": int(spec.production_steps or 5000),
        "sample_interval": int(spec.sample_interval or 10),
        "temperature_kelvin": float(spec.temperature_kelvin or 300.0),
        "timestep_fs": float(spec.timestep_fs or 1.0),
        "friction_per_ps": float(spec.friction_per_ps or 1.0),
        "random_seed": None if spec.random_seed is None else int(spec.random_seed),
        "coordinate_report_interval_ps": None
        if spec.coordinate_report_interval_ps is None
        else float(spec.coordinate_report_interval_ps),
        "prepared_runtime": dict(prepared_runtime),
    }


def _endpoint_summary_payload(
    result_payload: Mapping[str, object],
    *,
    prepared_runtime: Mapping[str, object],
) -> dict[str, object]:
    state_names = tuple(str(value) for value in result_payload["state_names"])
    runs_raw = result_payload["runs"]
    return {
        "workflow_kind": "endpoint_md_comparison",
        "state_names": state_names,
        "energy_unit": str(result_payload.get("energy_unit", "kJ/mol")),
        "protocol": dict(result_payload["protocol"]),
        "initial_conditions": dict(result_payload["initial_conditions"]),
        "prepared_runtime": dict(prepared_runtime),
        "runs": tuple(
            {
                "state_name": state_name,
                "seed": None if run["seed"] is None else int(run["seed"]),
                "evb_sample_count": len(run["evb_endpoint"]["samples"]),
                "reference_sample_count": len(run["reference_ff"]["samples"]),
                "evb_mean_energies_kj_mol": dict(run["evb_endpoint"]["mean_energies_kj_mol"]),
                "reference_mean_energies_kj_mol": dict(run["reference_ff"]["mean_energies_kj_mol"]),
                "mean_energy_delta_kj_mol": dict(run["mean_energy_delta_kj_mol"]),
            }
            for state_name, run in ((state_name, runs_raw[state_name]) for state_name in state_names)
        ),
    }


def _prepared_runtime_payload(
    *,
    evb_system: EVBSystem,
    prepare_manifest: PrepareWorkflowManifest,
    prepare_layout,
) -> dict[str, object]:
    bundle_path = _resolve_parameter_bundle_path(prepare_manifest, prepare_layout)
    bundle_payload = _read_prepare_parameter_bundle(prepare_manifest, prepare_layout)
    create_system_options: Mapping[str, object] = {}
    if bundle_payload is not None:
        force_field_payload = bundle_payload.get("force_field")
        if isinstance(force_field_payload, Mapping):
            raw_options = force_field_payload.get("create_system_options")
            if isinstance(raw_options, Mapping):
                create_system_options = raw_options
    uses_serialized_state_systems = bool(getattr(evb_system, "_state_systems_are_prebuilt", False))
    return _jsonify(
        {
            "mode": "serialized_state_systems" if uses_serialized_state_systems else "prepared_pdb_rebuild",
            "uses_serialized_state_systems": uses_serialized_state_systems,
            "serialized_state_systems_available": _serialized_state_systems_available(
                prepare_manifest.result.prepared_states
            ),
            "parameter_bundle_path": bundle_path,
            "frozen_create_system_options": dict(create_system_options),
        }
    )


def _serialized_state_systems_available(prepared_states) -> bool:
    return bool(prepared_states) and all(entry.system_xml_path is not None for entry in prepared_states)


def _write_completed_stage_resume_state(
    layout,
    manifest: RunWorkflowManifest,
    *,
    substages: list[SubstageStatusRecord],
    generated_files: tuple[StageFileRecord, ...],
) -> Path:
    return write_stage_resume_state(
        layout,
        stage_resume_state_from_manifest(
            manifest,
            status=STAGE_STATUS_COMPLETED,
            substages=tuple(substages),
            generated_files=generated_files,
            detail="Run workflow completed successfully.",
        ),
    )


def _write_failed_stage_resume_state(
    layout,
    manifest: RunWorkflowManifest,
    *,
    current_substage: str,
    completed_substages: list[SubstageStatusRecord],
    exc: Exception,
) -> Path:
    failed_substages = tuple(completed_substages) + (SubstageStatusRecord(current_substage, STAGE_STATUS_FAILED),)
    return write_stage_resume_state(
        layout,
        stage_resume_state_from_manifest(
            manifest,
            status=STAGE_STATUS_FAILED,
            active_substage=current_substage,
            substages=failed_substages,
            detail="Run workflow failed before completion.",
            error_message=f"{type(exc).__name__}: {exc}",
        ),
    )


def _run_stage_output_files(result: RunWorkflowResult) -> tuple[StageFileRecord, ...]:
    records: list[StageFileRecord] = []
    if result.manifest_path is not None:
        records.append(StageFileRecord("run_manifest_json", result.manifest_path))
    records.extend(result.coordinate_files)
    records.extend(result.sample_record_files)
    records.extend(result.summary_files)
    return tuple(records)


def _write_endpoint_summary_csv(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "state_name",
                "seed",
                "evb_sample_count",
                "reference_sample_count",
                "evb_potential_energy_kj_mol",
                "evb_kinetic_energy_kj_mol",
                "evb_total_energy_kj_mol",
                "reference_potential_energy_kj_mol",
                "reference_kinetic_energy_kj_mol",
                "reference_total_energy_kj_mol",
                "delta_potential_energy_kj_mol",
                "delta_kinetic_energy_kj_mol",
                "delta_total_energy_kj_mol",
            ]
        )
        for entry in payload["runs"]:
            evb_mean = entry["evb_mean_energies_kj_mol"]
            reference_mean = entry["reference_mean_energies_kj_mol"]
            delta_mean = entry["mean_energy_delta_kj_mol"]
            writer.writerow(
                [
                    entry["state_name"],
                    "" if entry["seed"] is None else int(entry["seed"]),
                    int(entry["evb_sample_count"]),
                    int(entry["reference_sample_count"]),
                    f"{float(evb_mean['potential_energy_kj_mol']):.8f}",
                    f"{float(evb_mean['kinetic_energy_kj_mol']):.8f}",
                    f"{float(evb_mean['total_energy_kj_mol']):.8f}",
                    f"{float(reference_mean['potential_energy_kj_mol']):.8f}",
                    f"{float(reference_mean['kinetic_energy_kj_mol']):.8f}",
                    f"{float(reference_mean['total_energy_kj_mol']):.8f}",
                    f"{float(delta_mean['potential_energy_kj_mol']):.8f}",
                    f"{float(delta_mean['kinetic_energy_kj_mol']):.8f}",
                    f"{float(delta_mean['total_energy_kj_mol']):.8f}",
                ]
            )


def _write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonify(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")


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


__all__ = [
    "build_evb_system_from_prepare_manifest",
    "run_workflow",
]
