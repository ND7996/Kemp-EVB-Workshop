from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any
import openmm

from openmm import app, unit

from .. import EVBSystem
from ..core import chain_identifier, residue_sequence_number
from ..input_reader.system_reader import read_pdb_state
from .layout import ensure_stage_layout, layout_from_prepare_manifest
from .models import (
    AtomMappingRecord,
    AtomReferenceRecord,
    PrepareWorkflowManifest,
    PrepareWorkflowResult,
    PrepareWorkflowSpec,
    PreparedStateArtifact,
    StageFileRecord,
)
from .resume import (
    STAGE_STATUS_COMPLETED,
    STAGE_STATUS_FAILED,
    STAGE_STATUS_RUNNING,
    SubstageStatusRecord,
    stage_resume_state_from_manifest,
    write_stage_resume_state,
)

PARAMETER_BUNDLE_TYPE = "evb_parameter_bundle"
PARAMETER_BUNDLE_VERSION = "evb_parameter_bundle_v1"


def _positions_with_nanometer_units(positions):
    if unit.is_quantity(positions):
        return positions
    return positions * unit.nanometer


def run_prepare_workflow(
    manifest: PrepareWorkflowManifest,
    *,
    evb_system_factory: Callable[..., EVBSystem] = EVBSystem,
) -> tuple[EVBSystem, PrepareWorkflowManifest]:
    """Execute the Phase B prepare stage and persist prepared state artifacts."""

    layout = ensure_stage_layout(layout_from_prepare_manifest(manifest))
    spec = manifest.spec
    current_substage = "build_evb_system"
    completed_substages: list[SubstageStatusRecord] = []
    write_stage_resume_state(
        layout,
        stage_resume_state_from_manifest(
            manifest,
            status=STAGE_STATUS_RUNNING,
            active_substage=current_substage,
            substages=(SubstageStatusRecord(current_substage, STAGE_STATUS_RUNNING),),
            detail="Prepare workflow started.",
        ),
    )

    try:
        evb_system = evb_system_factory(
            _state_input_mapping(spec),
            force_field=spec.force_field,
            atom_mapping=_atom_mapping_payload(spec.atom_mapping),
            ligand_charges=_none_if_empty(spec.ligand_charges),
            ligand_sdf_files=_ligand_sdf_payload(spec),
            ligand_smiles=_none_if_empty(spec.ligand_smiles),
            smirnoff_forcefield=spec.smirnoff_forcefield,
            charge_method=spec.charge_method,
            parameterization_output=layout.shared_artifacts_dir / "parameterization",
            verbose=bool(spec.verbose),
            evb_initial_state=spec.evb_initial_state or spec.reference_state,
            auto_build_forces=bool(spec.auto_build_forces),
        )
        completed_substages.append(SubstageStatusRecord(current_substage, STAGE_STATUS_COMPLETED))

        solvation_summary = None
        if _prepare_requires_solvation(spec):
            current_substage = "solvate_states"
            solvation_summary = dict(evb_system.solvate_states(**_solvation_kwargs(spec)))
            completed_substages.append(SubstageStatusRecord(current_substage, STAGE_STATUS_COMPLETED))

        scaas_summary = None
        if spec.scaas_enabled:
            current_substage = "enable_scaas"
            if spec.scaas_radius_angstrom is None:
                raise ValueError("'scaas_radius_angstrom' is required when 'scaas_enabled=True'.")
            scaas_summary = dict(
                evb_system.enable_scaas(
                    solvent_radius_angstrom=float(spec.scaas_radius_angstrom),
                    polarisation=bool(spec.scaas_polarisation),
                    backend=str(spec.scaas_backend),
                    solvent_residue_names=spec.solvent_residue_names or None,
                )
            )
            completed_substages.append(SubstageStatusRecord(current_substage, STAGE_STATUS_COMPLETED))

        current_substage = "write_prepare_outputs"
        prepared_states = tuple(_write_prepared_state_artifacts(evb_system, layout))
        prepared_atom_mapping = _prepared_atom_mapping_records(evb_system, prepared_states)
        summary_path = layout.shared_artifact_path("prepare_summary", ".json")
        bundle_path = layout.shared_artifact_path("parameter_bundle", ".json")
        summary_payload = _prepare_summary_payload(
            manifest=manifest,
            evb_system=evb_system,
            prepared_states=prepared_states,
            solvation_summary=solvation_summary,
            scaas_summary=scaas_summary,
        )
        parameter_bundle_payload = _parameter_bundle_payload(
            manifest=manifest,
            evb_system=evb_system,
            layout=layout,
            prepared_states=prepared_states,
            prepared_atom_mapping=prepared_atom_mapping,
            summary_path=summary_path,
            solvation_summary=solvation_summary,
            scaas_summary=scaas_summary,
        )
        _write_json(summary_path, summary_payload)
        _write_json(bundle_path, parameter_bundle_payload)

        generated_files = (
            StageFileRecord(label="prepare_summary", path=summary_path),
            StageFileRecord(label="parameter_bundle", path=bundle_path),
        )
        result = PrepareWorkflowResult(
            prepared_states=prepared_states,
            atom_mapping=prepared_atom_mapping,
            generated_files=generated_files,
            manifest_path=layout.manifest_path(),
        )
        output_manifest = PrepareWorkflowManifest(
            project=manifest.project,
            spec=manifest.spec,
            result=result,
        )
        output_manifest.write_json(result.manifest_path)
        stage_files = _prepare_stage_output_files(result)
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


def _state_input_mapping(spec: PrepareWorkflowSpec) -> dict[str, Path]:
    return {entry.state_name: entry.pdb_path for entry in spec.state_inputs}


def _atom_mapping_payload(
    atom_mapping: tuple[AtomMappingRecord, ...],
) -> Mapping[tuple[str, str, int, str], tuple[str, str, int, str]] | None:
    if not atom_mapping:
        return None
    return {source: target for source, target in (entry.to_mapping_pair() for entry in atom_mapping)}


def _ligand_sdf_payload(spec: PrepareWorkflowSpec) -> Mapping[str, object] | None:
    if not spec.ligand_sdf_inputs:
        return None
    return {entry.residue_name: entry.to_parameterizer_entry() for entry in spec.ligand_sdf_inputs}


def _prepare_requires_solvation(spec: PrepareWorkflowSpec) -> bool:
    return any(
        (
            spec.water_model is not None,
            spec.solvent_padding_nm is not None,
            spec.box_size_nm is not None,
            spec.ionic_strength_molar is not None,
            spec.neutralize,
            spec.num_added is not None,
            spec.solvent_sphere_radius_angstrom is not None,
        )
    )


def _solvation_kwargs(spec: PrepareWorkflowSpec) -> dict[str, object]:
    kwargs: dict[str, object] = {
        "model": spec.water_model or "tip3p",
        "rebuild_forces": True,
    }
    if spec.reference_state is not None:
        kwargs["reference_state"] = spec.reference_state
    if spec.solvent_padding_nm is not None:
        kwargs["padding"] = float(spec.solvent_padding_nm)
    elif spec.box_size_nm is not None or spec.num_added is not None:
        kwargs["padding"] = None
    if spec.box_size_nm is not None:
        kwargs["box_size"] = spec.box_size_nm
    if spec.ionic_strength_molar is not None:
        kwargs["ionic_strength"] = float(spec.ionic_strength_molar)
    if spec.neutralize:
        kwargs["neutralize"] = True
    if spec.positive_ion != "Na+":
        kwargs["positive_ion"] = spec.positive_ion
    if spec.negative_ion != "Cl-":
        kwargs["negative_ion"] = spec.negative_ion
    if spec.num_added is not None:
        kwargs["num_added"] = int(spec.num_added)
    if spec.solvent_sphere_radius_angstrom is not None:
        kwargs["sphere_radius_angstrom"] = float(spec.solvent_sphere_radius_angstrom)
    if spec.solvent_residue_names:
        kwargs["solvent_residue_names"] = spec.solvent_residue_names
    return kwargs


def _write_prepared_state_artifacts(
    evb_system: EVBSystem,
    layout,
) -> list[PreparedStateArtifact]:
    artifacts: list[PreparedStateArtifact] = []
    state_systems = evb_system.build_state_systems(**_prepare_state_system_create_system_options())
    for state_name, state in evb_system.states.items():
        topology_path = layout.state_artifact_path(state_name, "prepared", ".pdb")
        system_xml_path = layout.state_artifact_path(state_name, "state_system", ".xml")
        metadata_path = layout.state_artifact_path(state_name, "metadata", ".json")
        topology_path.parent.mkdir(parents=True, exist_ok=True)
        system_xml_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        with topology_path.open("w", encoding="utf-8") as handle:
            app.PDBFile.writeFile(
                state.topology,
                _positions_with_nanometer_units(state.positions),
                handle,
            )
        written_state = read_pdb_state(topology_path)
        system_xml_path.write_text(openmm.XmlSerializer.serialize(state_systems[state_name]), encoding="utf-8")
        _write_json(
            metadata_path,
            {
                "state_name": state_name,
                "source_pdb_path": str(state.path),
                "n_atoms": int(written_state.n_atoms),
                "residue_names": sorted(written_state.residue_names),
                "topology_signature_sha256": _topology_signature(written_state.topology),
                "prepared_pdb_path": str(topology_path),
                "system_xml_path": str(system_xml_path),
            },
        )
        artifacts.append(
            PreparedStateArtifact(
                state_name=state_name,
                topology_path=topology_path,
                system_xml_path=system_xml_path,
                metadata_path=metadata_path,
            )
        )
    return artifacts

def _prepared_atom_mapping_records(
    evb_system: EVBSystem,
    prepared_states: tuple[PreparedStateArtifact, ...],
) -> tuple[AtomMappingRecord, ...]:
    translated = _translated_written_atom_references(evb_system, prepared_states)
    records: list[AtomMappingRecord] = []
    grouped: dict[tuple[str, str, int, str], list[AtomReferenceRecord]] = {}
    for atom_reference, canonical_label in evb_system.atom_canonical_labels.items():
        grouped.setdefault(canonical_label, []).append(translated[atom_reference])

    for translated_refs in grouped.values():
        unique_refs = sorted(
            {record.to_atom_reference(): record for record in translated_refs}.values(),
            key=lambda record: record.to_atom_reference(),
        )
        if len(unique_refs) < 2:
            continue
        anchor = unique_refs[0]
        for other in unique_refs[1:]:
            if anchor.state_name == other.state_name:
                continue
            records.append(
                AtomMappingRecord(
                    source=anchor,
                    target=other,
                )
            )
    return tuple(records)


def _write_completed_stage_resume_state(
    layout,
    manifest: PrepareWorkflowManifest,
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
            detail="Prepare workflow completed successfully.",
        ),
    )


def _write_failed_stage_resume_state(
    layout,
    manifest: PrepareWorkflowManifest,
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
            detail="Prepare workflow failed before completion.",
            error_message=f"{type(exc).__name__}: {exc}",
        ),
    )


def _prepare_stage_output_files(result: PrepareWorkflowResult) -> tuple[StageFileRecord, ...]:
    records: list[StageFileRecord] = []
    if result.manifest_path is not None:
        records.append(StageFileRecord("prepare_manifest_json", result.manifest_path))
    records.extend(result.generated_files)
    for state in result.prepared_states:
        if state.topology_path is not None:
            records.append(StageFileRecord(f"{state.state_name}_prepared_pdb", state.topology_path))
        if state.metadata_path is not None:
            records.append(StageFileRecord(f"{state.state_name}_metadata_json", state.metadata_path))
        if state.system_xml_path is not None:
            records.append(StageFileRecord(f"{state.state_name}_system_xml", state.system_xml_path))
        if state.positions_path is not None:
            records.append(StageFileRecord(f"{state.state_name}_positions_file", state.positions_path))
    return tuple(records)


def _prepare_summary_payload(
    *,
    manifest: PrepareWorkflowManifest,
    evb_system: EVBSystem,
    prepared_states: tuple[PreparedStateArtifact, ...],
    solvation_summary: Mapping[str, object] | None,
    scaas_summary: Mapping[str, object] | None,
) -> dict[str, object]:
    return _jsonify(
        {
            "project": manifest.project.to_dict(),
            "spec": manifest.spec.to_dict(),
            "prepared_states": [asdict(entry) for entry in prepared_states],
            "state_atom_counts": {state_name: state.n_atoms for state_name, state in evb_system.states.items()},
            "pdb_paths": [str(path) for path in evb_system.pdb_paths],
            "missing_residues": dict(evb_system.missing_residues),
            "force_field_residue_names": list(evb_system.force_field_residue_names_list),
            "solvation": dict(solvation_summary) if solvation_summary is not None else None,
            "scaas": dict(scaas_summary) if scaas_summary is not None else None,
        }
    )


def _parameter_bundle_payload(
    *,
    manifest: PrepareWorkflowManifest,
    evb_system: EVBSystem,
    layout,
    prepared_states: tuple[PreparedStateArtifact, ...],
    prepared_atom_mapping: tuple[AtomMappingRecord, ...],
    summary_path: Path,
    solvation_summary: Mapping[str, object] | None,
    scaas_summary: Mapping[str, object] | None,
) -> dict[str, object]:
    version_module = getattr(openmm, "version", None)
    return _jsonify(
        {
            "bundle_type": PARAMETER_BUNDLE_TYPE,
            "bundle_version": PARAMETER_BUNDLE_VERSION,
            "project": manifest.project.to_dict(),
            "prepare_spec": manifest.spec.to_dict(),
            "state_names": list(evb_system.state_names),
            "artifacts": {
                "prepare_summary_path": summary_path,
                "parameterization_directory": layout.shared_artifacts_dir / "parameterization",
                "prepared_states": [asdict(entry) for entry in prepared_states],
            },
            "force_field": {
                "name": manifest.spec.force_field,
                "files": [str(path) for path in manifest.spec.force_field_files],
                "create_system_options": _prepare_state_system_create_system_options(),
                "residue_names": list(evb_system.force_field_residue_names_list),
                "missing_residues": dict(evb_system.missing_residues),
            },
            "parameterization": {
                "declared_parameterizer": manifest.spec.parameterizer,
                "parameterizer_inputs": [record.to_dict() for record in manifest.spec.parameterizer_inputs],
                "ligand_charges": dict(manifest.spec.ligand_charges),
                "ligand_sdf_inputs": [entry.to_dict() for entry in manifest.spec.ligand_sdf_inputs],
                "ligand_smiles": dict(manifest.spec.ligand_smiles),
                "smirnoff_forcefield": manifest.spec.smirnoff_forcefield,
                "charge_method": manifest.spec.charge_method,
            },
            "evb": {
                "atom_mapping": [record.to_dict() for record in prepared_atom_mapping],
                "reference_state": manifest.spec.reference_state,
                "evb_initial_state": manifest.spec.evb_initial_state or manifest.spec.reference_state,
                "auto_build_forces": bool(manifest.spec.auto_build_forces),
            },
            "solvation": dict(solvation_summary) if solvation_summary is not None else None,
            "scaas": dict(scaas_summary) if scaas_summary is not None else None,
            "software": {
                "openmm_version": getattr(openmm, "__version__", None),
                "openmm_full_version": getattr(version_module, "full_version", None),
                "openmm_git_revision": getattr(version_module, "git_revision", None),
            },
        }
    )


def _prepare_state_system_create_system_options() -> dict[str, object]:
    return {}


def _topology_signature(topology: app.Topology) -> str:
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


def _none_if_empty(mapping: Mapping[str, object]) -> Mapping[str, object] | None:
    return None if not mapping else dict(mapping)


def _translated_written_atom_references(
    evb_system: EVBSystem,
    prepared_states: tuple[PreparedStateArtifact, ...],
) -> dict[tuple[str, str, int, str], AtomReferenceRecord]:
    translated: dict[tuple[str, str, int, str], AtomReferenceRecord] = {}
    prepared_paths = {
        entry.state_name: entry.topology_path
        for entry in prepared_states
        if entry.topology_path is not None
    }
    for state_name in evb_system.state_names:
        topology_path = prepared_paths.get(state_name)
        if topology_path is None:
            raise ValueError(f"Prepared topology path missing for state '{state_name}'.")
        live_references = evb_system.state_atom_references(state_name)
        written_references = _state_atom_references_from_pdb(state_name, topology_path)
        if len(live_references) != len(written_references):
            raise ValueError(
                f"Prepared PDB for state '{state_name}' changed atom count during export "
                f"({len(live_references)} != {len(written_references)})."
            )
        for live_reference, written_reference in zip(live_references, written_references):
            translated[live_reference] = written_reference
    return translated


def _state_atom_references_from_pdb(state_name: str, pdb_path: Path) -> tuple[AtomReferenceRecord, ...]:
    state = read_pdb_state(pdb_path)
    references: list[AtomReferenceRecord] = []
    for atom in state.topology.atoms():
        residue = atom.residue
        references.append(
            AtomReferenceRecord(
                state_name=state_name,
                chain_id=chain_identifier(residue.chain),
                residue_id=residue_sequence_number(residue, state_name),
                atom_name=(atom.name or "").strip(),
            )
        )
    return tuple(references)


__all__ = ["run_prepare_workflow"]
