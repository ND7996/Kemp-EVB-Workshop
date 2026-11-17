from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import openmm
from openmm import app, unit

from .prepared_bundle import BUNDLE_VERSION, PreparedBundle, PreparedStateEntry

if TYPE_CHECKING:
    from ..models import MorseBondParameter
    from ..system import EVBSystem
    from ..types import AtomReference, CanonicalPairKey, PathLike


def _positions_with_nanometer_units(positions):
    if unit.is_quantity(positions):
        return positions
    try:
        first = positions[0]
    except (TypeError, IndexError, KeyError):
        return positions * unit.nanometer
    if unit.is_quantity(first):
        flattened = [entry.value_in_unit(unit.nanometer) for entry in positions]
        return flattened * unit.nanometer
    return positions * unit.nanometer


@dataclass
class EVBPreparation:
    """Drive the chemistry-heavy EVB preparation path and emit a PreparedBundle.

    Constructs an :class:`EVBSystem` from the given chemistry inputs, serializes
    per-state topologies and OpenMM ``System`` objects to ``output_directory``,
    and returns both the live system and a portable :class:`PreparedBundle` that
    a cluster-side job can consume via :meth:`EVBSystem.from_prepared_bundle`
    without the OpenFF/RDKit/AmberTools stack.

    This is a driver for the existing ``EVBSystem(...)`` + ``build_state_systems``
    + XML-serialization workflow. It does not replace ``EVBSystem.__init__``.
    """

    state_inputs: Mapping[str, "PathLike"]
    #: "amber14" preset, an OpenMM force-field XML path/name, or a sequence
    #: of them (e.g. ("amber14", "custom_residues.xml")).
    force_field: str | tuple[str, ...] = "amber14"
    atom_mapping: Mapping["AtomReference", "AtomReference"] | None = None
    ligand_charges: Mapping[str, float] | None = None
    ligand_sdf_files: Mapping[str, object] | None = None
    ligand_smiles: Mapping[str, str] | None = None
    #: OPLS-AA/M ffld_server files per ligand residue (path, or {state: path}); overrides the
    #: ligand's LJ + bonded with OPLS after build (charges stay on the ligand_charges/SDF axis).
    ligand_ffld_files: Mapping[str, object] | None = None
    smirnoff_forcefield: str = "openff-2.2.0.offxml"
    charge_method: str | None = None
    parameterization_output: Path | None = None
    verbose: bool = False
    evb_initial_state: str | None = None
    auto_build_forces: bool = True
    softpair_parameters_by_atom_type: Mapping[str, tuple[float, float]] | None = None
    reactive_morse_parameters: (
        "MorseBondParameter | Mapping[CanonicalPairKey, MorseBondParameter] | None"
    ) = None
    reactive_residue_names: "Iterable[str] | None" = None
    create_system_options: Mapping[str, Any] = field(default_factory=dict)

    def build_system(self) -> "EVBSystem":
        """Run the chemistry stack and return an in-memory :class:`EVBSystem`."""

        from ..system import EVBSystem

        return EVBSystem(
            self.state_inputs,
            force_field=self.force_field,
            atom_mapping=self.atom_mapping,
            ligand_charges=self.ligand_charges,
            ligand_sdf_files=self.ligand_sdf_files,
            ligand_smiles=self.ligand_smiles,
            ligand_ffld_files=self.ligand_ffld_files,
            smirnoff_forcefield=self.smirnoff_forcefield,
            charge_method=self.charge_method,
            parameterization_output=self.parameterization_output,
            verbose=self.verbose,
            evb_initial_state=self.evb_initial_state,
            auto_build_forces=self.auto_build_forces,
            softpair_parameters_by_atom_type=self.softpair_parameters_by_atom_type,
            reactive_morse_parameters=self.reactive_morse_parameters,
            reactive_residue_names=self.reactive_residue_names,
        )

    def prepare(
        self,
        output_directory: Path | str,
        *,
        write_bundle_json: bool = True,
        bundle_filename: str = "parameter_bundle.json",
    ) -> tuple["EVBSystem", PreparedBundle]:
        """Build the EVB system, serialize artifacts, and return the bundle."""

        evb_system = self.build_system()
        bundle = self.serialize(
            evb_system,
            output_directory,
            write_bundle_json=write_bundle_json,
            bundle_filename=bundle_filename,
        )
        return evb_system, bundle

    def serialize(
        self,
        evb_system: "EVBSystem",
        output_directory: Path | str,
        *,
        write_bundle_json: bool = True,
        bundle_filename: str = "parameter_bundle.json",
    ) -> PreparedBundle:
        """Serialize an already-built ``EVBSystem`` into a PreparedBundle.

        Use this to capture systems that were modified after construction
        (e.g. solvated via ``solvate_scaas_sphere_states``); :meth:`prepare`
        is a thin convenience wrapper that builds + serializes in one call.
        """

        output_root = Path(output_directory).expanduser().resolve()
        output_root.mkdir(parents=True, exist_ok=True)
        states_root = output_root / "states"
        states_root.mkdir(parents=True, exist_ok=True)

        per_state_systems = evb_system.build_state_systems(**dict(self.create_system_options))

        prepared_states: list[PreparedStateEntry] = []
        for state_name in evb_system.state_names:
            state = evb_system.states[state_name]
            state_dir = states_root / state_name
            state_dir.mkdir(parents=True, exist_ok=True)
            # mmCIF (not PDB): residue ids of large solvated systems exceed the PDB
            # 4-char limit (>9999) and do not survive a PDB round-trip, which breaks
            # the (chain, resid, name) atom-mapping references on bundle load.
            topology_path = state_dir / "prepared.cif"
            system_xml_path = state_dir / "state_system.xml"
            metadata_path = state_dir / "metadata.json"

            with topology_path.open("w", encoding="utf-8") as handle:
                app.PDBxFile.writeFile(
                    state.topology,
                    _positions_with_nanometer_units(state.positions),
                    handle,
                    keepIds=True,
                )
            system_xml_path.write_text(
                openmm.XmlSerializer.serialize(per_state_systems[state_name]),
                encoding="utf-8",
            )
            metadata_path.write_text(
                json.dumps(
                    {
                        "state_name": state_name,
                        "source_pdb_path": str(state.path),
                        "n_atoms": int(state.topology.getNumAtoms()),
                    },
                    indent=2,
                    sort_keys=True,
                ),
                encoding="utf-8",
            )

            prepared_states.append(
                PreparedStateEntry(
                    state_name=state_name,
                    topology_path=topology_path,
                    system_xml_path=system_xml_path,
                    metadata_path=metadata_path,
                )
            )

        atom_mapping_pairs = (
            tuple((source, target) for source, target in self.atom_mapping.items())
            if self.atom_mapping
            else ()
        )

        bundle = PreparedBundle(
            bundle_version=BUNDLE_VERSION,
            state_names=evb_system.state_names,
            prepared_states=tuple(prepared_states),
            atom_mapping=atom_mapping_pairs,
            force_field_name=self.force_field,
            evb_initial_state=self.evb_initial_state,
            reference_state=None,
            auto_build_forces=self.auto_build_forces,
            stored_project_root=output_root,
            parameterization_directory=self.parameterization_output,
            source_path=None,
        )

        if write_bundle_json:
            bundle_path = output_root / bundle_filename
            bundle_path.write_text(
                json.dumps(bundle.to_dict(), indent=2, sort_keys=True),
                encoding="utf-8",
            )

        return bundle
