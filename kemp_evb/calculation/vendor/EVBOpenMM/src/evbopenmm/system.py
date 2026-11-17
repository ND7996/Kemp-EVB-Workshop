from __future__ import annotations

import os

from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from numbers import Real
from typing import Dict, FrozenSet, Iterable, Optional, Set, Tuple
from types import MappingProxyType
import csv
import json
import math
import warnings

import openmm
from openmm import app, unit

from .analysis import (
    aggregate_evb_profiles,
    analyze_evb_profile_records,
    accumulate_energy_breakdown,
    analyze_evb_profile_replicas_result,
    analyze_mapping_fep_umbrella_records,
    build_evb_energy_breakdown,
    effective_energy_from_state_energies,
    energy_value_kj,
    scale_energy_breakdown,
    energy_scale_to_kj,
    fit_evb_parameters_nstate,
    fit_evb_parameters_replicas_nstate,
    plot_evb_profile_replicas,
    print_energy_debug_summary,
    reaction_coordinate_from_state_energies,
    reduce_reference_terms,
    sample_diabatic_state_energies_kj_mol,
    sample_diabatic_state_energies_via_force_groups_kj_mol,
    split_nonbonded_energy,
    write_evb_profile_outputs,
    write_evb_replica_outputs,
    write_mapping_fep_outputs,
)
from .constants import (
    KJ_PER_KCAL,
    ONE_4PI_EPS0,
    SOLUTE_RESTRAINT_GLOBAL_PARAMETER,
    STATE_ENERGY_TERMS,
)
from .core import (
    atom_reference_from_index,
    atom_reference_from_topology_atom,
    build_atom_canonical_labels,
    build_identity_atom_mapping,
    canonical_angle_key_from_indices,
    canonical_atom_label,
    canonical_bond_key_from_indices,
    canonical_pair_key,
    canonical_torsion_key_from_indices,
    chain_identifier,
    coerce_lambda_sequence,
    find_unmapped_atoms,
    next_auto_name,
    normalize_atom_reference,
    normalize_atom_mapping,
    normalize_residue_name,
    resolve_lambda_parameters,
    particle_indices_from_canonical_key,
    particle_indices_from_canonical_pair,
    residue_sequence_number,
    state_atom_references,
    state_canonical_label_index_map,
    state_particle_masses,
    validate_atom_mapping_atoms,
    validate_state_name,
)
from .forces import (
    add_shared_constraints_to_evb_system,
    adiabatic_force_group_masks,
    apply_context_lambdas,
    apply_indicator_lambdas,
    assign_reference_force_groups,
    build_channel_force,
    lj_sigma_combining_expr,
    collect_bonded_terms_by_state,
    collect_forming_breaking_bond_pairs_by_index,
    build_electro_bond_terms,
    build_lambda_scaled_electro_exception_force,
    build_lambda_scaled_vdw_exception_force,
    build_reactive_harmonic_bond_force,
    build_reactive_morse_bond_force,
    build_reactive_softpair_force,
    build_shared_angle_force,
    build_shared_bond_force,
    build_shared_exception_force,
    build_shared_torsion_force,
    build_state_combined_force,
    build_state_pme_force,
    build_state_angle_force,
    build_state_bond_force,
    build_state_torsion_force,
    build_vdw_bond_terms,
    build_evb_system_scaffold,
    categorize_reference_force,
    collect_group_energies,
    collect_exception_parameters_across_states,
    collect_state_exclusion_pairs,
    collect_global_exclusion_pairs,
    collect_electrostatics_terms_by_state,
    collect_vdw_terms_by_state,
    constraint_distance_nm,
    context_get_parameter_or_default,
    context_lambda_snapshot,
    context_parameter_exists,
    demote_atoms_for_exception_pairs,
    extract_harmonic_bond_terms_by_index,
    extract_electrostatics_exception_parameters,
    extract_electrostatics_particle_parameters,
    extract_vdw_exception_parameters,
    extract_vdw_particle_parameters,
    find_nonbonded_force,
    get_or_create_force_group,
    group_bonded_terms_by_pair,
    parameter_sequences_close,
    pair_is_universal,
    pair_in_qm_channel,
    partition_electrostatics_atoms,
    partition_electrostatics_exceptions,
    partition_vdw_atoms,
    partition_vdw_exceptions,
    partition_zero_exclusions,
    populate_state_combined_particles,
    restore_reference_force_groups,
    restrict_to_subset_interaction_group,
    split_shared_exception_terms,
    split_terms_by_qm_channel,
    state_constraint_map,
    summarize_bonded_terms,
)
from .input_reader.system_reader import EVBState, read_pdb_state
from .models import (
    BondedForceComparison,
    BondedTermEntry,
    ElectrostaticsPartition,
    FEPWindowRecord,
    MDEnergySampleRecord,
    MorseBondParameter,
    ResidueTarget,
    ResidueTargetIndex,
    VanDerWaalsPartition,
)
from .prep import (
    PreparedBundle,
    build_residue_targets,
    compute_missing_force_field_residues,
    ensure_residue_name_consistency,
    format_signature_difference,
    residue_bond_pairs,
    residue_signature,
)
from .scaas import (
    ANGSTROM_TO_NM,
    build_q6_radial_boundary_parameters,
    SCAASBackendController,
    SCAASConfig,
    SCAASRuntime,
    SCAASWaterMolecule,
    build_q6_shells,
    build_scaas_backend_controller,
    normalize_residue_names,
    normalize_scaas_center,
    resolve_scaas_backend,
    summarize_runtime,
)
from .sampling import (
    AdiabaticReporter,
    apply_lambda_stage_sequence,
    build_langevin_middle_adiabatic_integrator,
    build_nve_adiabatic_integrator,
    coerce_lambda_stage_sequence,
    collect_fep_window_records,
    collect_md_energy_samples,
    compute_adiabatic_mixing,
    default_equilibration_schedule,
    estimate_kinetic_dof,
    remove_coordinate_reporter,
    run_equilibration_schedule,
    run_endpoint_state_comparison,
    run_mapping_window_workflow,
    simulation_timestep_ps,
    setup_coordinate_reporting,
    steps_from_report_interval_ps,
    summarize_md_energy_samples,
    validate_velocity_mode,
)
from .types import AtomReference, CanonicalAtomKey, CanonicalPairKey, PathLike
from . import units as _units


def _read_adiabatic_trace_csv(path: Path) -> list[Dict[str, float]]:
    rows: list[Dict[str, float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for raw_row in reader:
            converted: Dict[str, float] = {}
            for key, value in raw_row.items():
                if key == "step":
                    converted[key] = int(value)
                else:
                    converted[key] = float(value)
            rows.append(converted)
    return rows


def _summarize_adiabatic_trace(rows: Sequence[Mapping[str, float]]) -> Dict[str, object]:
    if not rows:
        return {"n_samples": 0}
    n = len(rows)

    def _column(name: str) -> list[float]:
        return [float(row[name]) for row in rows]

    def _mean(values: Sequence[float]) -> float:
        return sum(values) / float(len(values))

    def _stddev(values: Sequence[float], mean: float) -> float:
        if len(values) <= 1:
            return 0.0
        variance = sum((value - mean) ** 2 for value in values) / float(len(values) - 1)
        return math.sqrt(max(variance, 0.0))

    summary: Dict[str, object] = {"n_samples": n}
    for column in ("E_minus_kj_mol", "dV_kj_mol", "c1sq", "c2sq", "H12_kj_mol"):
        values = _column(column)
        mean = _mean(values)
        std = _stddev(values, mean)
        summary[f"mean_{column}"] = mean
        summary[f"std_{column}"] = std
    summary["first_time_ps"] = float(rows[0]["time_ps"])
    summary["last_time_ps"] = float(rows[-1]["time_ps"])
    return summary


class EVBSystem:
    """
    Primary entry point for configuring and running EVB/OpenMM workflows.

    Notes
    -----
    Current EVB experiments are executed in the gas phase without periodic boundary
    conditions. Future work will revisit this construction when spherical boundary
    conditions are introduced.

    Reporting convention
    --------------------
    Any user-facing FEP/umbrella reporting in this library is expressed in
    ``kcal/mol``. Keep this convention when extending reporting APIs.
    """

    def __init__(
        self,
        *state_pdbs: PathLike | Mapping[str, PathLike],
        force_field: str = "amber14",
        lj_combining: str | None = None,
        atom_mapping: Mapping[AtomReference, AtomReference] | None = None,
        ligand_charges: Mapping[str, float] | None = None,
        ligand_sdf_files: Mapping[str, object] | None = None,
        ligand_smiles: Mapping[str, str] | None = None,
        ligand_ffld_files: Mapping[str, object] | None = None,
        smirnoff_forcefield: str = "openff-2.2.0.offxml",
        charge_method: str = "am1bcc",
        parameterization_output: PathLike | None = None,
        verbose: bool = False,
        evb_initial_state: str | None = None,
        auto_build_forces: bool = True,
        prebuilt_state_systems: Mapping[str, openmm.System] | None = None,
        softpair_parameters_by_atom_type: Mapping[str, tuple[float, float]] | None = None,
        reactive_morse_parameters: (
            "MorseBondParameter | Mapping[CanonicalPairKey, MorseBondParameter] | None"
        ) = None,
        reactive_residue_names: "frozenset[str] | set[str] | None" = None,
        **named_state_pdbs: PathLike,
    ) -> None:
        """
        Load EVB states supplied positionally, via mappings, or keyword arguments.

        Parameters
        ----------
        *state_pdbs :
            Positional arguments can be either raw PDB paths (auto-labeled
            ``state_1``, ``state_2``, …) or mapping objects of ``{name: path}``.
        **named_state_pdbs :
            Explicit mapping from user-defined state names (valid identifiers)
            to PDB paths. Each becomes accessible as ``system.<state_name>``.

        force_field :
            Keyword selecting the force field used to parameterize all states.
            Only ``"amber14"`` is currently implemented.
        atom_mapping :
            Optional atom-identity mapping. When provided it must cover every atom across
            all states; otherwise a ``ValueError`` is raised. When omitted, a summary of
            atoms requiring mapping is printed. All atoms appearing in the supplied
            mapping are automatically added to the QM region.
        ligand_charges :
            Optional mapping of residue name -> total charge used when parameterizing missing residues.
            Required if ``atom_mapping`` is provided and the active force field lacks residues.
        ligand_sdf_files :
            Optional mapping of residue name -> SDF/MOL definitions for the missing residues.
        ligand_smiles :
            Optional mapping of residue name -> SMILES strings (fallback if no SDF is given).
        smirnoff_forcefield :
            OFFXML force field used for OpenFF parameterization when needed.
        charge_method :
            Partial-charge assignment method for OpenFF parameterization.
        parameterization_output :
            Optional folder where OpenFF parameterization artifacts are written.
        verbose :
            When ``True``, emit informational messages (e.g., atoms lacking mappings).
        evb_initial_state :
            Name of the EVB state whose coordinates seed the umbrella system (particle ordering,
            reference positions). Defaults to the first registered state when omitted.
        auto_build_forces :
            When ``True`` (default), immediately attach nonbonded and bonded EVB forces after
            building the particle scaffold. Set to ``False`` to configure force-construction
            options manually via ``add_vdw_forces`` / ``add_*_forces`` before creating contexts.
        prebuilt_state_systems :
            Optional mapping of ``state_name -> OpenMM System`` used to initialize the pure-state
            systems directly from serialized artifacts. When provided, missing-residue
            parameterization is skipped and the supplied systems become the authoritative
            per-state force-field definition.

        Example
        -------
        >>> EVBSystem("state_a.pdb", "state_b.pdb")  # auto-labeled states
        >>> EVBSystem({"reactant": "a.pdb", "product": "b.pdb"})
        >>> EVBSystem(reactant="a.pdb", product="b.pdb", ts="c.pdb")

        """
        positional_paths: list[PathLike] = []
        mapped_states: Dict[str, PathLike] = {}

        self.force_field_name = force_field
        # LJ combining rule for the vdW CustomNonbondedForces: 'lorentz-berthelot' (Amber/CHARMM,
        # arithmetic sigma) or 'geometric' (OPLS-AA/M, sqrt sigma). Default (None) auto-selects
        # geometric for OPLS presets, else Lorentz-Berthelot. Validated eagerly.
        if lj_combining is None:
            lj_combining = "geometric" if self.spec_is_opls(force_field) else "lorentz-berthelot"
        lj_sigma_combining_expr(lj_combining)
        self._lj_combining = lj_combining
        # OPLS-AA/M ligand override files (ffld_server output), keyed by residue name; each value
        # is a single path (all states) or a {state_name: path} mapping (e.g. reactant/product).
        # Applied per state after createSystem so a novel ligand gets full OPLS LJ + bonded.
        self._ligand_ffld_files = dict(ligand_ffld_files) if ligand_ffld_files else {}
        if self._ligand_ffld_files and not self.spec_is_opls(force_field):
            raise ValueError(
                "ligand_ffld_files applies OPLS-AA/M ligand parameters and requires an OPLS "
                "force field (force_field='opls-aam'); got a non-OPLS spec, which would silently "
                "ignore the ffld overrides."
            )
        self.force_field = self._load_force_field(force_field)
        self._force_field_residue_names: Set[str] = set(self.force_field._templates.keys())
        # OpenMM normalizes some protonation variants in the topology to a single
        # standard name that is not a literal template key but is still matched by the
        # force field by atom content -- notably histidine, which OpenMM presents as
        # "HIS" while force fields template it under either the Amber (HID/HIE/HIP) or
        # CHARMM/OPLS (HSD/HSE/HSP) naming family. Register "HIS" whenever either family
        # is present so such protein residues are recognized as force-field
        # parameterized rather than mistaken for custom (OpenFF) ligands; OpenMM's
        # createSystem still graph-matches the actual protonation to the right template.
        if {"HID", "HIE", "HIP", "HSD", "HSE", "HSP"} & self._force_field_residue_names:
            self._force_field_residue_names.add("HIS")
        # Monatomic counterions: OpenMM's addSolvent (called during solvation) names them by
        # the plain/Amber convention (NA, CL, K), while CHARMM/OPLS template them under SOD/CLA/
        # POT. Register the plain names whenever the corresponding ion template is present so
        # neutralizing counterions are recognized as force-field-parameterized rather than
        # mistaken for custom ligands; createSystem still graph-matches the single ion atom.
        for _tmpl, _aliases in (("SOD", ("NA", "NA+")), ("CLA", ("CL", "CL-")), ("POT", ("K", "K+"))):
            if _tmpl in self._force_field_residue_names:
                self._force_field_residue_names.update(_aliases)
        self._normalized_force_field_residues: Set[str] = {
            self._normalize_residue_name(name) for name in self._force_field_residue_names
        }
        initial_generators = getattr(self.force_field, "_templateGenerators", None) or []
        self._base_template_generators: Tuple[object, ...] = tuple(initial_generators)
        self._template_generator_specs: Tuple[object, ...] = tuple()

        for arg in state_pdbs:
            if isinstance(arg, Mapping):
                for name, pdb_path in arg.items():
                    if name in mapped_states:
                        raise ValueError(f"State name '{name}' provided multiple times.")
                    mapped_states[name] = pdb_path
            else:
                positional_paths.append(arg)

        all_named: Dict[str, PathLike] = {}
        for source in (mapped_states, named_state_pdbs):
            for name, pdb_path in source.items():
                self._validate_state_name(name)
                if name in all_named:
                    raise ValueError(f"State name '{name}' provided multiple times.")
                all_named[name] = pdb_path

        total_states = len(positional_paths) + len(all_named)
        if total_states < 2:
            raise ValueError("EVBSystem requires at least two states (two PDB files).")

        self._states: Dict[str, EVBState] = {}
        reserved_names = set(all_named.keys())

        for pdb_path in positional_paths:
            name = self._next_auto_name(reserved_names | set(self._states.keys()))
            state = read_pdb_state(pdb_path)
            self._states[name] = state
            setattr(self, name, state)

        for name, pdb_path in all_named.items():
            if name in self._states:
                raise ValueError(f"State name '{name}' provided multiple times.")
            state = read_pdb_state(pdb_path)
            self._states[name] = state
            setattr(self, name, state)

        self._modellers: Dict[str, app.Modeller] = {}
        self._topology_bonds: Dict[str, Tuple[Tuple[app.Topology.Atom, app.Topology.Atom], ...]] = {}
        self._residue_targets: Dict[str, list[ResidueTarget]] = {}
        self._state_systems: Dict[str, openmm.System] = {}
        self._state_systems_cache: Mapping[str, openmm.System] | None = None
        self._state_system_signature: Tuple[Tuple[str, object], ...] | None = None
        self._state_systems_are_prebuilt: bool = False
        self._state_particle_mass_lists: Dict[str, Tuple[unit.Quantity, ...]] = {}
        self._openff_reference_molecules: Dict[str, object] = {}
        self._state_atom_reference_lists: Dict[str, Tuple[AtomReference, ...]] = {}
        self._state_atom_identifier_cache: Dict[str, Set[AtomReference]] = {}
        self._state_canonical_index_maps: Dict[str, Dict[AtomReference, int]] = {}
        self._evb_system_scaffold: openmm.System | None = None
        self._evb_system_source_state: str | None = None
        self._evb_initial_positions: Tuple[app.Vec3, ...] | None = None
        self._shared_bond_force: openmm.HarmonicBondForce | None = None
        self._state_specific_bond_forces: Dict[str, openmm.CustomBondForce] = {}
        self._shared_angle_force: openmm.HarmonicAngleForce | None = None
        self._state_specific_angle_forces: Dict[str, openmm.CustomAngleForce] = {}
        self._shared_torsion_force: openmm.PeriodicTorsionForce | None = None
        self._state_specific_torsion_forces: Dict[str, openmm.CustomTorsionForce] = {}
        self._shared_nonbonded_force: openmm.CustomNonbondedForce | None = None
        self._shared_nonbonded_qm_force: openmm.CustomNonbondedForce | None = None
        self._shared_nonbonded_exception_force: openmm.CustomBondForce | None = None
        self._shared_nonbonded_qm_exception_force: openmm.CustomBondForce | None = None
        self._state_specific_nonbonded_forces: Dict[str, openmm.CustomNonbondedForce] = {}
        self._state_specific_nonbonded_qm_forces: Dict[str, openmm.CustomNonbondedForce] = {}
        self._state_specific_vdw_exception_forces: Dict[str, openmm.CustomBondForce] = {}
        self._state_specific_vdw_qm_exception_forces: Dict[str, openmm.CustomBondForce] = {}
        self._state_specific_reactive_vdw_forces: Dict[str, openmm.CustomBondForce] = {}
        self._state_specific_reactive_morse_bond_forces: Dict[str, openmm.CustomBondForce] = {}
        # Endpoint-integration-only harmonic bonds matched to the Morse minimum.
        # Off (global param ``g_regular_ff`` = 0) for all energy/analysis; switched
        # on by ``restrict_integration_to_endpoint(regular_ff=True)`` and selected
        # for integration in place of the Morse (which is 2 fs-unstable). Each lives
        # in a dedicated force group that is NOT in the (state, term) energy registry,
        # so ``_state_force_group_masks`` / the diabatic-energy sampler never see it.
        self._endpoint_regular_ff_forces: Dict[str, openmm.CustomBondForce] = {}
        self._endpoint_regular_ff_groups: Dict[str, int] = {}
        self._regular_ff_global_param: str = "g_regular_ff"
        self._state_specific_electrostatics_exception_forces: Dict[str, openmm.CustomBondForce] = {}
        self._state_specific_electrostatics_qm_exception_forces: Dict[str, openmm.CustomBondForce] = {}
        # Periodic (PME) electrostatics: one Coulomb-only NonbondedForce per diabatic
        # state, plus the name of the global parameter (driven to sqrt(lambda_state))
        # that scales each state's charges.
        self._state_specific_pme_forces: Dict[str, openmm.NonbondedForce] = {}
        self._state_pme_g_param_names: Dict[str, str] = {}
        self._qm_region_atoms: Set[AtomReference] = set()
        self._mm_nonbonded_cutoff_nm: float | None = None
        self._qm_nonbonded_cutoff_nm: float | None = None
        # Boundary condition mode for the EVB scaffold. "nonperiodic" (gas phase or
        # SCAAS droplet, CutoffNonPeriodic) is the historical default; "periodic"
        # enables a real periodic box with PME electrostatics; "scaas_sphere" marks a
        # trimmed spherical droplet whose box vectors are intentionally cleared.
        self._boundary_mode: str = "nonperiodic"
        # Water model requested at the most recent solvation call. Selects which water
        # template the OPLS preset loads (tip3p.xml vs spce.xml) so the added solvent is
        # parameterized by the requested model rather than silently defaulting to TIP3P.
        self._solvation_water_model: str = "tip3p"
        # Electrostatics evaluation scheme in periodic mode: "per_state" (one PME force
        # per diabatic state, exact, N FFTs/step) or "mixed" (single mixed-charge PME
        # plus a sparse diabatic correction, ~1 FFT/step).
        self._electrostatics_mode: str = "per_state"
        self._ewald_error_tolerance: float = 5e-4
        self._scaas_runtime: SCAASRuntime | None = None
        self._scaas_backend_controller: SCAASBackendController | None = None
        self._force_group_ids: Dict[Tuple[str | None, str], int] = {}
        self._force_group_definitions: Dict[int, Tuple[str | None, str]] = {}
        self._next_force_group_id: int = 0
        self._vdw_partition: VanDerWaalsPartition | None = None
        self._electrostatics_partition: ElectrostaticsPartition | None = None
        self._vdw_partition_cache_key: Tuple[Tuple[Tuple[str, object], ...] | None, float] | None = None
        self._electrostatics_partition_cache_key: Tuple[Tuple[Tuple[str, object], ...] | None, float] | None = None
        self._state_lambda_parameter_names: Dict[str, str] = {}
        self.bonded_force_comparisons: Mapping[str, BondedForceComparison] | None = None
        # Residues whose intramolecular bonds legitimately differ between diabatic
        # states (intramolecular reactive residues). The OpenFF build path infers
        # these from per-state reference molecules, but the prepared-bundle path has
        # no OpenFF graphs, so callers declare them explicitly here.
        self._reactive_residue_names: frozenset[str] = frozenset(reactive_residue_names or ())
        # Populated by add_state_dependent_onethree_exclusion_patches: pairs whose
        # 1-4 classification differs between diabats (now restored with per-state
        # scaled 1-4 parameters rather than full interactions).
        self._state_dependent_onefour_pairs: frozenset = frozenset()
        # Whether add_state_dependent_onethree_exclusion_patches has run. Tracked
        # separately from _state_dependent_onefour_pairs because that set is
        # legitimately empty when no pair changes 1-4 class, so it cannot
        # distinguish "ran and found nothing" from "never ran".
        self._onethree_patch_applied: bool = False
        self._build_residue_targets()
        self._ensure_residue_name_consistency(include_bonds=True)
        self._missing_residues = self._compute_missing_force_field_residues()

        self._verbose = bool(verbose)

        self._atom_mapping: Dict[AtomReference, AtomReference] = {}
        qm_seed_atoms: Tuple[AtomReference, ...] = tuple()
        if atom_mapping is not None:
            self._validate_atom_mapping_atoms(atom_mapping)
            self.ensure_atom_mapping_complete(atom_mapping)
            self._atom_mapping = self._normalize_atom_mapping(atom_mapping)
            provided_refs: list[AtomReference] = []
            for left_ref, right_ref in atom_mapping.items():
                provided_refs.append(left_ref)
                provided_refs.append(right_ref)
            qm_seed_atoms = tuple(provided_refs)
        else:
            missing_atoms = self.find_unmapped_atoms({})
            if missing_atoms:
                lines = ["Atom mapping required for the following atoms:"]
                for state_name, atoms in missing_atoms.items():
                    lines.append(f"  {state_name}:")
                    for atom_ref in atoms:
                        lines.append(f"    {atom_ref!r}")
                self._emit("\n".join(lines))
            self._atom_mapping = self._build_identity_atom_mapping()

        self._atom_canonical_labels = self._build_atom_canonical_labels()
        if qm_seed_atoms:
            self.set_qm_atoms(qm_seed_atoms)

        if prebuilt_state_systems is None:
            self._auto_parameterize_missing_residues(
                ligand_charges=ligand_charges,
                ligand_sdf_files=ligand_sdf_files,
                ligand_smiles=ligand_smiles,
                smirnoff_forcefield=smirnoff_forcefield,
                charge_method=charge_method,
                output_directory=parameterization_output,
            )
        else:
            self._load_prebuilt_state_systems(prebuilt_state_systems)
        self._state_systems_cache = self.build_state_systems()
        self.bonded_force_comparisons = self.compare_bonded_terms()
        selected_initial_state = evb_initial_state or self.state_names[0]
        if selected_initial_state not in self._states:
            raise ValueError(
                f"EVB initial state '{selected_initial_state}' not found among defined states {self.state_names}."
            )
        self.build_evb_system(initial_state=selected_initial_state)
        # Retained so build_reactive_forces() can default to the constructor's
        # values instead of making callers repeat them.
        self._constructor_reactive_morse_parameters = reactive_morse_parameters
        self._constructor_softpair_parameters = softpair_parameters_by_atom_type

        if auto_build_forces:
            self.add_vdw_forces()
            self.add_bond_forces()
            self.add_angle_forces()
            self.add_torsion_forces()
            # Parity-validated EVB form: always uses the binary-accurate
            # soft-pair on forming/breaking reactive pairs, regardless of the
            # base force field used to parameterize the rest of the system.
            # Skipped silently when the system has no reactive pairs (e.g. a
            # static comparison setup with identical state topologies).
            if self.get_forming_breaking_bond_pairs():
                if reactive_morse_parameters is None:
                    raise ValueError(
                        "EVBSystem detected forming/breaking reactive pairs but "
                        "``reactive_morse_parameters`` was not provided. EVBOpenMM "
                        "requires Morse parameters (``r0_nm``, ``de_kj_mol``, "
                        "``alpha_per_nm``) for every reactive pair: a force-field "
                        "parameterisation step (OpenFF / antechamber / amber) only "
                        "produces a harmonic bond, which cannot capture bond "
                        "dissociation. Fit Morse parameters to a QM dissociation "
                        "curve (or import from a Q6 ``.fep`` file) and pass them "
                        "via this kwarg. Pass either a single ``MorseBondParameter`` "
                        "(applied to every reactive pair) or a "
                        "``Mapping[CanonicalPairKey, MorseBondParameter]``."
                    )
                self.replace_reactive_bonds_with_morse(reactive_morse_parameters)
                self.add_reactive_softpair_vdw_forces(
                    soft_pair_parameters_by_atom_type=softpair_parameters_by_atom_type,
                )

    @classmethod
    def from_prepared_state_systems(
        cls,
        state_inputs: Mapping[str, PathLike],
        *,
        state_system_xmls: Mapping[str, PathLike],
        force_field: str = "amber14",
        atom_mapping: Mapping[AtomReference, AtomReference] | None = None,
        verbose: bool = False,
        evb_initial_state: str | None = None,
        auto_build_forces: bool = True,
        softpair_parameters_by_atom_type: Mapping[str, tuple[float, float]] | None = None,
        reactive_morse_parameters: (
            "MorseBondParameter | Mapping[CanonicalPairKey, MorseBondParameter] | None"
        ) = None,
        reactive_residue_names: "frozenset[str] | set[str] | None" = None,
    ) -> "EVBSystem":
        """
        Construct an EVB system from prepared PDB topologies plus serialized state systems.

        This bypasses small-molecule parameterization and reuses the provided OpenMM
        ``System`` objects as the authoritative pure-state force-field definitions.

        ``reactive_residue_names`` declares residues whose intramolecular bonds differ
        between diabatic states (intramolecular reactive residues); without OpenFF
        graphs to infer them, the bundle path needs them spelled out to skip the
        bond-consistency check for those residues.
        """

        return cls(
            state_inputs,
            force_field=force_field,
            atom_mapping=atom_mapping,
            verbose=verbose,
            evb_initial_state=evb_initial_state,
            auto_build_forces=auto_build_forces,
            prebuilt_state_systems=cls._deserialize_state_system_xmls(state_system_xmls),
            softpair_parameters_by_atom_type=softpair_parameters_by_atom_type,
            reactive_morse_parameters=reactive_morse_parameters,
            reactive_residue_names=reactive_residue_names,
        )

    @classmethod
    def from_prepared_bundle(
        cls,
        bundle: "PreparedBundle",
        *,
        verbose: bool = False,
        softpair_parameters_by_atom_type: Mapping[str, tuple[float, float]] | None = None,
        reactive_morse_parameters: (
            "MorseBondParameter | Mapping[CanonicalPairKey, MorseBondParameter] | None"
        ) = None,
        reactive_residue_names: "frozenset[str] | set[str] | None" = None,
    ) -> "EVBSystem":
        """Construct an EVB system from a :class:`PreparedBundle` loaded from disk.

        ``reactive_residue_names`` declares intramolecular reactive residues whose
        bonds differ between diabatic states (e.g. ``{"SUB"}`` for a single-residue
        substrate that rearranges); the bundle path cannot infer them, so pass them
        explicitly to skip the per-residue bond-consistency check.
        """

        return cls.from_prepared_state_systems(
            bundle.state_topology_paths(),
            state_system_xmls=bundle.state_system_xml_paths(),
            force_field=bundle.force_field_name,
            atom_mapping=bundle.atom_mapping_payload(),
            verbose=verbose,
            evb_initial_state=bundle.evb_initial_state,
            auto_build_forces=bundle.auto_build_forces,
            softpair_parameters_by_atom_type=softpair_parameters_by_atom_type,
            reactive_morse_parameters=reactive_morse_parameters,
            reactive_residue_names=reactive_residue_names,
        )

    @property
    def evb_scaffold(self) -> openmm.System | None:
        """
        The current EVB scaffold ``openmm.System``, or ``None`` before
        :meth:`build_evb_system` has been called.

        When ``auto_build_forces=True`` is passed to :meth:`__init__` (or
        :meth:`from_prepared_bundle` / :meth:`from_prepared_state_systems`),
        the scaffold is built and populated with vdW, electrostatic, bonded,
        and torsion forces during construction. Calling
        :meth:`build_evb_system` later rebuilds the scaffold from scratch and
        clears any previously-attached forces.
        """
        return self._evb_system_scaffold

    @property
    def states(self) -> Mapping[str, EVBState]:
        """Read-only mapping of state name -> EVBState."""
        return MappingProxyType(self._states)

    @property
    def state_names(self) -> Tuple[str, ...]:
        return tuple(self._states.keys())

    @property
    def pdb_paths(self) -> Tuple[Path, ...]:
        return tuple(state.path for state in self._states.values())

    @property
    def residue_names_by_state(self) -> Mapping[str, FrozenSet[str]]:
        """Read-only mapping of state name -> unique residue names present in that state."""
        return MappingProxyType({name: state.residue_names for name, state in self._states.items()})

    @property
    def atom_canonical_labels(self) -> Mapping[AtomReference, AtomReference]:
        """Read-only mapping of atom reference -> canonical atom label."""
        return MappingProxyType(self._atom_canonical_labels)

    def state_atom_references(self, state_name: str) -> Tuple[AtomReference, ...]:
        """Return the ordered atom references for ``state_name``."""
        return self._state_atom_references(state_name)

    def set_qm_atoms(self, atoms: Sequence[AtomReference]) -> Tuple[AtomReference, ...]:
        """
        Define the QM region using atom references from any EVB state.

        Parameters
        ----------
        atoms :
            Sequence of atom references in ``(state_name, chain_id, residue_id, atom_name)`` form.
            References are normalized and mapped to canonical atom labels.

        Returns
        -------
        Tuple[AtomReference, ...]
            Sorted canonical labels currently assigned to the QM region.
        """

        if (
            self._shared_nonbonded_force is not None
            or self._shared_nonbonded_qm_force is not None
            or self._shared_nonbonded_exception_force is not None
            or self._shared_nonbonded_qm_exception_force is not None
            or self._state_specific_nonbonded_forces
            or self._state_specific_nonbonded_qm_forces
            or self._state_specific_vdw_exception_forces
            or self._state_specific_vdw_qm_exception_forces
            or self._state_specific_electrostatics_exception_forces
            or self._state_specific_electrostatics_qm_exception_forces
        ):
            raise RuntimeError(
                "QM region cannot be changed after nonbonded forces are added. "
                "Rebuild the EVB scaffold (or call 'solvate_states(..., rebuild_forces=True)') "
                "before changing QM atoms."
            )

        qm_labels: Set[AtomReference] = set()
        for atom_ref in atoms:
            normalized = self._normalize_atom_reference(atom_ref)
            canonical = self._atom_canonical_labels.get(normalized, normalized)
            qm_labels.add(canonical)
        self._qm_region_atoms = qm_labels
        return tuple(sorted(self._qm_region_atoms))

    def clear_qm_atoms(self) -> None:
        """Clear the currently configured QM atom set."""

        if (
            self._shared_nonbonded_force is not None
            or self._shared_nonbonded_qm_force is not None
            or self._shared_nonbonded_exception_force is not None
            or self._shared_nonbonded_qm_exception_force is not None
            or self._state_specific_nonbonded_forces
            or self._state_specific_nonbonded_qm_forces
            or self._state_specific_vdw_exception_forces
            or self._state_specific_vdw_qm_exception_forces
            or self._state_specific_electrostatics_exception_forces
            or self._state_specific_electrostatics_qm_exception_forces
        ):
            raise RuntimeError(
                "QM region cannot be changed after nonbonded forces are added. "
                "Rebuild the EVB scaffold before clearing QM atoms."
            )
        self._qm_region_atoms = set()

    def get_qm_atoms(self) -> Tuple[AtomReference, ...]:
        """Return sorted canonical atom labels currently assigned to the QM region."""

        return tuple(sorted(self._qm_region_atoms))

    def residue_bonds_by_state(
        self,
        residue_name: str | Tuple[str, str, int],
        *,
        state_name: str | None = None,
        chain_id: str | None = None,
        residue_id: int | str | None = None,
    ) -> Mapping[str, Tuple[CanonicalPairKey, ...]]:
        """Return the bonds associated with ``residue_name`` for each matching EVB state.

        Parameters
        ----------
        residue_name
            Either the residue name (case-insensitive) or a 3-tuple
            ``(state_name, chain_id, residue_id)`` identifying the residue.
        state_name
            Optional state filter. When omitted, all states containing the residue are reported.
        chain_id
            Optional chain identifier (matching the PDB chain id). Required when multiple
            chains contain the same residue name to avoid ambiguity.
        residue_id
            Optional residue sequence number (PDB ``resSeq``). Provide this when the same
            residue name appears multiple times within the filtered selection.

        Returns
        -------
        Mapping[str, Tuple[CanonicalPairKey, ...]]
            Mapping of state name -> sorted tuple of bonds. Each bond is expressed as
            ``((state, chain, residue_id, atom_name), (state, chain, residue_id, atom_name))``
            using canonical ordering.

        Raises
        ------
        ValueError
            If the residue cannot be found, if an invalid filter is supplied, or if the
            selection remains ambiguous across the requested states.
        """

        tuple_state_filter: str | None = None
        tuple_chain_filter: str | None = None
        tuple_residue_filter: int | None = None

        def _locate_context(state: str, chain: str, resid: int) -> Tuple[str, list[ResidueTarget]] | None:
            for name, entries in self._residue_targets.items():
                for target in entries:
                    if target.state_name != state:
                        continue
                    if self._chain_identifier(target.residue.chain) != chain:
                        continue
                    if self._residue_sequence_number(target.residue, target.state_name) != resid:
                        continue
                    return name, entries
            return None

        if isinstance(residue_name, tuple):
            if len(residue_name) != 3:
                raise ValueError(
                    "Tuple residue references must be (state_name, chain_id, residue_id)."
                )
            ref_state, ref_chain, ref_resid = residue_name
            if not isinstance(ref_state, str) or ref_state not in self._states:
                raise ValueError(f"State '{ref_state}' not found for residue reference {residue_name!r}.")
            if not isinstance(ref_chain, str) or not ref_chain.strip():
                raise ValueError("Tuple residue references must include a non-empty chain identifier.")
            cleaned_chain = ref_chain.strip()
            try:
                cleaned_resid = int(ref_resid)
            except (TypeError, ValueError) as exc:
                raise ValueError("Tuple residue references require an integer residue id.") from exc

            located = _locate_context(ref_state, cleaned_chain, cleaned_resid)
            if located is None:
                raise ValueError(
                    f"Residue reference {residue_name!r} does not match any residue in the EVB states."
                )
            normalized_name, contexts = located
            tuple_state_filter = ref_state
            tuple_chain_filter = cleaned_chain
            tuple_residue_filter = cleaned_resid
        else:
            normalized_name = self._normalize_residue_name(residue_name)
            contexts = self._residue_targets.get(normalized_name)
            if not contexts:
                raise ValueError(f"Residue '{residue_name}' not present in any EVB state.")

        if state_name is not None:
            if state_name not in self._states:
                raise ValueError(
                    f"State '{state_name}' not found. Available states: {', '.join(self.state_names)}."
                )
            state_filter = state_name
            if tuple_state_filter and state_filter != tuple_state_filter:
                raise ValueError("Conflicting state filters supplied via tuple reference and 'state_name'.")
        else:
            state_filter = tuple_state_filter

        if chain_id is not None:
            if not isinstance(chain_id, str) or not chain_id.strip():
                raise ValueError("'chain_id' must be a non-empty string when provided.")
            chain_filter = chain_id.strip()
            if tuple_chain_filter and chain_filter != tuple_chain_filter:
                raise ValueError("Conflicting chain filters supplied via tuple reference and 'chain_id'.")
        else:
            chain_filter = tuple_chain_filter

        if residue_id is not None:
            try:
                residue_filter = int(residue_id)
            except (TypeError, ValueError) as exc:
                raise ValueError("'residue_id' must be an integer or integer-like string.") from exc
            if tuple_residue_filter is not None and residue_filter != tuple_residue_filter:
                raise ValueError("Conflicting residue ids supplied via tuple reference and 'residue_id'.")
        else:
            residue_filter = tuple_residue_filter

        matches: Dict[str, list[Tuple[ResidueTarget, str, int]]] = {}
        for target in contexts:
            if state_filter and target.state_name != state_filter:
                continue
            target_chain = self._chain_identifier(target.residue.chain)
            if chain_filter and target_chain != chain_filter:
                continue
            target_resid = self._residue_sequence_number(target.residue, target.state_name)
            if residue_filter is not None and target_resid != residue_filter:
                continue
            matches.setdefault(target.state_name, []).append((target, target_chain, target_resid))

        if not matches:
            parts = []
            if state_filter:
                parts.append(f"state='{state_filter}'")
            if chain_filter:
                parts.append(f"chain='{chain_filter}'")
            if residue_filter is not None:
                parts.append(f"residue_id={residue_filter}")
            qualifiers = f" ({', '.join(parts)})" if parts else ""
            raise ValueError(
                f"Residue '{residue_name}' matching the supplied qualifiers{qualifiers or ''} was not found."
            )

        ambiguous: Dict[str, list[Tuple[str, int]]] = {}
        for state, entries in matches.items():
            if len(entries) > 1:
                ambiguous[state] = [(chain, resid) for _, chain, resid in entries]
        if ambiguous:
            descriptions = []
            for state, options in sorted(ambiguous.items()):
                formatted = ", ".join(f"{chain}:{resid}" for chain, resid in options)
                descriptions.append(f"{state} -> [{formatted}]")
            details = "; ".join(descriptions)
            raise ValueError(
                "Residue selection is ambiguous. Provide 'chain_id' and/or 'residue_id' to "
                f"disambiguate. Candidates: {details}"
            )

        ordering = {name: idx for idx, name in enumerate(self.state_names)}
        ordered_matches = sorted(matches.items(), key=lambda item: ordering[item[0]])
        bonds_by_state: Dict[str, Tuple[CanonicalPairKey, ...]] = {}
        for state, entries in ordered_matches:
            target, _, _ = entries[0]
            bonds_by_state[state] = self._residue_bond_pairs(target)

        return MappingProxyType(bonds_by_state)

    @property
    def missing_residues(self) -> Mapping[str, Tuple[str, ...]]:
        """Residues lacking parameters in the currently selected force field."""
        return MappingProxyType(self._missing_residues)

    @property
    def force_field_residue_names_list(self) -> Tuple[str, ...]:
        """Sorted list of residue names available in the active force field."""
        return tuple(sorted(self._force_field_residue_names))

    @property
    def evb_initial_positions(self) -> Tuple[app.Vec3, ...]:
        """
        Reference coordinates used when constructing the EVB umbrella system.

        Returns
        -------
        Tuple[app.Vec3, ...]
            Positions copied from the EVB state passed to :meth:`build_evb_system`.

        Raises
        ------
        ValueError
            If :meth:`build_evb_system` has not been called yet.
        """

        if self._evb_initial_positions is None:
            raise ValueError("EVB base system not constructed yet; call 'build_evb_system'.")
        return self._evb_initial_positions

    @staticmethod
    def _normalize_residue_name(name: str) -> str:
        return normalize_residue_name(name)

    def find_unmapped_atoms(self, mapping: Mapping[AtomReference, AtomReference]) -> Dict[str, Tuple[AtomReference, ...]]:
        """
        Identify atoms that are absent from a user-provided atom mapping.

        Parameters
        ----------
        mapping :
            Mapping describing how atoms change identity between EVB states.
            Each key/value is a 4-tuple ``(state_name, chain_id, residue_id, atom_name)``.

        Returns
        -------
        Dict[str, Tuple[AtomReference, ...]]
            Mapping of state name -> sorted tuple of unmapped atom identifiers.
        """
        return find_unmapped_atoms(
            state_names=self.state_names,
            state_atom_references=self._state_atom_references,
            normalize_atom_reference=self._normalize_atom_reference,
            mapping=mapping,
        )

    def ensure_atom_mapping_complete(self, mapping: Mapping[AtomReference, AtomReference]) -> None:
        """
        Validate that a user-provided mapping covers every atom in each EVB state.

        Parameters
        ----------
        mapping :
            Atom mapping supplied by the user.

        Raises
        ------
        ValueError
            If any atom was omitted from the mapping.
        """
        missing = self.find_unmapped_atoms(mapping)
        if not missing:
            return

        details: list[str] = []
        for state_name, atoms in missing.items():
            preview = ", ".join(repr(atom) for atom in atoms[:5])
            suffix = "" if len(atoms) <= 5 else ", ..."
            details.append(f"{state_name}: {preview}{suffix}")
        detail_str = "; ".join(details)
        raise ValueError(
            "Atom mapping is incomplete. Provide entries for every atom. Missing examples -> " + detail_str
        )

    def _validate_atom_mapping_atoms(self, mapping: Mapping[AtomReference, AtomReference]) -> None:
        """
        Ensure every atom reference mentioned in ``mapping`` exists in the loaded PDB states.
        """
        validate_atom_mapping_atoms(mapping, normalize_atom_reference=self._normalize_atom_reference)

    def _normalize_atom_mapping(self, mapping: Mapping[AtomReference, AtomReference]) -> Dict[AtomReference, AtomReference]:
        """
        Create a bidirectional normalized atom mapping covering every state atom.
        """
        return normalize_atom_mapping(
            mapping,
            state_names=self.state_names,
            state_atom_references=self._state_atom_references,
            normalize_atom_reference=self._normalize_atom_reference,
        )

    def _build_identity_atom_mapping(self) -> Dict[AtomReference, AtomReference]:
        return build_identity_atom_mapping(
            state_names=self.state_names,
            state_atom_references=self._state_atom_references,
        )

    def _build_atom_canonical_labels(self) -> Dict[AtomReference, AtomReference]:
        return build_atom_canonical_labels(self._atom_mapping)

    def _normalize_atom_reference(self, ref: AtomReference) -> AtomReference:
        return normalize_atom_reference(
            ref,
            states=self._states,
            state_atom_identifier_cache=self._state_atom_identifier_cache,
            state_atom_references=self._state_atom_references,
        )

    def _atom_reference_from_index(self, state_name: str, atom_index: int) -> AtomReference:
        return atom_reference_from_index(
            state_name=state_name,
            atom_index=atom_index,
            state_atom_references=self._state_atom_references,
        )

    def _canonical_atom_label(self, state_name: str, atom_index: int) -> AtomReference:
        return canonical_atom_label(
            state_name=state_name,
            atom_index=atom_index,
            atom_reference_from_index=self._atom_reference_from_index,
            atom_canonical_labels=self._atom_canonical_labels,
        )

    def _canonical_bond_key_from_indices(self, state_name: str, atom1: int, atom2: int) -> CanonicalAtomKey:
        return canonical_bond_key_from_indices(
            state_name=state_name,
            atom1=atom1,
            atom2=atom2,
            canonical_atom_label=self._canonical_atom_label,
        )

    def _canonical_angle_key_from_indices(self, state_name: str, atom1: int, atom2: int, atom3: int) -> CanonicalAtomKey:
        return canonical_angle_key_from_indices(
            state_name=state_name,
            atom1=atom1,
            atom2=atom2,
            atom3=atom3,
            canonical_atom_label=self._canonical_atom_label,
        )

    def _canonical_torsion_key_from_indices(
        self, state_name: str, atom1: int, atom2: int, atom3: int, atom4: int
    ) -> CanonicalAtomKey:
        return canonical_torsion_key_from_indices(
            state_name=state_name,
            atom1=atom1,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            canonical_atom_label=self._canonical_atom_label,
        )

    def _canonical_pair_key(self, first: AtomReference, second: AtomReference) -> CanonicalPairKey:
        return canonical_pair_key(first, second)

    def _atom_reference_from_topology_atom(
        self, state_name: str, atom: app.Topology.Atom
    ) -> AtomReference:
        return atom_reference_from_topology_atom(
            state_name=state_name,
            atom=atom,
            chain_identifier=self._chain_identifier,
            residue_sequence_number=self._residue_sequence_number,
        )

    def _state_atom_references(self, state_name: str) -> Tuple[AtomReference, ...]:
        return state_atom_references(
            state_name=state_name,
            states=self._states,
            state_atom_reference_lists=self._state_atom_reference_lists,
            state_atom_identifier_cache=self._state_atom_identifier_cache,
            atom_reference_from_topology_atom=self._atom_reference_from_topology_atom,
        )

    def _state_particle_masses(self, state_name: str) -> Tuple[unit.Quantity, ...]:
        return state_particle_masses(
            state_name=state_name,
            modellers=self._modellers,
            state_particle_mass_lists=self._state_particle_mass_lists,
        )

    def _state_canonical_label_index_map(self, state_name: str) -> Mapping[AtomReference, int]:
        return state_canonical_label_index_map(
            state_name=state_name,
            state_canonical_index_maps=self._state_canonical_index_maps,
            state_atom_references=self._state_atom_references,
            canonical_atom_label=self._canonical_atom_label,
        )

    def _particle_indices_from_canonical_key(
        self, state_name: str, canonical_key: CanonicalAtomKey
    ) -> Tuple[int, ...]:
        return particle_indices_from_canonical_key(
            state_name=state_name,
            canonical_key=canonical_key,
            state_canonical_label_index_map=self._state_canonical_label_index_map,
        )

    def _particle_indices_from_canonical_pair(
        self, state_name: str, canonical_pair: CanonicalPairKey
    ) -> Tuple[int, int]:
        return particle_indices_from_canonical_pair(
            state_name=state_name,
            canonical_pair=canonical_pair,
            state_canonical_label_index_map=self._state_canonical_label_index_map,
        )

    def _resolve_lambda_parameters(self, overrides: Mapping[str, str] | None) -> Mapping[str, str]:
        mapping = resolve_lambda_parameters(
            self._state_lambda_parameter_names,
            state_names=self.state_names,
            overrides=overrides,
        )
        self._state_lambda_parameter_names = mapping
        return mapping

    def _coerce_lambda_sequence(self, lambda_values: Iterable[float] | float) -> Tuple[float, ...]:
        """
        Validate user lambda values and append the residual state weight.

        Parameters
        ----------
        lambda_values
            Iterable containing ``N-1`` values ordered as :attr:`state_names`,
            or a single float for the common two-state case.

        Returns
        -------
        Tuple[float, ...]
            Sequence of length ``N`` containing lambda weights for all states.
        """

        return coerce_lambda_sequence(lambda_values, state_names=self.state_names)

    def _shared_and_state_terms(
        self, force_label: str
    ) -> Tuple[Mapping[CanonicalAtomKey, BondedTermEntry], Dict[str, Dict[CanonicalAtomKey, BondedTermEntry]]]:
        if self.bonded_force_comparisons is None:
            self.bonded_force_comparisons = self.compare_bonded_terms()
        comparison = self.bonded_force_comparisons[force_label]
        shared = comparison.universal_identical
        per_state: Dict[str, Dict[CanonicalAtomKey, BondedTermEntry]] = {}
        for category in (comparison.universal_conflicting, comparison.partial):
            for canonical_key, state_entries in category.items():
                for state_name, entry in state_entries.items():
                    per_state.setdefault(state_name, {})[canonical_key] = entry
        return shared, per_state

    @staticmethod
    def _torsion_components(parameters: Tuple[object, ...]) -> Tuple[Tuple[int, float, float], ...]:
        if not parameters:
            return tuple()
        first = parameters[0]
        if isinstance(first, tuple) and len(first) == 3:
            return parameters
        if len(parameters) == 3:
            periodicity, phase, k = parameters  # type: ignore[misc]
            return ((int(periodicity), float(phase), float(k)),)
        raise ValueError("Periodic torsion parameters must be sequences of (periodicity, phase, k).")

    @staticmethod
    def _chain_identifier(chain: app.Topology.Chain) -> str:
        return chain_identifier(chain)

    @staticmethod
    def _residue_sequence_number(residue: app.Topology.Residue, state_name: str) -> int:
        return residue_sequence_number(residue, state_name)

    def _auto_parameterize_missing_residues(
        self,
        *,
        ligand_charges: Mapping[str, float] | None,
        ligand_sdf_files: Mapping[str, object] | None,
        ligand_smiles: Mapping[str, str] | None,
        smirnoff_forcefield: str,
        charge_method: str,
        output_directory: PathLike | None,
    ) -> None:
        if not self._missing_residues:
            return
        if ligand_charges is None:
            missing = ", ".join(
                sorted({res for residues in self._missing_residues.values() for res in residues})
            )
            raise ValueError(
                "Force field lacks parameters for residues: "
                f"{missing}. Provide expected ligand total charges via 'ligand_charges' "
                "along with the corresponding ligand definitions to enable automatic parameterization."
            )
        self.parameterize_missing_residues(
            ligand_charges=ligand_charges,
            ligand_sdf_files=ligand_sdf_files,
            ligand_smiles=ligand_smiles,
            smirnoff_forcefield=smirnoff_forcefield,
            charge_method=charge_method,
            output_directory=output_directory,
        )

    def _emit(self, message: str) -> None:
        if not self._verbose:
            return
        print(message)

    def build_evb_system(self, *, initial_state: str) -> openmm.System:
        """
        Construct the OpenMM scaffold that will host all EVB forces.

        Parameters
        ----------
        initial_state :
            Name of the EVB state whose coordinates should seed the umbrella system.
            This determines the particle ordering and the reference positions stored
            in :attr:`evb_initial_positions`.

        Returns
        -------
        openmm.System
            Fresh System instance containing particles/masses and constraints that
            are common to all EVB states, but no forces. EVB-specific forces can be
            attached later without re-reading topology data.
        """

        if initial_state not in self._states:
            raise ValueError(
                f"Unknown EVB state '{initial_state}'. Available states: {self.state_names}."
            )

        system = build_evb_system_scaffold(
            initial_state=initial_state,
            state_particle_masses=self._state_particle_masses,
            add_shared_constraints_to_evb_system=lambda scaffold, state: self._add_shared_constraints_to_evb_system(
                scaffold,
                initial_state=state,
            ),
        )

        self._evb_system_scaffold = system
        self._evb_system_source_state = initial_state
        self._evb_initial_positions = self._states[initial_state].positions
        self._scaas_runtime = None
        self._scaas_backend_controller = None
        # Guard state for attach_scaas: distinguish "never enabled" from "enabled, then silently
        # invalidated by a force rebuild" so a droplet run cannot start with SCAAS quietly inactive.
        self._scaas_enable_called = False
        self._scaas_invalidated_by_rebuild = False
        if self._boundary_mode == "periodic":
            self._apply_default_box_to_scaffold(state_name=initial_state)
        return system

    def _apply_default_box_to_scaffold(self, *, state_name: str) -> None:
        """Copy a solvated state's periodic box vectors onto the EVB scaffold System.

        ``app.Simulation`` builds its ``Context`` from the System's *default* periodic
        box vectors, not from the topology's. The scaffold produced by
        :func:`build_evb_system_scaffold` carries none, so periodic (PME) electrostatics
        would not engage without this step. Called automatically by
        :meth:`build_evb_system` while ``_boundary_mode == "periodic"``.
        """

        if self._evb_system_scaffold is None:
            raise ValueError("Build the EVB scaffold before applying box vectors.")
        topology = self._states[state_name].topology
        vectors = topology.getPeriodicBoxVectors()
        if vectors is None:
            raise ValueError(
                f"State '{state_name}' has no periodic box vectors; periodic mode "
                "requires a solvated periodic box (call 'solvate_periodic_box_states')."
            )
        self._evb_system_scaffold.setDefaultPeriodicBoxVectors(*vectors)

    @staticmethod
    def _min_image_width_nm(vectors) -> float:
        """Minimum-image width of a (possibly triclinic) periodic cell, in nm.

        OpenMM's minimum-image cutoff limit is half the smallest *perpendicular* distance
        between opposite cell faces, not the shortest box-vector length. For a truncated
        octahedron or rhombic dodecahedron the box vectors are non-orthogonal, so the naive
        diagonal ``vectors[i][i]`` overestimates the usable cutoff and would wave through a
        box OpenMM later rejects deep in ``Context`` creation. The perpendicular widths are
        ``V / |b x c|`` (and cyclic), where ``V`` is the cell volume; for an orthorhombic
        box these reduce to the edge lengths, so the cube path is unchanged."""
        a = [float(vectors[0][i].value_in_unit(unit.nanometer)) for i in range(3)]
        b = [float(vectors[1][i].value_in_unit(unit.nanometer)) for i in range(3)]
        c = [float(vectors[2][i].value_in_unit(unit.nanometer)) for i in range(3)]

        def cross(u, v):
            return [u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0]]

        def norm(u):
            return math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])

        bc, ca, ab = cross(b, c), cross(c, a), cross(a, b)
        volume = abs(a[0] * bc[0] + a[1] * bc[1] + a[2] * bc[2])
        return min(volume / norm(bc), volume / norm(ca), volume / norm(ab))

    def _validate_periodic_cutoff(self, cutoff_nm: float) -> None:
        """Ensure the scaffold's box is large enough for the minimum-image cutoff.

        OpenMM requires the cutoff to be at most half the cell's minimum-image width;
        otherwise ``Context`` creation throws deep inside the platform. For non-orthorhombic
        boxes (truncated octahedron / rhombic dodecahedron) that width is the perpendicular
        face distance, not the box-vector length -- see :meth:`_min_image_width_nm`.
        """

        if self._evb_system_scaffold is None:
            raise ValueError("Build the EVB scaffold before validating periodic cutoffs.")
        vectors = self._evb_system_scaffold.getDefaultPeriodicBoxVectors()
        if vectors is None:
            raise RuntimeError(
                "Periodic mode requires box vectors on the scaffold; "
                "call 'solvate_periodic_box_states' first."
            )
        min_width = self._min_image_width_nm(vectors)
        if min_width < 2.0 * cutoff_nm:
            raise ValueError(
                f"Periodic box minimum-image width {min_width:.3f} nm is smaller than twice "
                f"the nonbonded cutoff (2 x {cutoff_nm:.3f} nm). Increase padding/box_size, "
                "choose a larger box shape, or reduce the cutoff."
            )

    @staticmethod
    def _constraint_distance_nm(value: unit.Quantity | float) -> float:
        return constraint_distance_nm(value)

    def _state_constraint_map(self, *, state_name: str, system: openmm.System) -> Dict[CanonicalPairKey, float]:
        return state_constraint_map(
            state_name=state_name,
            system=system,
            canonical_pair_key_from_indices=lambda name, atom1, atom2: self._canonical_pair_key(
                self._canonical_atom_label(name, atom1),
                self._canonical_atom_label(name, atom2),
            ),
            constraint_distance_nm=self._constraint_distance_nm,
        )

    def _add_shared_constraints_to_evb_system(
        self,
        system: openmm.System,
        *,
        initial_state: str,
        tolerance_nm: float = 1e-6,
    ) -> int:
        return add_shared_constraints_to_evb_system(
            system=system,
            initial_state=initial_state,
            state_names=self.state_names,
            reference_systems=self._get_state_systems(),
            state_constraint_map=lambda state_name, state_system: self._state_constraint_map(
                state_name=state_name,
                system=state_system,
            ),
            particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
            tolerance_nm=tolerance_nm,
        )

    def add_vdw_forces(
        self,
        *,
        lambda_parameter_names: Mapping[str, str] | None = None,
        absolute_tolerance: float = 1e-6,
        mm_cutoff_distance: unit.Quantity | float = 2.0,
        qm_cutoff_distance: unit.Quantity | float = 9.9,
    ) -> Mapping[str, openmm.Force]:
        """
        Attach shared and state-specific nonbonded interactions to the EVB scaffold.

        Notes
        -----
        Nonbonded interactions are evaluated with two pair channels:
        ``MM-MM`` pairs are computed with ``CutoffNonPeriodic`` at ``mm_cutoff_distance``,
        while ``QM-*`` pairs (QM-QM and QM-MM) are computed with ``CutoffNonPeriodic``
        at ``qm_cutoff_distance``.
        The QM region is defined through :meth:`set_qm_atoms`.

        Both cutoffs default to 2.0 nm. Keeping them equal lets OpenMM's CUDA
        platform share a single neighbor list across every CustomNonbondedForce
        on the EVB scaffold (a small but free speed win) and is a prerequisite
        for collapsing the MM/QM channel masks in the per-state overlays.
        Shortening ``mm_cutoff_distance`` below ``qm_cutoff_distance`` is
        allowed but loses the shared neighbor list; shortening
        ``qm_cutoff_distance`` is the dangerous direction for EVB profiles
        because charged reactive atoms see their long-range water Coulomb
        differently in each diabatic state and the truncation no longer
        cancels in the state-1-minus-state-0 gap.

        All state-specific CustomNonbondedForce objects are multiplied by their state's
        lambda parameter; callers must ensure the chosen EVB protocol drives those
        lambdas to physically meaningful values (typically forming a convex
        combination whose weights sum to one). Deviations from that constraint will
        scale the nonbonded energy unphysically.

        Van der Waals and electrostatics partitions are cached per ``build_state_systems`` signature.
        Rebuilding OpenMM state systems with new keywords (e.g., when periodic cutoffs
        are introduced) requires invalidating and recreating the EVBSystem to avoid
        stale cached partitions.

        Caution
        -------
        Atoms participating in any state-specific exception are removed from the shared
        NonbondedForce and handled via CustomNonbondedForce overlays; extensive reactive
        regions may therefore shift most nonbonded work to the slower custom kernels.
        In addition, the present exclusion bookkeeping assumes zeroed pairs are consistent
        across states—if that assumption is violated, double-counting checks should be
        added before enabling cutoffs or PME-style protocols.
        """

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("Call 'build_evb_system' before adding van der Waals forces.")
        if (
            self._shared_nonbonded_force is not None
            or self._shared_nonbonded_qm_force is not None
            or self._shared_nonbonded_exception_force is not None
            or self._shared_nonbonded_qm_exception_force is not None
            or self._state_specific_nonbonded_forces
            or self._state_specific_nonbonded_qm_forces
            or self._state_specific_vdw_exception_forces
            or self._state_specific_vdw_qm_exception_forces
            or self._state_specific_electrostatics_exception_forces
            or self._state_specific_electrostatics_qm_exception_forces
        ):
            raise RuntimeError("Nonbonded forces have already been added to this EVB system.")

        lambda_map = self._resolve_lambda_parameters(lambda_parameter_names)
        partition = self._vdw_partition_terms(tolerance=absolute_tolerance)
        electrostatics_partition = self._electrostatics_partition_terms(tolerance=absolute_tolerance)
        mm_cutoff_q = _units.coerce_length_nm(mm_cutoff_distance, label="mm_cutoff_distance")
        mm_cutoff_nm = float(mm_cutoff_q.value_in_unit(unit.nanometer))
        qm_cutoff_q = _units.coerce_length_nm(qm_cutoff_distance, label="qm_cutoff_distance")
        qm_cutoff_nm = float(qm_cutoff_q.value_in_unit(unit.nanometer))
        self._mm_nonbonded_cutoff_nm = mm_cutoff_nm
        self._qm_nonbonded_cutoff_nm = qm_cutoff_nm
        # When MM and QM cutoffs match, the MM/QM channel masks are
        # complementary and we can collapse the per-channel duplicate
        # CustomNonbondedForce instances into single forces. OpenMM also
        # shares one neighbor list across same-cutoff forces in this case.
        cutoffs_unified = abs(mm_cutoff_nm - qm_cutoff_nm) < 1e-12
        periodic = self._boundary_mode == "periodic"
        if periodic:
            self._validate_periodic_cutoff(max(mm_cutoff_nm, qm_cutoff_nm))
        evb_system = self._evb_system_scaffold
        initial_state = self._evb_system_source_state
        reference_systems = self._get_state_systems()
        all_exception_parameters = collect_exception_parameters_across_states(
            reference_systems=reference_systems,
            find_nonbonded_force=self._find_nonbonded_force,
            extract_vdw_exception_parameters=self._extract_vdw_exception_parameters,
        )
        reference_labels = tuple(
            self._canonical_atom_label(initial_state, idx) for idx in range(evb_system.getNumParticles())
        )
        reference_label_set = set(reference_labels)
        unknown_qm_labels = sorted(self._qm_region_atoms - reference_label_set)
        # Guard the reactive-atom cutoff HERE, where the QM labels are resolved to indices.
        self._assert_qm_cutoff_covers_system(
            qm_cutoff_nm,
            [i for i, lab in enumerate(reference_labels) if lab in self._qm_region_atoms],
        )
        if unknown_qm_labels:
            raise ValueError(
                "Configured QM atom labels are not present in the active EVB scaffold. "
                "Reapply 'set_qm_atoms(...)' after topology/mapping changes."
            )
        qm_labels = set(self._qm_region_atoms)
        has_qm_region = bool(qm_labels)

        ordered_global_exclusions = collect_global_exclusion_pairs(
            partition=partition,
            electrostatics_partition=electrostatics_partition,
        )
        global_exclusion_pair_set = set(ordered_global_exclusions)

        # In periodic mode electrostatics is handled by dedicated per-state PME
        # NonbondedForce objects, so the shared custom force carries van der Waals
        # only. The non-periodic expressions are kept byte-for-byte to preserve the
        # default code path.
        _sig = lj_sigma_combining_expr(self._lj_combining)
        _lj = f"4*sqrt(epsilon1*epsilon2) * (({_sig}/r)^12 - ({_sig}/r)^6)"
        _coul = " + ONE_4PI_EPS0 * charge1 * charge2 / r"
        shared_body = _lj if periodic else _lj + _coul
        if periodic:
            shared_expression_merged = _lj
        else:
            shared_expression_merged = _lj + _coul
        shared_expression_mm = f"(1 - qm1)*(1 - qm2) * ({shared_body})"
        shared_expression_qm = f"(1 - (1 - qm1)*(1 - qm2)) * ({shared_body})"
        if cutoffs_unified:
            # Single shared force; the MM and QM channels are complementary
            # at equal cutoffs, so dropping the channel mask is equivalent to
            # summing the two old shared forces.
            shared_force = build_channel_force(
                expression=shared_expression_merged,
                cutoff_distance=mm_cutoff_q,
                one_4pi_eps0=ONE_4PI_EPS0,
                periodic=periodic,
            )
            shared_qm_force = None
        else:
            shared_force = build_channel_force(
                expression=shared_expression_mm,
                cutoff_distance=mm_cutoff_q,
                one_4pi_eps0=ONE_4PI_EPS0,
                periodic=periodic,
            )
            shared_qm_force = (
                build_channel_force(
                    expression=shared_expression_qm,
                    cutoff_distance=qm_cutoff_q,
                    one_4pi_eps0=ONE_4PI_EPS0,
                    periodic=periodic,
                )
                if has_qm_region
                else None
            )

        qm_indices: set[int] = set()
        for index, label in enumerate(reference_labels):
            sigma, epsilon = partition.universal_atoms.get(label, (0.0, 0.0))
            charge = electrostatics_partition.universal_atoms.get(label, 0.0)
            qm_flag = 1.0 if label in qm_labels else 0.0
            if qm_flag > 0.5:
                qm_indices.add(index)
            shared_force.addParticle((sigma, epsilon, charge, qm_flag))
            if shared_qm_force is not None:
                shared_qm_force.addParticle((sigma, epsilon, charge, qm_flag))

        for pair in ordered_global_exclusions:
            idx_a, idx_b = self._particle_indices_from_canonical_pair(initial_state, pair)
            shared_force.addExclusion(idx_a, idx_b)
            if shared_qm_force is not None:
                shared_qm_force.addExclusion(idx_a, idx_b)

        if shared_qm_force is not None:
            restrict_to_subset_interaction_group(shared_qm_force, qm_indices)

        evb_system.addForce(shared_force)
        self._register_force_to_group(shared_force, term="shared", state=None)
        self._shared_nonbonded_force = shared_force

        added_forces: Dict[str, openmm.Force] = {
            "shared_mm": shared_force,
            "shared": shared_force,
        }
        if shared_qm_force is not None:
            evb_system.addForce(shared_qm_force)
            self._register_force_to_group(shared_qm_force, term="shared_qm", state=None)
            added_forces["shared_qm"] = shared_qm_force
        self._shared_nonbonded_qm_force = shared_qm_force

        mm_shared_exception_terms, qm_shared_exception_terms = split_shared_exception_terms(
            partition=partition,
            electrostatics_partition=electrostatics_partition,
            all_exception_parameters=all_exception_parameters,
            qm_labels=qm_labels,
            pair_in_qm_channel=lambda pair: pair_in_qm_channel(pair, qm_labels=qm_labels),
        )

        shared_exception_force = build_shared_exception_force(
            terms=mm_shared_exception_terms,
            one_4pi_eps0=ONE_4PI_EPS0,
            initial_state=initial_state,
            particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
            periodic=periodic,
        )
        if shared_exception_force is not None:
            evb_system.addForce(shared_exception_force)
            self._register_force_to_group(shared_exception_force, term="shared", state=None)
            added_forces["shared_exceptions_mm"] = shared_exception_force
        self._shared_nonbonded_exception_force = shared_exception_force

        shared_qm_exception_force = build_shared_exception_force(
            terms=qm_shared_exception_terms,
            one_4pi_eps0=ONE_4PI_EPS0,
            initial_state=initial_state,
            particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
            periodic=periodic,
        )
        if shared_qm_exception_force is not None:
            evb_system.addForce(shared_qm_exception_force)
            self._register_force_to_group(shared_qm_exception_force, term="shared", state=None)
            added_forces["shared_exceptions_qm"] = shared_qm_exception_force
        self._shared_nonbonded_qm_exception_force = shared_qm_exception_force

        state_forces: Dict[str, openmm.CustomNonbondedForce] = {}
        state_qm_forces: Dict[str, openmm.CustomNonbondedForce] = {}
        exception_forces: Dict[str, openmm.CustomBondForce] = {}
        exception_qm_forces: Dict[str, openmm.CustomBondForce] = {}
        electro_exception_forces: Dict[str, openmm.CustomBondForce] = {}
        electro_exception_qm_forces: Dict[str, openmm.CustomBondForce] = {}

        for state_name in self.state_names:
            atom_parameters = partition.state_specific_atoms.get(state_name, {})
            exception_parameters = partition.state_specific_exceptions.get(state_name, {})
            electro_atom_parameters = electrostatics_partition.state_specific_atoms.get(state_name, {})
            electro_exception_parameters = electrostatics_partition.state_specific_exceptions.get(state_name, {})
            state_exclusion_pairs = collect_state_exclusion_pairs(
                partition=partition,
                electrostatics_partition=electrostatics_partition,
                state_name=state_name,
                exception_parameters=exception_parameters,
                electro_exception_parameters=electro_exception_parameters,
            )
            restored_pairs = tuple(sorted(global_exclusion_pair_set - state_exclusion_pairs))
            lambda_param = lambda_map[state_name]

            if atom_parameters or electro_atom_parameters:
                if cutoffs_unified:
                    force = build_state_combined_force(
                        lambda_param=lambda_param,
                        cutoff_distance=mm_cutoff_q,
                        qm_channel=None,
                        one_4pi_eps0=ONE_4PI_EPS0,
                        periodic=periodic,
                        include_coulomb=not periodic,
                        lj_combining=self._lj_combining,
                    )
                    force_qm = None
                else:
                    force = build_state_combined_force(
                        lambda_param=lambda_param,
                        cutoff_distance=mm_cutoff_q,
                        qm_channel=False,
                        one_4pi_eps0=ONE_4PI_EPS0,
                        periodic=periodic,
                        include_coulomb=not periodic,
                        lj_combining=self._lj_combining,
                    )
                    force_qm = (
                        build_state_combined_force(
                            lambda_param=lambda_param,
                            cutoff_distance=qm_cutoff_q,
                            qm_channel=True,
                            one_4pi_eps0=ONE_4PI_EPS0,
                            periodic=periodic,
                            include_coulomb=not periodic,
                        lj_combining=self._lj_combining,
                        )
                        if has_qm_region
                        else None
                    )
                active_indices = populate_state_combined_particles(
                    force=force,
                    force_qm=force_qm,
                    reference_labels=reference_labels,
                    vdw_atom_parameters=atom_parameters,
                    vdw_universal_atoms=partition.universal_atoms,
                    electro_atom_parameters=electro_atom_parameters,
                    electro_universal_atoms=electrostatics_partition.universal_atoms,
                    qm_labels=qm_labels,
                )

                for pair in ordered_global_exclusions:
                    idx_a, idx_b = self._particle_indices_from_canonical_pair(initial_state, pair)
                    force.addExclusion(idx_a, idx_b)
                    if force_qm is not None:
                        force_qm.addExclusion(idx_a, idx_b)

                restrict_to_subset_interaction_group(force, active_indices)
                if force_qm is not None:
                    restrict_to_subset_interaction_group(force_qm, active_indices)

                evb_system.addForce(force)
                self._register_force_to_group(force, term="nonbonded", state=state_name)
                state_forces[state_name] = force
                added_forces[state_name] = force

                if force_qm is not None:
                    evb_system.addForce(force_qm)
                    self._register_force_to_group(force_qm, term="nonbonded_qm", state=state_name)
                    state_qm_forces[state_name] = force_qm
                    added_forces[f"{state_name}_qm"] = force_qm

            vdw_bond_terms = build_vdw_bond_terms(
                exception_parameters=exception_parameters,
                restored_pairs=restored_pairs,
                atom_parameters=atom_parameters,
                universal_atoms=partition.universal_atoms,
            )
            vdw_bond_terms_mm, vdw_bond_terms_qm = split_terms_by_qm_channel(
                vdw_bond_terms,
                pair_in_qm_channel=lambda pair: pair_in_qm_channel(pair, qm_labels=qm_labels),
            )

            bond_force = build_lambda_scaled_vdw_exception_force(
                lambda_param=lambda_param,
                terms=vdw_bond_terms_mm,
                initial_state=initial_state,
                particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
                periodic=periodic,
            )
            if bond_force is not None:
                evb_system.addForce(bond_force)
                self._register_force_to_group(bond_force, term="nonbonded", state=state_name)
                exception_forces[state_name] = bond_force
                added_forces[f"{state_name}_exceptions"] = bond_force

            bond_qm_force = build_lambda_scaled_vdw_exception_force(
                lambda_param=lambda_param,
                terms=vdw_bond_terms_qm,
                initial_state=initial_state,
                particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
                periodic=periodic,
            )
            if bond_qm_force is not None:
                evb_system.addForce(bond_qm_force)
                self._register_force_to_group(bond_qm_force, term="nonbonded_qm", state=state_name)
                exception_qm_forces[state_name] = bond_qm_force
                added_forces[f"{state_name}_exceptions_qm"] = bond_qm_force

            electro_bond_terms = build_electro_bond_terms(
                electro_exception_parameters=electro_exception_parameters,
                restored_pairs=restored_pairs,
                atom_parameters=electro_atom_parameters,
                universal_atoms=electrostatics_partition.universal_atoms,
            )
            electro_bond_terms_mm, electro_bond_terms_qm = split_terms_by_qm_channel(
                electro_bond_terms,
                pair_in_qm_channel=lambda pair: pair_in_qm_channel(pair, qm_labels=qm_labels),
            )

            elec_bond_force = build_lambda_scaled_electro_exception_force(
                lambda_param=lambda_param,
                terms=electro_bond_terms_mm,
                initial_state=initial_state,
                particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
                one_4pi_eps0=ONE_4PI_EPS0,
                periodic=periodic,
            )
            if elec_bond_force is not None:
                evb_system.addForce(elec_bond_force)
                self._register_force_to_group(elec_bond_force, term="nonbonded", state=state_name)
                electro_exception_forces[state_name] = elec_bond_force
                added_forces[f"{state_name}_electrostatics_exceptions"] = elec_bond_force

            elec_bond_qm_force = build_lambda_scaled_electro_exception_force(
                lambda_param=lambda_param,
                terms=electro_bond_terms_qm,
                initial_state=initial_state,
                particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
                one_4pi_eps0=ONE_4PI_EPS0,
                periodic=periodic,
            )
            if elec_bond_qm_force is not None:
                evb_system.addForce(elec_bond_qm_force)
                self._register_force_to_group(elec_bond_qm_force, term="nonbonded_qm", state=state_name)
                electro_exception_qm_forces[state_name] = elec_bond_qm_force
                added_forces[f"{state_name}_electrostatics_exceptions_qm"] = elec_bond_qm_force

        if periodic:
            self._add_periodic_pme_forces(
                evb_system=evb_system,
                reference_systems=reference_systems,
                reference_labels=reference_labels,
                ordered_global_exclusions=ordered_global_exclusions,
                cutoff_distance=mm_cutoff_q,
                lambda_map=lambda_map,
                added_forces=added_forces,
            )

        self._state_specific_nonbonded_forces = state_forces
        self._state_specific_nonbonded_qm_forces = state_qm_forces
        self._state_specific_vdw_exception_forces = exception_forces
        self._state_specific_vdw_qm_exception_forces = exception_qm_forces
        self._state_specific_electrostatics_exception_forces = electro_exception_forces
        self._state_specific_electrostatics_qm_exception_forces = electro_exception_qm_forces
        return MappingProxyType(added_forces)

    def _add_periodic_pme_forces(
        self,
        *,
        evb_system: openmm.System,
        reference_systems: Mapping[str, openmm.System],
        reference_labels: Sequence[AtomReference],
        ordered_global_exclusions: Sequence[CanonicalPairKey],
        cutoff_distance: unit.Quantity,
        lambda_map: Mapping[str, str],
        added_forces: Dict[str, openmm.Force],
    ) -> None:
        """Attach one Coulomb-only PME ``NonbondedForce`` per diabatic state.

        Each force carries that state's *complete* charge vector (read from the
        reference state System) and is scaled by ``g_state = sqrt(lambda_state)``
        through :func:`build_state_pme_force`, so its full Ewald energy contributes
        ``lambda_state * E_elec(q_state)`` to the EVB diagonal. Intramolecular
        exception/exclusion pairs are excluded here and their direct Coulomb is
        restored, lambda-scaled, by the electrostatics exception bond forces;
        reactive soft-pair exclusions are added later by the soft-pair installer
        (which already patches every ``NonbondedForce`` in the system).
        """

        initial_state = self._evb_system_source_state
        num_particles = evb_system.getNumParticles()
        self._state_specific_pme_forces = {}
        self._state_pme_g_param_names = {}

        for state_name in self.state_names:
            state_system = reference_systems[state_name]
            state_nb = self._find_nonbonded_force(state_system)
            label_to_state_index = self._state_canonical_label_index_map(state_name)
            state_charges: list[float] = []
            for index in range(num_particles):
                label = reference_labels[index]
                state_index = label_to_state_index[label]
                charge, _sigma, _epsilon = state_nb.getParticleParameters(state_index)
                state_charges.append(float(charge.value_in_unit(unit.elementary_charge)))

            g_param = f"evb_pme_g_{state_name}"
            lambda_param = lambda_map[state_name]
            pme_force = build_state_pme_force(
                state_charges=state_charges,
                g_param=g_param,
                lambda_param=lambda_param,
                cutoff_distance=cutoff_distance,
                ewald_error_tolerance=self._ewald_error_tolerance,
            )
            for pair in ordered_global_exclusions:
                idx_a, idx_b = self._particle_indices_from_canonical_pair(initial_state, pair)
                pme_force.addException(idx_a, idx_b, 0.0, 1.0, 0.0)

            evb_system.addForce(pme_force)
            self._register_force_to_group(pme_force, term="nonbonded_pme", state=state_name)
            self._state_specific_pme_forces[state_name] = pme_force
            self._state_pme_g_param_names[state_name] = g_param
            added_forces[f"{state_name}_pme"] = pme_force

    def get_forming_breaking_bond_pairs(
        self,
        *,
        reference_state: str | None = None,
        comparison_states: Iterable[str] | None = None,
        prefer_index_pairs: bool = True,
    ) -> Mapping[CanonicalPairKey, Mapping[str, BondedTermEntry]]:
        """
        Return atom pairs that form/break across EVB states.

        Notes
        -----
        By default this method first compares bond *presence* by atom-index pairs
        between ``reference_state`` and ``comparison_states``. This preserves a
        directional EVB view (useful for symmetric reactions where canonical atom
        mapping can hide bond swaps).

        If no index-based differences are detected (or index comparison is disabled),
        it falls back to the canonical bonded-term ``partial`` category from
        :meth:`compare_bonded_terms`.
        """

        if reference_state is None:
            reference = self.state_names[0]
        else:
            reference = str(reference_state)
        if reference not in self._states:
            raise ValueError(
                f"State '{reference}' not found. Available states: {', '.join(self.state_names)}."
            )

        if comparison_states is None:
            selected_states = tuple(state for state in self.state_names if state != reference)
        else:
            normalized: list[str] = []
            seen: Set[str] = set()
            for state in comparison_states:
                state_name = str(state)
                if state_name == reference:
                    continue
                if state_name not in self._states:
                    raise ValueError(
                        f"State '{state_name}' not found. Available states: {', '.join(self.state_names)}."
                    )
                if state_name in seen:
                    continue
                seen.add(state_name)
                normalized.append(state_name)
            selected_states = tuple(normalized)
        if not selected_states:
            raise ValueError("At least one comparison state is required to detect forming/breaking bonds.")

        if prefer_index_pairs:
            directional_pairs = self._forming_breaking_bond_pairs_by_index(
                reference_state=reference,
                comparison_states=selected_states,
            )
            if directional_pairs:
                return directional_pairs

        if self.bonded_force_comparisons is None:
            self.bonded_force_comparisons = self.compare_bonded_terms()
        comparison = self.bonded_force_comparisons["bonds"]
        pairs: Dict[CanonicalPairKey, Mapping[str, BondedTermEntry]] = {}
        for canonical_key, state_entries in comparison.partial.items():
            if len(canonical_key) != 2:
                continue
            pair = self._canonical_pair_key(canonical_key[0], canonical_key[1])
            pairs[pair] = MappingProxyType(dict(state_entries))
        return MappingProxyType(pairs)

    def write_vmd_session(
        self,
        output_tcl,
        topology_file,
        dcd_file,
        *,
        reference_state: str | None = None,
        sphere_center_nm=None,
        sphere_radius_a: float | None = None,
        **kwargs,
    ) -> str:
        """Write a VMD TCL session for a saved EVB trajectory.

        Thin convenience wrapper over :func:`evbopenmm.analysis.vmd.write_vmd_session`: the
        reactive-core atoms and forming/breaking bond pairs are filled in automatically from
        :meth:`get_forming_breaking_bond_pairs` (resolved to topology indices). The SCAAS
        sphere centre/radius are NOT persisted on the system after solvation, so the caller
        supplies them (via ``sphere_center_nm``/``sphere_radius_a``) when ``show_sphere=True``.
        Every visualisation feature is a keyword forwarded to the free function, so callers
        toggle each one individually.
        """
        from .analysis.vmd import write_vmd_session as _write

        ref = reference_state or self.state_names[0]
        topology = self.states[ref].topology
        cidx = {
            (a.residue.chain.id, str(a.residue.id), a.name): a.index
            for a in topology.atoms()
        }

        def _resolve(refkey):
            if isinstance(refkey, int):
                return int(refkey)
            _state, chain, resid, name = refkey
            return cidx[(chain, str(resid), name)]

        pairs = []
        for key in self.get_forming_breaking_bond_pairs(reference_state=ref):
            a, b = key
            pairs.append((_resolve(a), _resolve(b)))

        idx_to_name = {a.index: a.name for a in topology.atoms()}
        reactive = {
            idx_to_name.get(i, f"idx{i}"): i
            for i in sorted({x for p in pairs for x in p})
        }

        return _write(
            output_tcl,
            topology_file,
            dcd_file,
            reactive_atom_indices=reactive,
            forming_breaking_pairs=pairs,
            sphere_center_nm=sphere_center_nm,
            sphere_radius_a=sphere_radius_a,
            **kwargs,
        )

    def _harmonic_bond_terms_by_index(
        self,
        state_name: str,
    ) -> Dict[Tuple[int, int], BondedTermEntry]:
        systems = self._get_state_systems()
        return extract_harmonic_bond_terms_by_index(systems[state_name])

    def _forming_breaking_bond_pairs_by_index(
        self,
        *,
        reference_state: str,
        comparison_states: Sequence[str],
    ) -> Mapping[CanonicalPairKey, Mapping[str, BondedTermEntry]]:
        state_group = (reference_state, *comparison_states)
        state_atom_references = {state_name: self._state_atom_references(state_name) for state_name in state_group}
        terms_by_state: Dict[str, Dict[Tuple[int, int], BondedTermEntry]] = {
            state_name: self._harmonic_bond_terms_by_index(state_name) for state_name in state_group
        }
        return collect_forming_breaking_bond_pairs_by_index(
            reference_state=reference_state,
            comparison_states=comparison_states,
            state_atom_references=state_atom_references,
            state_bond_terms_by_index=terms_by_state,
            atom_reference_from_index=self._atom_reference_from_index,
            # Map each raw atom reference to its cross-state canonical label
            # (representative) before forming the pair key, so the keys match the
            # canonical-label index maps used downstream. Other code paths apply
            # this; the by-index path (taken when states share an atom layout,
            # e.g. an intramolecular reactant/product) previously did not.
            canonical_pair_key=lambda a, b: self._canonical_pair_key(
                self._atom_canonical_labels.get(a, a),
                self._atom_canonical_labels.get(b, b),
            ),
        )

    def _bonded_terms_by_pair(self) -> Mapping[CanonicalPairKey, Mapping[str, BondedTermEntry]]:
        """
        Return harmonic bond terms grouped by canonical atom pair across states.
        """

        systems = self._get_state_systems()
        by_force = self._collect_bonded_terms_by_state(systems)
        bonds_by_state = by_force.get("bonds", {})
        return group_bonded_terms_by_pair(
            bonds_by_state=bonds_by_state,
            canonical_pair_key=self._canonical_pair_key,
        )

    def _normalize_reactive_pair(self, pair: CanonicalPairKey) -> CanonicalPairKey:
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise TypeError(
                "Reactive pairs must be 2-tuples of atom references "
                "((state, chain, resid, atom), (state, chain, resid, atom))."
            )
        left = self._normalize_atom_reference(pair[0])
        right = self._normalize_atom_reference(pair[1])
        left_label = self._atom_canonical_labels.get(left, left)
        right_label = self._atom_canonical_labels.get(right, right)
        return self._canonical_pair_key(left_label, right_label)

    def _resolve_reactive_morse_parameters(
        self,
        reactive_morse_parameters: (
            "MorseBondParameter | Mapping[CanonicalPairKey, MorseBondParameter]"
        ),
        reactive_pairs: Sequence[CanonicalPairKey],
    ) -> Dict[CanonicalPairKey, MorseBondParameter]:
        """Normalise a single MorseBondParameter or per-pair mapping into a dict."""
        if isinstance(reactive_morse_parameters, MorseBondParameter):
            return {pair: reactive_morse_parameters for pair in reactive_pairs}
        if isinstance(reactive_morse_parameters, Mapping):
            normalized: Dict[CanonicalPairKey, MorseBondParameter] = {}
            for pair, parameter in reactive_morse_parameters.items():
                normalized[self._normalize_reactive_pair(pair)] = parameter
            missing = [pair for pair in reactive_pairs if pair not in normalized]
            if missing:
                raise ValueError(
                    "Missing Morse parameters for reactive pairs: "
                    f"{missing}. Provide a MorseBondParameter for every reactive "
                    "pair or pass a single MorseBondParameter to apply to all."
                )
            return normalized
        raise TypeError(
            "reactive_morse_parameters must be a MorseBondParameter or "
            "Mapping[CanonicalPairKey, MorseBondParameter]; got "
            f"{type(reactive_morse_parameters).__name__}"
        )

    def replace_reactive_bonds_with_morse(
        self,
        reactive_morse_parameters: (
            "MorseBondParameter | Mapping[CanonicalPairKey, MorseBondParameter]"
        ),
    ) -> Mapping[str, openmm.CustomBondForce]:
        """Swap the FF-derived harmonic entry on each reactive pair for Morse.

        Energy form per state, per reactive pair::

            lambda_state * de * (1 - exp(-alpha * (r - r0)))^2

        Pre-conditions: ``add_bond_forces`` has been called and produced a
        per-state harmonic ``CustomBondForce`` per state. The reactive pair's
        harmonic entry is left in the force list with ``k`` set to 0 (a no-op
        bond), since OpenMM's ``CustomBondForce`` has no ``removeBond``; the
        Morse force is added alongside.

        Returns the new per-state Morse forces keyed by state name (only states
        in which the pair is 1-2 bonded).
        """
        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("Call 'build_evb_system' before replacing reactive bonds.")
        if not self._state_specific_bond_forces:
            raise ValueError("Call 'add_bond_forces' before replacing reactive bonds.")

        forming_breaking = self.get_forming_breaking_bond_pairs()
        if not forming_breaking:
            return MappingProxyType({})
        reactive_pair_list = list(forming_breaking.keys())
        morse_by_pair = self._resolve_reactive_morse_parameters(
            reactive_morse_parameters, reactive_pair_list
        )

        # Zero out the reactive entries' harmonic spring constants in the per-state
        # bond forces. The energy form is ``lambda * 0.5 * k * (r-r0)^2``, so k=0
        # makes the bond a no-op while leaving the entry slot stable (CustomBondForce
        # has no removeBond).
        initial_state = self._evb_system_source_state
        for state_name, force in self._state_specific_bond_forces.items():
            for bond_idx in range(force.getNumBonds()):
                atom1, atom2, params = force.getBondParameters(bond_idx)
                ref_a = self._atom_reference_from_index(initial_state, atom1)
                ref_b = self._atom_reference_from_index(initial_state, atom2)
                pair = self._canonical_pair_key(
                    self._atom_canonical_labels.get(ref_a, ref_a),
                    self._atom_canonical_labels.get(ref_b, ref_b),
                )
                if pair in morse_by_pair:
                    # Preserve r0 (params[0]); zero k (params[1])
                    force.setBondParameters(bond_idx, atom1, atom2, [params[0], 0.0])

        # Build per-state Morse bond force for the reactive pairs that are
        # 1-2 bonded in that state.
        added: Dict[str, openmm.CustomBondForce] = {}
        for state_name in self.state_names:
            state_morse_entries: Dict[CanonicalPairKey, tuple[float, float, float]] = {}
            for pair, state_to_term in forming_breaking.items():
                if state_name in state_to_term:
                    parameter = morse_by_pair[pair]
                    state_morse_entries[pair] = (
                        float(parameter.r0_nm),
                        float(parameter.de_kj_mol),
                        float(parameter.alpha_per_nm),
                    )
            if not state_morse_entries:
                continue
            lambda_param = self._state_lambda_parameter_names[state_name]
            morse_force = build_reactive_morse_bond_force(
                lambda_param=lambda_param,
                state_parameters=state_morse_entries,
                initial_state=initial_state,
                particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
            )
            self._evb_system_scaffold.addForce(morse_force)
            self._register_force_to_group(morse_force, term="reactive_morse", state=state_name)
            self._state_specific_reactive_morse_bond_forces[state_name] = morse_force
            added[state_name] = morse_force

            # Endpoint-integration-only harmonic, matched to the Morse minimum
            # (same r0; curvature k = 2*de*alpha^2). Lives in a dedicated force
            # group OUTSIDE the (state, term) energy registry so the diabatic-
            # energy sampler never sees it, and is gated by g_regular_ff (default
            # 0) so it is a no-op until endpoint regular-FF integration enables it.
            harmonic_entries = {
                pair: (r0_nm, 2.0 * de_kj_mol * alpha_per_nm * alpha_per_nm)
                for pair, (r0_nm, de_kj_mol, alpha_per_nm) in state_morse_entries.items()
            }
            harmonic_force = build_reactive_harmonic_bond_force(
                lambda_param=lambda_param,
                regular_ff_param=self._regular_ff_global_param,
                state_parameters=harmonic_entries,
                initial_state=initial_state,
                particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
            )
            self._evb_system_scaffold.addForce(harmonic_force)
            if self._next_force_group_id >= 32:
                raise RuntimeError(
                    "EVBSystem ran out of OpenMM force groups while allocating an "
                    "endpoint regular-FF group (32 max)."
                )
            regular_group = self._next_force_group_id
            self._next_force_group_id += 1  # reserved, intentionally NOT in _force_group_definitions
            harmonic_force.setForceGroup(regular_group)
            self._endpoint_regular_ff_forces[state_name] = harmonic_force
            self._endpoint_regular_ff_groups[state_name] = regular_group

        return MappingProxyType(added)

    def add_reactive_softpair_vdw_forces(
        self,
        *,
        lambda_parameter_names: Mapping[str, str] | None = None,
        reactive_pairs: Iterable[CanonicalPairKey] | None = None,
        sp_c_kcal_mol_per_atom: float = 90.0,
        sp_a_per_angstrom_per_atom: float = 2.5,
        soft_pair_parameters_by_atom_type: Mapping[str, tuple[float, float]] | None = None,
        softpair_ab_by_pair: Mapping[CanonicalPairKey, tuple[float, float]] | None = None,
        exclude_from_shared_nonbonded: bool = True,
    ) -> Mapping[str, openmm.CustomBondForce]:
        """Q6-binary-accurate soft-pair for forming/breaking atom pairs.

        Installs per-state CustomBondForces implementing the exact form used
        by Q6 (``nonbonded.f90`` / ``simprep.f90`` in the Q6 repository):

        ``V(r) = lambda_state * [ A * exp(-B * r) + k_coul * q_i * q_j / r ]``

        A and B are the reactive soft-pair parameters combined geometrically:
        ``A_cross = C_i * C_j`` and ``B_cross = a_i * a_j``. Defaults follow
        Q6's conventional per-atom values (C ≈ 90 kcal/mol, a ≈ 2.5 /Å) which
        give bounded behaviour at short r without any softcore lookup. Per-state
        ``q_i * q_j`` is taken from the state's own NonbondedForce, so each
        diabatic state contributes its physical charge product.

        ``soft_pair_parameters_by_atom_type`` lets callers override the scalar
        defaults on a per-element basis (matching Q6's atom-type parameterization
        in e.g. ``probr_cl.prm``). Keyed by element symbol, values are
        ``(SP_C_kcal_mol, SP_a_per_angstrom)``. Example for Cl⁻ + CH₃Cl matching
        Q6's ``ccl_cl.fep``::

            {"C": (91.0, 2.5), "Cl": (90.0, 2.5)}

        For each reactive pair the per-atom C and a are looked up by element and
        combined geometrically (``A = C_i·C_j``, ``B = a_i·a_j``); elements absent
        from the dict fall back to ``sp_c_kcal_mol_per_atom`` /
        ``sp_a_per_angstrom_per_atom``.

        ``softpair_ab_by_pair`` gives an explicit **per-pair** override of the
        combined ``(A, B)`` in OpenMM units — ``A`` in kJ/mol, ``B`` in 1/nm —
        keyed by :class:`CanonicalPairKey`. It takes precedence over the
        per-element combining for the listed pairs (other pairs still combine by
        element). This is the faithful route for importing published Q/Åqvist EVB
        setups whose soft-repulsion is tabulated per pair and does *not* factor
        into a single per-element ``C_i`` (e.g. two atoms of the same element
        carrying different soft-pair amplitudes). Every key must also be present
        in ``reactive_pairs`` so the pair is selected. Values must be
        ``(A_kj_mol >= 0, B_per_nm > 0)``.

        When ``exclude_from_shared_nonbonded=True`` (default), the reactive
        pairs are added as exclusions in the shared scaffold nonbonded force so
        the reactive soft-pair is the *sole* nonbonded contribution for that pair —
        matching the parity-validated Aqvist/Q6 EVB form.

        This method is normally invoked automatically by ``EVBSystem.__init__``
        when ``auto_build_forces=True`` and reactive pairs exist. Callers using
        ``auto_build_forces=False`` should invoke it themselves after
        ``add_vdw_forces`` / ``add_bond_forces`` / ``add_angle_forces`` /
        ``add_torsion_forces``.
        """
        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("Call 'build_evb_system' before adding reactive soft-pair forces.")
        if self._shared_nonbonded_force is None:
            raise ValueError("Call 'add_vdw_forces' before adding reactive soft-pair forces.")
        if self._state_specific_reactive_vdw_forces:
            raise RuntimeError(
                "Reactive softpair VdW forces were already added; "
                "pick only one soft-pair scheme per scaffold."
            )

        lambda_map = self._resolve_lambda_parameters(lambda_parameter_names)
        if reactive_pairs is None:
            selected_pairs = set(self.get_forming_breaking_bond_pairs().keys())
        else:
            selected_pairs = {self._normalize_reactive_pair(pair) for pair in reactive_pairs}
        if not selected_pairs:
            raise ValueError("No reactive atom pairs selected for Q6 softpair.")

        # Explicit per-pair (A, B) overrides in OpenMM units (A kJ/mol, B 1/nm).
        explicit_ab: Dict[CanonicalPairKey, tuple[float, float]] = {}
        if softpair_ab_by_pair:
            for raw_pair, ab in softpair_ab_by_pair.items():
                norm = self._normalize_reactive_pair(raw_pair)
                if norm not in selected_pairs:
                    raise ValueError(
                        f"softpair_ab_by_pair pair {raw_pair!r} is not among the selected "
                        "reactive pairs; add it to 'reactive_pairs' so it is included."
                    )
                try:
                    a_kj, b_per_nm = float(ab[0]), float(ab[1])
                except (TypeError, ValueError, IndexError) as exc:
                    raise ValueError(
                        f"softpair_ab_by_pair[{raw_pair!r}] must be a 2-tuple "
                        f"(A_kj_mol, B_per_nm); got {ab!r}."
                    ) from exc
                if not (a_kj >= 0.0 and b_per_nm > 0.0):
                    raise ValueError(
                        f"softpair_ab_by_pair[{raw_pair!r}] must be (A>=0 kJ/mol, B>0 /nm); "
                        f"got {ab!r}."
                    )
                explicit_ab[norm] = (a_kj, b_per_nm)

        initial_state = self._evb_system_source_state
        system = self._evb_system_scaffold

        # Build per-state map of pairs bonded in that state. Q6 applies softpair
        # ONLY to reactive pairs that are NOT 1-2 bonded in the current state;
        # bonded reactive pairs are handled by Morse alone (their Coulomb is
        # implicitly part of the bond energy, not separately accounted for).
        # Without this filter, we would double-count q_i*q_j / r on the bonded
        # reactive pair, adding a spurious positive offset that shifts E1 and
        # E2 away from Q6 by ~3-4 kcal/mol each. See gas-phase parity validation
        # against Q6's Q-SUM at eq.log step 0.
        bonded_reactive_by_state: Dict[str, set] = {}
        try:
            reactive_map = self.get_forming_breaking_bond_pairs()
        except Exception:
            reactive_map = {}
        for raw_pair, state_entries in reactive_map.items():
            normalized = self._normalize_reactive_pair(raw_pair)
            if normalized not in selected_pairs:
                continue
            for state_name in state_entries.keys():
                bonded_reactive_by_state.setdefault(state_name, set()).add(normalized)

        # Resolve per-atom (C, a) via element lookup; fall back to scalar defaults.
        default_c = float(sp_c_kcal_mol_per_atom)
        default_a = float(sp_a_per_angstrom_per_atom)
        per_type = soft_pair_parameters_by_atom_type or {}
        init_atoms = list(self._modellers[initial_state].topology.atoms())

        def _per_atom_ca(scaffold_idx: int) -> tuple[float, float]:
            element = init_atoms[scaffold_idx].element
            symbol = element.symbol if element is not None else None
            if symbol is not None and symbol in per_type:
                c_val, a_val = per_type[symbol]
                return float(c_val), float(a_val)
            return default_c, default_a

        # Precompute per-pair (A_cross, B_cross) from geometric combining of
        # the per-atom (C, a). OpenMM units: A in kJ/mol, B in 1/nm.
        KJ_PER_KCAL = 4.184
        per_pair_ab: Dict[CanonicalPairKey, tuple[float, float]] = {}
        for pair in selected_pairs:
            if pair in explicit_ab:
                per_pair_ab[pair] = explicit_ab[pair]  # explicit per-pair (A, B) wins
                continue
            idx_a, idx_b = self._particle_indices_from_canonical_pair(initial_state, pair)
            c_a, a_a = _per_atom_ca(idx_a)
            c_b, a_b = _per_atom_ca(idx_b)
            A_cross_kj = c_a * c_b * KJ_PER_KCAL
            B_cross_per_nm = a_a * a_b / 0.1
            per_pair_ab[pair] = (A_cross_kj, B_cross_per_nm)

        # Extract state-specific charge products from each per-state System.
        state_systems = self._state_systems
        charge_products_by_state: Dict[str, Dict[CanonicalPairKey, float]] = {}
        for state_name in self.state_names:
            omm_sys = state_systems[state_name]
            nb = None
            for fi in range(omm_sys.getNumForces()):
                f = omm_sys.getForce(fi)
                if isinstance(f, openmm.NonbondedForce):
                    nb = f
                    break
            if nb is None:
                continue
            per_state: Dict[CanonicalPairKey, float] = {}
            for pair in selected_pairs:
                idx_a, idx_b = self._particle_indices_from_canonical_pair(state_name, pair)
                q_a, _, _ = nb.getParticleParameters(idx_a)
                q_b, _, _ = nb.getParticleParameters(idx_b)
                q_a_e = float(q_a.value_in_unit(unit.elementary_charge))
                q_b_e = float(q_b.value_in_unit(unit.elementary_charge))
                per_state[pair] = q_a_e * q_b_e
            charge_products_by_state[state_name] = per_state

        added: Dict[str, openmm.CustomBondForce] = {}
        for state_name in self.state_names:
            per_state_charges = charge_products_by_state.get(state_name, {})
            bonded_here = bonded_reactive_by_state.get(state_name, set())
            state_parameters: Dict[CanonicalPairKey, tuple[float, float, float]] = {}
            for pair in sorted(selected_pairs):
                if pair in bonded_here:
                    continue  # Q6 parity: Morse handles bond energy; no softpair Coul/V_a.
                A_cross_kj, B_cross_per_nm = per_pair_ab[pair]
                qiqj = per_state_charges.get(pair, 0.0)
                state_parameters[pair] = (A_cross_kj, B_cross_per_nm, qiqj)
            if not state_parameters:
                continue

            lambda_param = lambda_map[state_name]
            force = build_reactive_softpair_force(
                lambda_param=lambda_param,
                state_parameters=state_parameters,
                initial_state=initial_state,
                particle_indices_from_canonical_pair=self._particle_indices_from_canonical_pair,
            )
            system.addForce(force)
            self._register_force_to_group(force, term="nonbonded", state=state_name)
            self._state_specific_reactive_vdw_forces[state_name] = force
            added[state_name] = force

        if not added:
            raise ValueError("No reactive Q6 softpair forces were created.")

        # Exclude the reactive pairs from ALL other nonbonded / exception
        # forces so the reactive soft-pair is the SOLE nonbonded contribution for
        # that pair. The scaffold has three flavours to handle:
        #
        # 1. CustomNonbondedForce (shared + per-state) — use addExclusion.
        # 2. NonbondedForce (standard OpenMM) — use addException with zeros.
        # 3. Per-state LJ/Coulomb "exception" CustomBondForces (that encode
        #    the bond-based vdW/electrostatic contribution for pairs that
        #    WOULD otherwise be excluded from NB) — zero the per-bond
        #    epsilon / charge_prod entries for the reactive pair so they
        #    contribute zero energy.
        if exclude_from_shared_nonbonded:
            added_force_ids = {id(f) for f in added.values()}
            for pair in sorted(selected_pairs):
                idx_a, idx_b = self._particle_indices_from_canonical_pair(initial_state, pair)
                key = (min(idx_a, idx_b), max(idx_a, idx_b))
                for fi in range(system.getNumForces()):
                    f = system.getForce(fi)
                    if id(f) in added_force_ids:
                        continue  # Don't touch our own reactive soft-pair force.
                    if isinstance(f, openmm.CustomNonbondedForce):
                        existing = set()
                        for ei in range(f.getNumExclusions()):
                            ea, eb = f.getExclusionParticles(ei)
                            existing.add((min(ea, eb), max(ea, eb)))
                        if key not in existing:
                            f.addExclusion(idx_a, idx_b)
                    elif isinstance(f, openmm.NonbondedForce):
                        existing = set()
                        for ei in range(f.getNumExceptions()):
                            params = f.getExceptionParameters(ei)
                            ea, eb = params[0], params[1]
                            existing.add((min(ea, eb), max(ea, eb)))
                        if key not in existing:
                            f.addException(idx_a, idx_b, 0.0, 1.0, 0.0)
                    elif isinstance(f, openmm.CustomBondForce):
                        # State-specific LJ/Coulomb exception bond forces carry
                        # per-bond params like epsilon/sigma or charge_prod.
                        # Zero the epsilon and charge_prod for the reactive
                        # pair to suppress their contribution.
                        expr = f.getEnergyFunction()
                        has_lj  = "epsilon" in expr and "sigma" in expr
                        has_coul = "charge_prod" in expr
                        if not (has_lj or has_coul):
                            continue
                        n_per = f.getNumPerBondParameters()
                        param_names = [f.getPerBondParameterName(pi) for pi in range(n_per)]
                        for bi in range(f.getNumBonds()):
                            a, b, params = f.getBondParameters(bi)
                            if (min(a, b), max(a, b)) != key:
                                continue
                            new_params = list(params)
                            for pi, name in enumerate(param_names):
                                if name in ("epsilon", "charge_prod"):
                                    new_params[pi] = 0.0
                            f.setBondParameters(bi, a, b, new_params)

        return MappingProxyType(added)

    def add_state_dependent_onethree_exclusion_patches(
        self,
        *,
        lambda_parameter_names: Mapping[str, str] | None = None,
    ) -> Mapping[str, openmm.CustomBondForce]:
        """Patch state-dependent 1-3 exclusion mismatches across EVB diabats.

        The scaffold's NonbondedForce exclusion list is built once from the
        ``initial_state`` bond topology — so 1-3 pairs like (H, Cl_leaving)
        that are excluded in the reactant (connected via C-Cl_leaving) but
        should be fully interacting in the product (where C-Cl_leaving bond
        is absent) have their Coulomb/vdW silently suppressed. For identity
        SN2, this breaks the PMF symmetry by ~30 kcal/mol.

        For each atom pair whose 1-3 (or 1-2) status **differs** between the
        initial state and any other state, this method:

        1. Ensures the pair is excluded from the shared nonbonded force(s)
           (and from state-specific exception bonds) in both states.
        2. Adds a state-λ-weighted ``CustomBondForce`` restoring the correct
           effective LJ + Coulomb interaction per state, using each state's
           own per-atom parameters.

        Net effect at the active state's indicator λ=1: the pair contributes
        exactly what the state's topology would have, matching Q6's own
        per-state Q-atom energy evaluation.

        Call AFTER :meth:`add_vdw_forces`, :meth:`add_bond_forces`, and any
        reactive-pair soft-pair installer (Q6 or soft-exp) — but BEFORE any
        Context is created.
        """
        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("Call 'build_evb_system' before patching state-dependent exclusions.")
        if self._shared_nonbonded_force is None:
            raise ValueError("Call 'add_vdw_forces' before patching state-dependent exclusions.")

        lambda_map = self._resolve_lambda_parameters(lambda_parameter_names)
        initial_state = self._evb_system_source_state
        system = self._evb_system_scaffold

        # Build a helper that translates a state's own local atom index to
        # the scaffold's (initial_state's) atom index via canonical labels.
        def _state_to_scaffold_index(state_name: str, local_idx: int) -> int:
            label = self._canonical_atom_label(state_name, local_idx)
            # The initial_state's canonical label map is an identity, so we
            # look the label up to find the scaffold-space atom index.
            init_map = self._state_canonical_label_index_map(initial_state)
            return int(init_map[label])

        # 1. Enumerate per-state 1-2 and 1-3 pairs from each state's
        # HarmonicBondForce (the authoritative per-state bond topology),
        # then translate atom indices to scaffold space.
        def _harmonic_bonds(state_name: str) -> set[tuple[int, int]]:
            omm_sys = self._state_systems[state_name]
            pairs: set[tuple[int, int]] = set()
            for fi in range(omm_sys.getNumForces()):
                f = omm_sys.getForce(fi)
                if isinstance(f, openmm.HarmonicBondForce):
                    for bi in range(f.getNumBonds()):
                        i, j, _r0, _k = f.getBondParameters(bi)
                        si = _state_to_scaffold_index(state_name, i)
                        sj = _state_to_scaffold_index(state_name, j)
                        pairs.add((min(si, sj), max(si, sj)))
            return pairs

        def _one_three_from_bonds(bond_pairs: set[tuple[int, int]]) -> set[tuple[int, int]]:
            adj: Dict[int, set[int]] = {}
            for i, j in bond_pairs:
                adj.setdefault(i, set()).add(j)
                adj.setdefault(j, set()).add(i)
            onethree: set[tuple[int, int]] = set()
            for k, nbrs in adj.items():
                nbr_list = sorted(nbrs)
                for a_idx in range(len(nbr_list)):
                    for b_idx in range(a_idx + 1, len(nbr_list)):
                        a, b = nbr_list[a_idx], nbr_list[b_idx]
                        # Only count as 1-3 if a and b are not themselves 1-2 neighbors
                        if (min(a, b), max(a, b)) not in bond_pairs:
                            onethree.add((min(a, b), max(a, b)))
            return onethree

        def _one_four_from_bonds(
            bond_pairs: set[tuple[int, int]],
            one_three: set[tuple[int, int]],
        ) -> set[tuple[int, int]]:
            adj: Dict[int, set[int]] = {}
            for i, j in bond_pairs:
                adj.setdefault(i, set()).add(j)
                adj.setdefault(j, set()).add(i)
            onefour: set[tuple[int, int]] = set()
            for b, c in bond_pairs:  # central bond b-c of an a-b-c-d path
                for a in adj.get(b, ()):  # noqa: E741
                    if a == c:
                        continue
                    for d in adj.get(c, ()):
                        if d == b or d == a:
                            continue
                        key = (min(a, d), max(a, d))
                        if key in bond_pairs or key in one_three:
                            continue  # closer than 1-4 (ring); skip
                        onefour.add(key)
            return onefour

        state_12: Dict[str, set[tuple[int, int]]] = {}
        state_13: Dict[str, set[tuple[int, int]]] = {}
        state_14: Dict[str, set[tuple[int, int]]] = {}
        for s in self.state_names:
            bonds = _harmonic_bonds(s)
            state_12[s] = bonds
            state_13[s] = _one_three_from_bonds(bonds)
            state_14[s] = _one_four_from_bonds(bonds, state_13[s])

        # Per-state effective exclusion set (1-2 ∪ 1-3) in scaffold indices.
        state_excl = {s: state_12[s] | state_13[s] for s in self.state_names}
        init_excl = state_excl[initial_state]

        # Collect pairs already handled by a reactive softpair installer
        # (Q6 or soft-exp) so we don't double-count them here.
        reactive_already_handled: set[tuple[int, int]] = set()
        for rf in self._state_specific_reactive_vdw_forces.values():
            if isinstance(rf, openmm.CustomBondForce):
                for bi in range(rf.getNumBonds()):
                    a, b, _p = rf.getBondParameters(bi)
                    reactive_already_handled.add((min(a, b), max(a, b)))

        # 2. Find ALL pairs whose exclusion status differs between ANY two
        # states (excluding reactive pairs already handled).
        # and construct the per-state correction bond force.
        KJ_PER_KCAL = 4.184
        ONE_4PI_EPS0 = 138.935456  # kJ/mol·nm/e²
        added: Dict[str, openmm.CustomBondForce] = {}

        # Standard state-specific Coulomb+LJ form, no softcore — 1-3 pair
        # distances are moderate so raw LJ is fine here. `sign` is +1 when
        # the state expects the pair to be FREE but initial has it excluded
        # (we ADD the interaction in that state). `sign` is -1 when the state
        # expects the pair EXCLUDED but initial has it free (we SUBTRACT).
        def _scaffold_to_state_index(state_name: str, scaffold_idx: int) -> int:
            # Find the canonical label at this scaffold index via the initial
            # state's map, then look it up in the target state's map.
            init_map = self._state_canonical_label_index_map(initial_state)
            # Reverse initial map: scaffold_idx → canonical label
            label = None
            for lab, sidx in init_map.items():
                if int(sidx) == int(scaffold_idx):
                    label = lab
                    break
            if label is None:
                raise RuntimeError(f"Scaffold index {scaffold_idx} has no canonical label.")
            target_map = self._state_canonical_label_index_map(state_name)
            return int(target_map[label])

        def _get_state_lj_q(state_name: str, pair: tuple[int, int]) -> tuple[float, float, float]:
            omm_sys = self._state_systems[state_name]
            nb = None
            for fi in range(omm_sys.getNumForces()):
                f = omm_sys.getForce(fi)
                if isinstance(f, openmm.NonbondedForce):
                    nb = f
                    break
            if nb is None:
                return (0.0, 0.0, 0.0)
            # Translate scaffold-space indices to this state's local indices.
            li = _scaffold_to_state_index(state_name, pair[0])
            lj = _scaffold_to_state_index(state_name, pair[1])
            q_i, sig_i, eps_i = nb.getParticleParameters(li)
            q_j, sig_j, eps_j = nb.getParticleParameters(lj)
            # Geometric combining rule (matches Q6's default + OpenFF/OPLS-AA).
            sig_ij = math.sqrt(
                float(sig_i.value_in_unit(unit.nanometer))
                * float(sig_j.value_in_unit(unit.nanometer))
            )
            eps_ij = math.sqrt(
                float(eps_i.value_in_unit(unit.kilojoules_per_mole))
                * float(eps_j.value_in_unit(unit.kilojoules_per_mole))
            )
            qiqj = (
                float(q_i.value_in_unit(unit.elementary_charge))
                * float(q_j.value_in_unit(unit.elementary_charge))
            )
            return sig_ij, eps_ij, qiqj

        def _get_state_onefour(state_name: str, pair: tuple[int, int]):
            """FF-scaled 1-4 ``(sigma_nm, eps_kj, qiqj)`` from the state's own
            NonbondedForce exception (the authoritative per-state 1-4 scaling),
            or ``None`` if the pair has no 1-4 exception in that state."""
            omm_sys = self._state_systems[state_name]
            nb = None
            for fi in range(omm_sys.getNumForces()):
                f = omm_sys.getForce(fi)
                if isinstance(f, openmm.NonbondedForce):
                    nb = f
                    break
            if nb is None:
                return None
            li = _scaffold_to_state_index(state_name, pair[0])
            lj = _scaffold_to_state_index(state_name, pair[1])
            key = (min(li, lj), max(li, lj))
            for ei in range(nb.getNumExceptions()):
                a, b, qp, sig, eps = nb.getExceptionParameters(ei)
                if (min(a, b), max(a, b)) == key:
                    return (
                        float(sig.value_in_unit(unit.nanometer)),
                        float(eps.value_in_unit(unit.kilojoules_per_mole)),
                        float(qp.value_in_unit(unit.elementary_charge**2)),
                    )
            return None

        shared_nb_exclusions = set()
        shared_nb = self._shared_nonbonded_force
        for ei in range(shared_nb.getNumExclusions()):
            a, b = shared_nb.getExclusionParticles(ei)
            shared_nb_exclusions.add((min(a, b), max(a, b)))

        # A pair is "mismatched" if its exclusion set membership differs
        # across any pair of states. For each state, that state's patch
        # lists pairs the state expects to be FREE but which exist in the
        # union of all states' exclusions (so at least one other state
        # excludes it → the scaffold's topology might reflect that exclusion).
        all_excl_union = set().union(*[state_excl[s] for s in self.state_names])
        all_excl_inter = set.intersection(*[state_excl[s] for s in self.state_names])
        mismatched_pairs = {p for p in (all_excl_union - all_excl_inter)
                            if p not in reactive_already_handled}
        all_mismatched_pairs = mismatched_pairs

        # State-dependent 1-4 pairs: a pair that is a (scaled) 1-4 in one diabat
        # but a different class in another -- a full >=1-5 pair, a 1-2/1-3
        # exclusion, or simply not 1-4. Restoring these as FULL interactions (the
        # legacy behaviour) mis-scales them, because the active state expects the
        # FF's 1-4 scale factor (OpenFF 0.5 LJ / 0.833 Coulomb). Q6 handles this
        # natively: it rebuilds ``qconn`` per state and carries a per-state 1-4
        # list, so each diabat applies its own 1-4 scaling. We reproduce that by
        # folding these pairs into the per-state correction below with each
        # state's OWN scaled 1-4 parameters (read from that state's NonbondedForce
        # exception) where the pair is 1-4, and the full interaction where it is
        # >=1-5.
        onefour_union = set().union(*[state_14[s] for s in self.state_names])
        onefour_inter = set.intersection(*[state_14[s] for s in self.state_names])
        onefour_mismatch = {
            p for p in (onefour_union - onefour_inter) if p not in reactive_already_handled
        }
        onefour_in_mismatched = {p for p in all_mismatched_pairs if p in onefour_union}
        problematic_14 = onefour_mismatch | onefour_in_mismatched
        self._state_dependent_onefour_pairs = frozenset(problematic_14)
        self._onethree_patch_applied = True

        # Every pair whose per-state nonbonded class is non-uniform across diabats.
        # These are excluded from all scaffold nonbonded forces (step 3) and rebuilt
        # by per-state corrections (step 4) so each state sees exactly its own
        # topology's interaction (excluded / scaled-1-4 / full).
        pairs_to_patch = set(all_mismatched_pairs) | problematic_14

        # For each state, the patched pairs it does NOT exclude as 1-2/1-3: each
        # gets a per-state correction bond -- scaled 1-4 if the pair is 1-4 in that
        # state (handled in step 4), otherwise the full >=1-5 interaction.
        per_state_free_patches: Dict[str, list[tuple[int, int]]] = {}
        for state_name in self.state_names:
            per_state_free_patches[state_name] = sorted(pairs_to_patch - state_excl[state_name])

        if not pairs_to_patch:
            return MappingProxyType({})  # nothing to patch

        # 3. Exclude every patched pair from shared/state nonbonded forces (and
        # zero them in existing exception bond forces) so only our per-state
        # correction bond forces contribute for these pairs. Existing exceptions
        # are force-overwritten to a full exclusion (chargeProd=0, epsilon=0) --
        # critical for a pair that is a scaled 1-4 in the initial state, whose
        # exception would otherwise survive and double-count against the
        # correction.
        new_softpair_force_ids: set[int] = set()
        for fi in range(system.getNumForces()):
            f = system.getForce(fi)
            if isinstance(f, openmm.CustomNonbondedForce):
                existing = set()
                for ei in range(f.getNumExclusions()):
                    a, b = f.getExclusionParticles(ei)
                    existing.add((min(a, b), max(a, b)))
                for pair in sorted(pairs_to_patch):
                    if pair not in existing:
                        f.addExclusion(pair[0], pair[1])
            elif isinstance(f, openmm.NonbondedForce):
                existing_idx: Dict[tuple[int, int], int] = {}
                for ei in range(f.getNumExceptions()):
                    params = f.getExceptionParameters(ei)
                    a, b = params[0], params[1]
                    existing_idx[(min(a, b), max(a, b))] = ei
                for pair in sorted(pairs_to_patch):
                    if pair in existing_idx:
                        f.setExceptionParameters(existing_idx[pair], pair[0], pair[1], 0.0, 1.0, 0.0)
                    else:
                        f.addException(pair[0], pair[1], 0.0, 1.0, 0.0)
            elif isinstance(f, openmm.CustomBondForce):
                expr = f.getEnergyFunction()
                has_lj = ("epsilon" in expr) and ("sigma" in expr)
                has_coul = "charge_prod" in expr
                if not (has_lj or has_coul):
                    continue
                n_per = f.getNumPerBondParameters()
                param_names = [f.getPerBondParameterName(pi) for pi in range(n_per)]
                for bi in range(f.getNumBonds()):
                    a, b, params = f.getBondParameters(bi)
                    if (min(a, b), max(a, b)) not in pairs_to_patch:
                        continue
                    new_params = list(params)
                    for pi, name in enumerate(param_names):
                        if name in ("epsilon", "charge_prod"):
                            new_params[pi] = 0.0
                    f.setBondParameters(bi, a, b, new_params)

        # 4. For each state, add a correction CustomBondForce that
        # re-introduces the LJ+Coulomb interaction for pairs the state
        # considers FREE (using the state's own parameters). Each patch is
        # λ_state-scaled so it contributes fully only when that state is
        # the active indicator, and 0 when any other state is active. This
        # correctly implements "each state's topology defines its own
        # exclusions" for mapping/indicator-λ workflows.
        for state_name, free_pairs in per_state_free_patches.items():
            if not free_pairs:
                continue
            lambda_param = lambda_map[state_name]
            expr = (
                f"{lambda_param} * ("
                f" 4*epsilon*((sigma/r)^12 - (sigma/r)^6) "
                f" + k_coul * qiqj / r )"
            )
            force = openmm.CustomBondForce(expr)
            force.addGlobalParameter(lambda_param, 0.0)
            force.addGlobalParameter("k_coul", ONE_4PI_EPS0)
            force.addPerBondParameter("sigma")
            force.addPerBondParameter("epsilon")
            force.addPerBondParameter("qiqj")
            for pair in free_pairs:
                if pair in state_14[state_name]:
                    # 1-4 in this state -> use the FF's own scaled 1-4 parameters
                    # (Q6-equivalent per-state 1-4 handling), not the full LJ+Coul.
                    scaled = _get_state_onefour(state_name, pair)
                    sigma_nm, eps_kj, qiqj = scaled if scaled is not None else _get_state_lj_q(state_name, pair)
                else:
                    sigma_nm, eps_kj, qiqj = _get_state_lj_q(state_name, pair)
                force.addBond(pair[0], pair[1], (sigma_nm, eps_kj, qiqj))
            system.addForce(force)
            self._register_force_to_group(force, term="nonbonded", state=state_name)
            added[state_name] = force

        return MappingProxyType(added)

    def add_bond_forces(
        self,
        *,
        lambda_parameter_names: Mapping[str, str] | None = None,
    ) -> Mapping[str, openmm.Force]:
        """
        Attach shared and state-specific bond forces to the EVB system scaffold.

        Parameters
        ----------
        lambda_parameter_names :
            Optional mapping of ``state_name -> global parameter name`` controlling the
            CustomBondForce scaling for that state's unique/conflicting bonds. Defaults
            to ``lambda_<state_name>`` when omitted.

        Returns
        -------
        Mapping[str, openmm.Force]
            Mapping containing ``\"shared\"`` for the universal bond force and one entry
            per state for the lambda-scaled bond forces that were added.
        """

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("Call 'build_evb_system' before adding bond forces.")
        if self._shared_bond_force is not None or self._state_specific_bond_forces:
            raise RuntimeError("Bond forces have already been added to this EVB system.")

        lambda_map = self._resolve_lambda_parameters(lambda_parameter_names)
        shared_entries, per_state_entries = self._shared_and_state_terms("bonds")
        system = self._evb_system_scaffold
        initial_state = self._evb_system_source_state

        shared_force = build_shared_bond_force(
            shared_entries=shared_entries,
            initial_state=initial_state,
            particle_indices_from_canonical_key=self._particle_indices_from_canonical_key,
        )
        system.addForce(shared_force)
        self._register_force_to_group(shared_force, term="shared", state=None)
        self._shared_bond_force = shared_force

        state_forces: Dict[str, openmm.CustomBondForce] = {}
        for state_name, entries in per_state_entries.items():
            lambda_param = lambda_map[state_name]
            force = build_state_bond_force(
                lambda_param=lambda_param,
                entries=entries,
                initial_state=initial_state,
                particle_indices_from_canonical_key=self._particle_indices_from_canonical_key,
            )
            system.addForce(force)
            self._register_force_to_group(force, term="bond", state=state_name)
            state_forces[state_name] = force
        self._state_specific_bond_forces = state_forces

        added_forces: Dict[str, openmm.Force] = {"shared": shared_force}
        added_forces.update(state_forces)
        return MappingProxyType(added_forces)

    def add_angle_forces(
        self,
        *,
        lambda_parameter_names: Mapping[str, str] | None = None,
    ) -> Mapping[str, openmm.Force]:
        """
        Attach shared and state-specific angle forces to the EVB system scaffold.
        """

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("Call 'build_evb_system' before adding angle forces.")
        if self._shared_angle_force is not None or self._state_specific_angle_forces:
            raise RuntimeError("Angle forces have already been added to this EVB system.")

        lambda_map = self._resolve_lambda_parameters(lambda_parameter_names)
        shared_entries, per_state_entries = self._shared_and_state_terms("angles")
        system = self._evb_system_scaffold
        initial_state = self._evb_system_source_state

        shared_force = build_shared_angle_force(
            shared_entries=shared_entries,
            initial_state=initial_state,
            particle_indices_from_canonical_key=self._particle_indices_from_canonical_key,
        )
        system.addForce(shared_force)
        self._register_force_to_group(shared_force, term="shared", state=None)
        self._shared_angle_force = shared_force

        state_forces: Dict[str, openmm.CustomAngleForce] = {}
        for state_name, entries in per_state_entries.items():
            lambda_param = lambda_map[state_name]
            force = build_state_angle_force(
                lambda_param=lambda_param,
                entries=entries,
                initial_state=initial_state,
                particle_indices_from_canonical_key=self._particle_indices_from_canonical_key,
            )
            system.addForce(force)
            self._register_force_to_group(force, term="angle", state=state_name)
            state_forces[state_name] = force
        self._state_specific_angle_forces = state_forces

        added_forces: Dict[str, openmm.Force] = {"shared": shared_force}
        added_forces.update(state_forces)
        return MappingProxyType(added_forces)

    def add_torsion_forces(
        self,
        *,
        lambda_parameter_names: Mapping[str, str] | None = None,
    ) -> Mapping[str, openmm.Force]:
        """
        Attach shared and state-specific torsion forces to the EVB system scaffold.
        """

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("Call 'build_evb_system' before adding torsion forces.")
        if self._shared_torsion_force is not None or self._state_specific_torsion_forces:
            raise RuntimeError("Torsion forces have already been added to this EVB system.")

        lambda_map = self._resolve_lambda_parameters(lambda_parameter_names)
        shared_entries, per_state_entries = self._shared_and_state_terms("torsions")
        system = self._evb_system_scaffold
        initial_state = self._evb_system_source_state

        shared_force = build_shared_torsion_force(
            shared_entries=shared_entries,
            initial_state=initial_state,
            particle_indices_from_canonical_key=self._particle_indices_from_canonical_key,
            torsion_components=self._torsion_components,
        )
        system.addForce(shared_force)
        self._register_force_to_group(shared_force, term="shared", state=None)
        self._shared_torsion_force = shared_force

        state_forces: Dict[str, openmm.CustomTorsionForce] = {}
        for state_name, entries in per_state_entries.items():
            lambda_param = lambda_map[state_name]
            force = build_state_torsion_force(
                lambda_param=lambda_param,
                entries=entries,
                initial_state=initial_state,
                particle_indices_from_canonical_key=self._particle_indices_from_canonical_key,
                torsion_components=self._torsion_components,
            )
            system.addForce(force)
            self._register_force_to_group(force, term="torsion", state=state_name)
            state_forces[state_name] = force
        self._state_specific_torsion_forces = state_forces

        added_forces: Dict[str, openmm.Force] = {"shared": shared_force}
        added_forces.update(state_forces)
        return MappingProxyType(added_forces)

    def get_evb_forces(self) -> Mapping[str, openmm.Force]:
        """
        Return a mapping of all EVB-specific forces currently attached to the scaffold.
        """

        if self._evb_system_scaffold is None:
            raise ValueError("EVB base system not constructed yet; call 'build_evb_system'.")

        forces: Dict[str, openmm.Force] = {}

        if self._shared_nonbonded_force is not None:
            forces["nonbonded_shared"] = self._shared_nonbonded_force
            forces["nonbonded_shared_mm"] = self._shared_nonbonded_force
        if self._shared_nonbonded_qm_force is not None:
            forces["nonbonded_shared_qm"] = self._shared_nonbonded_qm_force
        if self._shared_nonbonded_exception_force is not None:
            forces["nonbonded_shared_exceptions_mm"] = self._shared_nonbonded_exception_force
        if self._shared_nonbonded_qm_exception_force is not None:
            forces["nonbonded_shared_exceptions_qm"] = self._shared_nonbonded_qm_exception_force
        for state_name, force in self._state_specific_nonbonded_forces.items():
            forces[f"nonbonded_{state_name}"] = force
            forces[f"nonbonded_{state_name}_mm"] = force
        for state_name, force in self._state_specific_nonbonded_qm_forces.items():
            forces[f"nonbonded_{state_name}_qm"] = force
        for state_name, force in self._state_specific_vdw_exception_forces.items():
            forces[f"vdw_{state_name}_exceptions"] = force
            forces[f"vdw_{state_name}_exceptions_mm"] = force
        for state_name, force in self._state_specific_vdw_qm_exception_forces.items():
            forces[f"vdw_{state_name}_exceptions_qm"] = force
        for state_name, force in self._state_specific_reactive_vdw_forces.items():
            forces[f"vdw_{state_name}_softpair"] = force
        for state_name, force in self._state_specific_reactive_morse_bond_forces.items():
            forces[f"bond_{state_name}_reactive_morse"] = force
        for state_name, force in self._state_specific_electrostatics_exception_forces.items():
            forces[f"elec_{state_name}_exceptions"] = force
            forces[f"elec_{state_name}_exceptions_mm"] = force
        for state_name, force in self._state_specific_electrostatics_qm_exception_forces.items():
            forces[f"elec_{state_name}_exceptions_qm"] = force

        if self._shared_bond_force is not None:
            forces["bond_shared"] = self._shared_bond_force
        for state_name, force in self._state_specific_bond_forces.items():
            forces[f"bond_{state_name}"] = force

        if self._shared_angle_force is not None:
            forces["angle_shared"] = self._shared_angle_force
        for state_name, force in self._state_specific_angle_forces.items():
            forces[f"angle_{state_name}"] = force

        if self._shared_torsion_force is not None:
            forces["torsion_shared"] = self._shared_torsion_force
        for state_name, force in self._state_specific_torsion_forces.items():
            forces[f"torsion_{state_name}"] = force

        if self._scaas_runtime is not None:
            forces["scaas_solvent_radial"] = self._scaas_runtime.radial_force
            if self._scaas_runtime.shell_fix_force is not None:
                forces["scaas_shell_fix"] = self._scaas_runtime.shell_fix_force
            if self._scaas_runtime.polarisation_force is not None:
                forces["scaas_polarisation"] = self._scaas_runtime.polarisation_force

        return MappingProxyType(forces)

    def enable_scaas(
        self,
        *,
        solvent_radius_angstrom: float,
        center: unit.Quantity | Sequence[float] | openmm.Vec3 | None = None,
        temperature: unit.Quantity | float = 300.0,
        shell_radius_fraction: float = 0.85,
        shell_force_kcal_mol_a2: float = 10.0,
        excluded_force_kcal_mol_a2: float = 200.0,
        radial_force_kcal_mol_a2: float = 60.0,
        morse_depth_kcal_mol: float | None = None,
        morse_width_per_angstrom: float | None = None,
        radial_wall_type: str = "q6_morse",
        polarisation: bool = True,
        charge_correction: bool | None = None,
        polarisation_force_kcal_mol_rad2: float = 20.0,
        pol_update_steps: int = 100,
        solvent_residue_names: Sequence[str] | None = None,
        solvent_oxygen_name: str | None = None,
        solvent_hydrogen_names: Sequence[str] | None = None,
        solvent_density_per_a3: float = 0.0334,
        water_dipole_ea: float = 0.49,
        dielectric_constant: float = 78.4,
        solute_total_charge_e: float | None = None,
        excluded_atom_indices: Sequence[int] | None = None,
        shell_atom_indices: Sequence[int] | None = None,
        excluded_freeze: bool = False,
        exclude_bonded: bool = False,
        solute_com_restraint_kcal_mol_a2: float = 0.0,
        solute_com_restraint_radius_angstrom: float = 0.0,
        solute_com_restraint_atom_indices: Sequence[int] | None = None,
        backend: str = "cpp",
    ) -> Mapping[str, object]:
        """
        Enable Q6-style SCAAS restraints on the EVB scaffold.

        Solute centering
        ----------------
        The SCAAS boundary (radial wall, shell-fix) restrains only *solvent*;
        the solute is otherwise free to drift toward the droplet boundary, where
        it becomes under-solvated and the per-window solvation (hence the EVB
        gap) varies. Set ``solute_com_restraint_kcal_mol_a2 > 0`` to add a
        flat-bottom harmonic restraint of the solute centre of mass to the sphere
        centre: ``0.5 k (r - r0)^2`` for ``r > r0``, zero inside ``r0``. With
        ``solute_com_restraint_radius_angstrom = 0`` it is a pure harmonic to the
        centre. ``solute_com_restraint_atom_indices`` selects the restrained
        atoms (default: all non-solvent atoms, i.e. the solute complex).

        Notes
        -----
        Backend selection:
        - ``backend="cpp"`` is the only supported public SCAAS backend and the current default.
        - Passing any other backend value raises ``ValueError``.
        - ``backend="cpp"`` requires optional native backend support and raises
          if that support is unavailable.
        - On CUDA, the conservative host-staged native path remains the default.
          The faster hybrid device-update path can be enabled experimentally with
          ``EVBOPENMM_SCAAS_EXPERIMENTAL_DEVICE_UPDATE=1`` for exploratory or
          benchmarking runs. Keep it disabled for conservative/reference runs.
        """

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("EVB system not initialized; call 'build_evb_system' and add forces first.")
        if self._boundary_mode == "periodic":
            raise RuntimeError(
                "SCAAS spherical restraints are incompatible with periodic boundary mode. "
                "Use 'solvate_scaas_sphere_states' for droplet runs, not "
                "'solvate_periodic_box_states'."
            )
        if self._scaas_runtime is not None:
            raise RuntimeError(
                "SCAAS is already enabled on this EVBSystem instance. Rebuild the EVB scaffold to change settings."
            )

        backend_resolution = resolve_scaas_backend(str(backend))
        backend_key = backend_resolution.selected_backend
        if excluded_freeze:
            raise NotImplementedError(
                "'excluded_freeze=True' requires per-step coordinate pinning and is not implemented yet."
            )
        if exclude_bonded:
            raise NotImplementedError(
                "'exclude_bonded=True' is not implemented in this OpenMM integration yet."
            )

        radius_a = float(solvent_radius_angstrom)
        if radius_a <= 0.0:
            raise ValueError("'solvent_radius_angstrom' must be > 0.")
        _units.warn_if_scaas_radius_looks_like_nm(radius_a, label="solvent_radius_angstrom")

        temperature_q = _units.coerce_temperature(temperature)
        temperature_kelvin = float(temperature_q.value_in_unit(unit.kelvin))
        residue_names = normalize_residue_names(solvent_residue_names)
        water_molecules = self._create_scaas_water_list(
            state_name=self._evb_system_source_state,
            residue_names=residue_names,
            oxygen_name=solvent_oxygen_name,
            hydrogen_names=solvent_hydrogen_names,
        )
        if not water_molecules:
            raise ValueError(
                "No compatible 3-site solvent molecules were found for SCAAS. "
                f"Residue filter: {residue_names!r}."
            )

        if center is None:
            center_nm = self._estimate_scaas_center_nm(water_molecules)
        else:
            _units.warn_if_center_unitless(center, context="enable_scaas")
            center_nm = normalize_scaas_center(center)

        charge_correction_enabled = bool(polarisation) if charge_correction is None else bool(charge_correction)
        shell_radius_nm = float(shell_radius_fraction) * radius_a * ANGSTROM_TO_NM
        if shell_radius_nm < 0.0:
            raise ValueError("'shell_radius_fraction' must be >= 0.")

        radial_params = build_q6_radial_boundary_parameters(
            solvent_radius_angstrom=radius_a,
            temperature_kelvin=temperature_kelvin,
            radial_force_kcal_mol_a2=float(radial_force_kcal_mol_a2),
            morse_depth_kcal_mol=(
                None if morse_depth_kcal_mol is None else float(morse_depth_kcal_mol)
            ),
            morse_width_per_angstrom=(
                None if morse_width_per_angstrom is None else float(morse_width_per_angstrom)
            ),
        )

        morse_depth_kcal = float(radial_params.morse_depth_kcal_mol)
        morse_width_per_a = float(radial_params.morse_width_per_angstrom)
        radial_force_kj_mol_nm2 = float(radial_params.radial_force_kcal_mol_a2) * KJ_PER_KCAL / (ANGSTROM_TO_NM**2)
        shell_force_kj_mol_nm2 = float(shell_force_kcal_mol_a2) * KJ_PER_KCAL / (ANGSTROM_TO_NM**2)
        excluded_force_kj_mol_nm2 = float(excluded_force_kcal_mol_a2) * KJ_PER_KCAL / (ANGSTROM_TO_NM**2)
        morse_depth_kj_mol = morse_depth_kcal * KJ_PER_KCAL
        morse_width_per_nm = morse_width_per_a / ANGSTROM_TO_NM
        polarisation_force_kj_mol_rad2 = float(polarisation_force_kcal_mol_rad2) * KJ_PER_KCAL

        if radial_force_kj_mol_nm2 <= 0.0:
            raise ValueError("'radial_force_kcal_mol_a2' must be > 0.")
        if polarisation and polarisation_force_kj_mol_rad2 <= 0.0:
            raise ValueError("'polarisation_force_kcal_mol_rad2' must be > 0 when polarisation is enabled.")
        if pol_update_steps <= 0:
            raise ValueError("'pol_update_steps' must be > 0.")

        solvent_radius_nm = radius_a * ANGSTROM_TO_NM
        radial_boundary_nm = float(radial_params.effective_boundary_angstrom) * ANGSTROM_TO_NM

        radial_delta = "sqrt((x-xc)^2+(y-yc)^2+(z-zc)^2)-r0"
        
        if radial_wall_type == "flat_harmonic":
            # Correct units: kr is [energy/distance^2], delta is [distance].
            # energy is kr * delta^2 which is [energy].
            radial_energy = f"scaas_enable*step({radial_delta})*(0.5*kr*({radial_delta})^2)"
        elif radial_wall_type == "flat_exponential":
            # For exponential, we use kr as an energy scale (e.g. 1 kT) if we want a soft wall.
            # But the user provided 60 kcal/mol/A^2. If we treat it as an energy prefactor, it explodes.
            # Fix: use a fixed small prefactor (1.0 kJ/mol) and alpha to define the shape.
            # Energy = 1.0 * exp(alpha * delta)
            radial_delta_pos = f"0.5*(({radial_delta})+abs({radial_delta}))"
            radial_energy = f"scaas_enable*(exp(alpha*({radial_delta_pos})) - 1)"
        elif radial_wall_type == "q6_morse":
            # Keep Morse exponent numerically safe: OpenMM evaluates all branches,
            # so clamp to non-positive values to avoid overflow far outside r0.
            radial_delta_neg = f"0.5*(({radial_delta})-abs({radial_delta}))"
            radial_exp = f"exp(alpha*({radial_delta_neg}))"
            radial_energy = (
                f"scaas_enable*(step({radial_delta})*(0.5*kr*({radial_delta})^2 - de) "
                f"+ (1-step({radial_delta}))*de*(({radial_exp})*({radial_exp}) - 2*({radial_exp})))"
            )
        else:
            raise ValueError(f"Unsupported radial_wall_type '{radial_wall_type}'")
        radial_force = openmm.CustomExternalForce(radial_energy)
        radial_force.addGlobalParameter("scaas_enable", 1.0)
        radial_force.addGlobalParameter("xc", center_nm[0])
        radial_force.addGlobalParameter("yc", center_nm[1])
        radial_force.addGlobalParameter("zc", center_nm[2])
        radial_force.addGlobalParameter("r0", radial_boundary_nm)
        radial_force.addGlobalParameter("kr", radial_force_kj_mol_nm2)
        radial_force.addGlobalParameter("de", morse_depth_kj_mol)
        radial_force.addGlobalParameter("alpha", morse_width_per_nm)
        for molecule in water_molecules:
            radial_force.addParticle(molecule.oxygen_index, [])
        self._evb_system_scaffold.addForce(radial_force)
        self._register_force_to_group(radial_force, term="scaas_radial", state=None)

        source_state = self._evb_system_source_state
        topology = self._modellers[source_state].topology
        initial_positions = self.evb_initial_positions
        solvent_atom_indices = {
            idx
            for molecule in water_molecules
            for idx in (molecule.oxygen_index, molecule.hydrogen1_index, molecule.hydrogen2_index)
        }

        if excluded_atom_indices is None or shell_atom_indices is None:
            auto_excluded: list[int] = []
            auto_shell: list[int] = []
            cx, cy, cz = center_nm
            for atom in topology.atoms():
                atom_index = atom.index
                if atom_index in solvent_atom_indices:
                    continue
                position = initial_positions[atom_index]
                dx = float(position.x) - cx
                dy = float(position.y) - cy
                dz = float(position.z) - cz
                radius = math.sqrt(dx * dx + dy * dy + dz * dz)
                if radius > solvent_radius_nm:
                    auto_excluded.append(atom_index)
                    continue
                atomic_number = atom.element.atomic_number if atom.element is not None else 0
                if atomic_number != 1 and radius > shell_radius_nm:
                    auto_shell.append(atom_index)
            if excluded_atom_indices is None:
                excluded_atom_indices = tuple(sorted(auto_excluded))
            if shell_atom_indices is None:
                shell_atom_indices = tuple(sorted(auto_shell))

        normalized_excluded = self._normalize_atom_index_sequence(
            excluded_atom_indices, expected_particles=self._evb_system_scaffold.getNumParticles()
        )
        normalized_shell = self._normalize_atom_index_sequence(
            shell_atom_indices, expected_particles=self._evb_system_scaffold.getNumParticles()
        )

        shell_fix_force: openmm.CustomExternalForce | None = None
        if normalized_excluded or normalized_shell:
            shell_fix_force = openmm.CustomExternalForce(
                "scaas_enable*0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
            )
            shell_fix_force.addGlobalParameter("scaas_enable", 1.0)
            shell_fix_force.addPerParticleParameter("x0")
            shell_fix_force.addPerParticleParameter("y0")
            shell_fix_force.addPerParticleParameter("z0")
            shell_fix_force.addPerParticleParameter("k")

            for atom_index in normalized_excluded:
                reference = initial_positions[atom_index]
                shell_fix_force.addParticle(
                    atom_index,
                    (
                        float(reference.x),
                        float(reference.y),
                        float(reference.z),
                        excluded_force_kj_mol_nm2,
                    ),
                )
            for atom_index in normalized_shell:
                if atom_index in normalized_excluded:
                    continue
                reference = initial_positions[atom_index]
                shell_fix_force.addParticle(
                    atom_index,
                    (
                        float(reference.x),
                        float(reference.y),
                        float(reference.z),
                        shell_force_kj_mol_nm2,
                    ),
                )

            self._evb_system_scaffold.addForce(shell_fix_force)
            self._register_force_to_group(shell_fix_force, term="scaas_shell_fix", state=None)

        polarisation_force: openmm.CustomCompoundBondForce | None = None
        shells = tuple()
        if polarisation:
            if charge_correction_enabled:
                total_charge_e = (
                    float(solute_total_charge_e)
                    if solute_total_charge_e is not None
                    else self._estimate_scaas_solute_charge(
                        state_name=source_state,
                        solvent_atom_indices=solvent_atom_indices,
                    )
                )
            else:
                total_charge_e = 0.0

            shells = build_q6_shells(
                solvent_radius_angstrom=radius_a,
                total_charge_e=total_charge_e,
                water_density_per_a3=float(solvent_density_per_a3),
                water_dipole_ea=float(water_dipole_ea),
                dielectric_constant=float(dielectric_constant),
            )

            pol_expression = (
                "0.5*scaas_enable*active*kpol*delta*delta;"
                "delta=theta-theta_target+theta_corr;"
                "theta=acos(cos_theta);"
                # Avoid acos' singular derivative at |cos_theta|=1.
                "cos_theta=max(-0.999999,min(0.999999,cos_raw));"
                "cos_raw=(mux*rcx+muy*rcy+muz*rcz)/den;"
                # Keep denominator away from zero when dipole/radial vectors become tiny.
                "den=sqrt(mux*mux+muy*muy+muz*muz)*sqrt(rcx*rcx+rcy*rcy+rcz*rcz)+1e-8;"
                "mux=x2+x3-2*x1;"
                "muy=y2+y3-2*y1;"
                "muz=z2+z3-2*z1;"
                "rcx=x1-xc;"
                "rcy=y1-yc;"
                "rcz=z1-zc"
            )
            polarisation_force = openmm.CustomCompoundBondForce(3, pol_expression)
            polarisation_force.addGlobalParameter("scaas_enable", 1.0)
            polarisation_force.addGlobalParameter("xc", center_nm[0])
            polarisation_force.addGlobalParameter("yc", center_nm[1])
            polarisation_force.addGlobalParameter("zc", center_nm[2])
            polarisation_force.addGlobalParameter("kpol", polarisation_force_kj_mol_rad2)
            polarisation_force.addGlobalParameter(
                "scaas_native_runtime",
                1.0 if backend_key == "cpp" else 0.0,
            )
            polarisation_force.addGlobalParameter("scaas_update_interval", float(pol_update_steps))
            polarisation_force.addGlobalParameter(
                "scaas_charge_correction",
                1.0 if charge_correction_enabled else 0.0,
            )
            polarisation_force.addGlobalParameter("scaas_shell_count", float(len(shells)))
            for shell_idx, shell in enumerate(shells):
                polarisation_force.addGlobalParameter(
                    f"scaas_shell_outer_{shell_idx}",
                    float(shell.outer_radius_nm),
                )
                polarisation_force.addGlobalParameter(
                    f"scaas_shell_thickness_{shell_idx}",
                    float(shell.thickness_nm),
                )
                polarisation_force.addGlobalParameter(
                    f"scaas_shell_cstb_{shell_idx}",
                    float(shell.cstb),
                )
            polarisation_force.addPerBondParameter("active")
            polarisation_force.addPerBondParameter("theta_target")
            polarisation_force.addPerBondParameter("theta_corr")

            for molecule in water_molecules:
                polarisation_force.addBond(
                    (molecule.oxygen_index, molecule.hydrogen1_index, molecule.hydrogen2_index),
                    (0.0, 0.0, 0.0),
                )

            self._evb_system_scaffold.addForce(polarisation_force)
            self._register_force_to_group(polarisation_force, term="scaas_polarisation", state=None)

        # ---- optional solute COM restraint (keep the reacting solute near the sphere centre) ----
        com_k_kcal = float(solute_com_restraint_kcal_mol_a2)
        if com_k_kcal < 0.0:
            raise ValueError("'solute_com_restraint_kcal_mol_a2' must be >= 0.")
        if com_k_kcal > 0.0:
            com_r0_a = float(solute_com_restraint_radius_angstrom)
            if com_r0_a < 0.0:
                raise ValueError("'solute_com_restraint_radius_angstrom' must be >= 0.")
            if solute_com_restraint_atom_indices is None:
                solvent_atom_set: Set[int] = set()
                for mol in water_molecules:
                    solvent_atom_set.update(
                        (mol.oxygen_index, mol.hydrogen1_index, mol.hydrogen2_index))
                com_atoms = [i for i in range(self._evb_system_scaffold.getNumParticles())
                             if i not in solvent_atom_set]
            else:
                com_atoms = [int(i) for i in solute_com_restraint_atom_indices]
            if not com_atoms:
                raise ValueError("solute COM restraint requested but no solute atoms were resolved.")
            com_k_kj = com_k_kcal * KJ_PER_KCAL / (ANGSTROM_TO_NM ** 2)
            com_r0_nm = com_r0_a * ANGSTROM_TO_NM
            # flat-bottom harmonic to the sphere centre: 0.5 k (r-r0)^2 for r>r0, else 0
            com_restraint_force = openmm.CustomCentroidBondForce(
                1,
                "0.5*kcom*sd^2; sd=max(0, r-r0com); "
                "r=sqrt((x1-xcc)^2+(y1-ycc)^2+(z1-zcc)^2)",
            )
            com_restraint_force.addGlobalParameter("kcom", com_k_kj)
            com_restraint_force.addGlobalParameter("r0com", com_r0_nm)
            com_restraint_force.addGlobalParameter("xcc", center_nm[0])
            com_restraint_force.addGlobalParameter("ycc", center_nm[1])
            com_restraint_force.addGlobalParameter("zcc", center_nm[2])
            group = com_restraint_force.addGroup(com_atoms)   # mass-weighted centroid (default)
            com_restraint_force.addBond([group], [])
            com_restraint_force.setUsesPeriodicBoundaryConditions(False)
            self._evb_system_scaffold.addForce(com_restraint_force)
            self._register_force_to_group(com_restraint_force, term="scaas_solute_com", state=None)

        config = SCAASConfig(
            center_nm=center_nm,
            solvent_radius_nm=solvent_radius_nm,
            temperature_kelvin=temperature_kelvin,
            shell_radius_nm=shell_radius_nm,
            shell_force_kj_mol_nm2=shell_force_kj_mol_nm2,
            excluded_force_kj_mol_nm2=excluded_force_kj_mol_nm2,
            radial_force_kj_mol_nm2=radial_force_kj_mol_nm2,
            morse_depth_kj_mol=morse_depth_kj_mol,
            morse_width_per_nm=morse_width_per_nm,
            polarisation_enabled=bool(polarisation),
            charge_correction_enabled=bool(charge_correction_enabled),
            polarisation_force_kj_mol_rad2=polarisation_force_kj_mol_rad2,
            pol_update_steps=int(pol_update_steps),
            water_density_per_a3=float(solvent_density_per_a3),
            water_dipole_ea=float(water_dipole_ea),
            dielectric_constant=float(dielectric_constant),
            solvent_residue_names=residue_names,
            radial_wall_type=str(radial_wall_type),
        )
        self._scaas_enable_called = True
        self._scaas_invalidated_by_rebuild = False
        self._scaas_runtime = SCAASRuntime(
            config=config,
            solvent_molecules=water_molecules,
            shells=shells,
            radial_force=radial_force,
            shell_fix_force=shell_fix_force,
            polarisation_force=polarisation_force,
            excluded_atom_indices=normalized_excluded,
            shell_atom_indices=normalized_shell,
            frozen_atom_indices=tuple(),
            backend_name=backend_key,
        )
        self._scaas_backend_controller = build_scaas_backend_controller(backend_key)
        return MappingProxyType(dict(self.get_scaas_diagnostics()))

    def disable_scaas(self, simulation: app.Simulation | None = None) -> None:
        """
        Disable SCAAS terms by setting the shared ``scaas_enable`` parameter to 0.
        """

        runtime = self._scaas_runtime
        if runtime is None:
            return

        for force in (runtime.radial_force, runtime.shell_fix_force, runtime.polarisation_force):
            if force is None:
                continue
            param_index = self._force_global_parameter_index(force, "scaas_enable")
            if param_index is None:
                continue
            force.setGlobalParameterDefaultValue(param_index, 0.0)

        if simulation is not None and self._context_parameter_exists(simulation.context, "scaas_enable"):
            simulation.context.setParameter("scaas_enable", 0.0)

    def attach_scaas(
        self,
        simulation: app.Simulation,
        *,
        initialize: bool = True,
        required: bool = True,
    ) -> None:
        """
        Attach SCAAS runtime updates to an already created EVB ``Simulation``.

        Raises ``RuntimeError`` when there is no SCAAS runtime to attach — i.e. when
        :meth:`enable_scaas` was never called, or was called *before* a force rebuild
        (``solvate_*`` / ``build_reactive_forces``) discarded its forces. Without this guard the
        call is a silent no-op and the run proceeds with **no boundary potential at all**: the
        droplet evaporates, the solute desolvates, and nothing in the output says so. Pass
        ``required=False`` for pipelines that legitimately run without SCAAS (e.g. periodic PME).
        """

        if self._scaas_runtime is None and required:
            if getattr(self, "_scaas_invalidated_by_rebuild", False):
                raise RuntimeError(
                    "SCAAS was enabled but its forces were then discarded by a force rebuild "
                    "(solvate_* / build_reactive_forces reset the scaffold). Call enable_scaas(...) "
                    "AFTER the last force-building step, then attach_scaas(). "
                    "Pass required=False if this run is intentionally without SCAAS."
                )
            raise RuntimeError(
                "attach_scaas() has nothing to attach: enable_scaas(...) was never called on this "
                "EVBSystem, so the droplet has NO boundary potential (solvate_scaas_sphere_states "
                "only adds the waters; it does not enable the boundary). Call "
                "enable_scaas(solvent_radius_angstrom=..., center=...) after build_reactive_forces(), "
                "then attach_scaas(). Pass required=False if this run is intentionally without SCAAS."
            )

        self._attach_scaas_to_simulation(simulation, initialize=initialize)

    def get_scaas_diagnostics(self, simulation: app.Simulation | None = None) -> Mapping[str, object]:
        runtime = self._scaas_runtime
        if runtime is None:
            return MappingProxyType({"enabled": False})

        backend = self._scaas_backend_controller
        backend_diag = {} if backend is None else dict(backend.diagnostics(simulation))
        payload = dict(summarize_runtime(runtime, backend_diag))
        return MappingProxyType(payload)

    def _attach_scaas_to_simulation(self, simulation: app.Simulation, *, initialize: bool = True) -> None:
        runtime = self._scaas_runtime
        backend = self._scaas_backend_controller
        if runtime is None or backend is None:
            return
        backend.attach(
            simulation,
            runtime=runtime,
            evb_scaffold=self._evb_system_scaffold,
            initialize=initialize,
        )

    @staticmethod
    def _force_global_parameter_index(force: object, parameter_name: str) -> int | None:
        get_count = getattr(force, "getNumGlobalParameters", None)
        get_name = getattr(force, "getGlobalParameterName", None)
        if get_count is None or get_name is None:
            return None
        try:
            count = int(get_count())
        except (TypeError, openmm.OpenMMException):
            return None
        for idx in range(count):
            if str(get_name(idx)) == parameter_name:
                return idx
        return None

    @staticmethod
    def _normalize_atom_index_sequence(
        values: Sequence[int] | None,
        *,
        expected_particles: int,
    ) -> Tuple[int, ...]:
        if values is None:
            return tuple()
        normalized: list[int] = []
        seen: Set[int] = set()
        for raw in values:
            idx = int(raw)
            if idx < 0 or idx >= expected_particles:
                raise ValueError(
                    f"Atom index {idx} is outside the valid range [0, {expected_particles - 1}]."
                )
            if idx in seen:
                continue
            seen.add(idx)
            normalized.append(idx)
        return tuple(sorted(normalized))

    def _default_solute_restraint_atom_indices(
        self,
        *,
        topology: app.Topology,
        system: openmm.System,
        include_hydrogens: bool,
        solvent_residue_names: Sequence[str] | None,
    ) -> Tuple[int, ...]:
        solvent_names = set(normalize_residue_names(solvent_residue_names))
        selected: list[int] = []
        for atom in topology.atoms():
            atom_index = int(atom.index)
            if atom_index < 0 or atom_index >= system.getNumParticles():
                continue

            residue_name = (atom.residue.name or "").strip().upper()
            if residue_name in solvent_names:
                continue

            mass_da = float(system.getParticleMass(atom_index).value_in_unit(unit.dalton))
            if mass_da <= 0.0:
                continue

            if not include_hydrogens and self._is_hydrogen_like_atom(atom):
                continue

            selected.append(atom_index)

        return tuple(sorted(selected))

    def _create_scaas_water_list(
        self,
        *,
        state_name: str,
        residue_names: Sequence[str],
        oxygen_name: str | None,
        hydrogen_names: Sequence[str] | None,
    ) -> Tuple[SCAASWaterMolecule, ...]:
        topology = self._modellers[state_name].topology
        residue_set = {name.strip().upper() for name in residue_names if str(name).strip()}
        oxygen_token = None if oxygen_name is None else oxygen_name.strip().upper()
        hydrogen_tokens: Tuple[str, ...]
        if hydrogen_names is None:
            hydrogen_tokens = tuple()
        else:
            hydrogen_tokens = tuple(name.strip().upper() for name in hydrogen_names if str(name).strip())

        molecules: list[SCAASWaterMolecule] = []
        for residue in topology.residues():
            residue_name = (residue.name or "").strip().upper()
            if residue_name not in residue_set:
                continue
            atoms = tuple(residue.atoms())
            if len(atoms) < 3:
                continue

            oxygen_atom = None
            if oxygen_token is not None:
                for atom in atoms:
                    if atom.name.strip().upper() == oxygen_token:
                        oxygen_atom = atom
                        break
            if oxygen_atom is None:
                oxygen_candidates = [
                    atom
                    for atom in atoms
                    if atom.element is not None and atom.element.atomic_number == 8
                ]
                if oxygen_candidates:
                    oxygen_atom = oxygen_candidates[0]
            if oxygen_atom is None:
                continue

            hydrogens: list[app.Topology.Atom] = []
            if hydrogen_tokens:
                for token in hydrogen_tokens:
                    for atom in atoms:
                        if atom.index == oxygen_atom.index:
                            continue
                        if atom.name.strip().upper() == token:
                            hydrogens.append(atom)
                            break
                    if len(hydrogens) == 2:
                        break
            if len(hydrogens) < 2:
                hydrogen_candidates = [
                    atom
                    for atom in atoms
                    if atom.index != oxygen_atom.index
                    and atom.element is not None
                    and atom.element.atomic_number == 1
                ]
                for atom in hydrogen_candidates:
                    if atom in hydrogens:
                        continue
                    hydrogens.append(atom)
                    if len(hydrogens) == 2:
                        break
            if len(hydrogens) != 2:
                continue

            molecules.append(
                SCAASWaterMolecule(
                    oxygen_index=oxygen_atom.index,
                    hydrogen1_index=hydrogens[0].index,
                    hydrogen2_index=hydrogens[1].index,
                    residue_name=residue_name,
                    residue_id=self._residue_sequence_number(residue, state_name),
                )
            )
        return tuple(molecules)

    def _estimate_scaas_center_nm(self, water_molecules: Sequence[SCAASWaterMolecule]) -> Tuple[float, float, float]:
        if self._evb_system_source_state is None:
            raise ValueError("EVB source state is undefined.")
        positions = self.evb_initial_positions
        topology = self._modellers[self._evb_system_source_state].topology
        solvent_atoms = {
            idx
            for molecule in water_molecules
            for idx in (molecule.oxygen_index, molecule.hydrogen1_index, molecule.hydrogen2_index)
        }

        selected_indices = [atom.index for atom in topology.atoms() if atom.index not in solvent_atoms]
        if not selected_indices:
            selected_indices = [atom.index for atom in topology.atoms()]
        if not selected_indices:
            raise ValueError("Cannot estimate SCAAS center because no atoms are available in the source topology.")

        sx = 0.0
        sy = 0.0
        sz = 0.0
        for idx in selected_indices:
            vec = positions[idx]
            sx += float(vec.x)
            sy += float(vec.y)
            sz += float(vec.z)
        inv_n = 1.0 / float(len(selected_indices))
        return (sx * inv_n, sy * inv_n, sz * inv_n)

    def _estimate_scaas_solute_charge(
        self,
        *,
        state_name: str,
        solvent_atom_indices: Set[int],
    ) -> float:
        systems = self._get_state_systems()
        system = systems[state_name]
        nonbonded = self._find_nonbonded_force(system)
        total_charge = 0.0
        for idx in range(nonbonded.getNumParticles()):
            if idx in solvent_atom_indices:
                continue
            charge, _, _ = nonbonded.getParticleParameters(idx)
            total_charge += float(charge.value_in_unit(unit.elementary_charge))
        return total_charge

    _warn_if_scaas_radius_looks_like_nm = staticmethod(_units.warn_if_scaas_radius_looks_like_nm)
    _warn_if_center_unitless = staticmethod(_units.warn_if_center_unitless)

    @staticmethod
    def _normalize_solvation_model(model: str) -> str:
        if not isinstance(model, str):
            raise TypeError("'model' must be provided as a string.")
        normalized = model.strip().lower()
        if not normalized:
            raise ValueError("'model' must be a non-empty string.")
        supported_models = ("tip3p", "spce")
        if normalized not in supported_models:
            supported_text = ", ".join(supported_models)
            raise ValueError(
                f"Water model '{model}' is not supported by the current force-field setup. "
                f"Use one of: {supported_text}."
            )
        return normalized

    def solvate_states(
        self,
        *,
        reference_state: str | None = None,
        model: str = "tip3p",
        padding: unit.Quantity | float | None = 1.0,
        box_size: unit.Quantity | Sequence[float] | None = None,
        box_shape: str = "cube",
        ionic_strength: unit.Quantity | float | None = None,
        neutralize: bool = False,
        positive_ion: str = "Na+",
        negative_ion: str = "Cl-",
        num_added: int | None = None,
        sphere_radius_angstrom: float | None = None,
        sphere_center: unit.Quantity | Sequence[float] | openmm.Vec3 | None = None,
        solvent_residue_names: Sequence[str] | None = None,
        rebuild_forces: bool = True,
        preserve_solute_parameters: bool = True,
    ) -> Mapping[str, object]:
        """
        Solvate all EVB states consistently and refresh EVB scaffold bookkeeping.

        Notes
        -----
        Solvent molecules are generated once on ``reference_state`` and then appended
        to all other states with identical ordering so EVB atom mapping remains valid.

        When ``preserve_solute_parameters=True`` (default), per-state
        ``NonbondedForce`` particle parameters (charge, sigma, epsilon) for the
        pre-solvation solute atoms — plus any solute-solute exceptions — are
        snapshotted before the internal state-system rebuild and reapplied
        afterwards. This lets callers apply parameter overrides (e.g.
        forcefield-specific charges) before solvation without having them
        silently discarded by the solvent re-parameterization. Overrides on
        non-NonbondedForce forces are NOT preserved; apply those after
        solvation.
        """

        if reference_state is None:
            reference_state = self._evb_system_source_state or self.state_names[0]
        if reference_state not in self._states:
            raise ValueError(
                f"Unknown reference state '{reference_state}'. Available states: {self.state_names!r}."
            )
        if box_size is not None and padding is not None:
            raise ValueError("Specify only one of 'box_size' or 'padding'.")
        if num_added is not None and int(num_added) <= 0:
            raise ValueError("'num_added' must be > 0 when provided.")

        allowed_shapes = ("cube", "dodecahedron", "octahedron")
        normalized_shape = str(box_shape).strip().lower()
        if normalized_shape not in allowed_shapes:
            raise ValueError(
                f"box_shape '{box_shape}' is not supported; use one of {allowed_shapes}."
            )
        if box_size is not None and normalized_shape != "cube":
            # OpenMM's addSolvent only honours boxShape with padding/numAdded; an explicit
            # box_size gives a rectangular box, so a non-cube shape would be silently ignored.
            raise ValueError(
                "box_shape is only honoured with 'padding' (or 'num_added'); it cannot be "
                "combined with an explicit 'box_size'."
            )

        if sphere_radius_angstrom is not None and float(sphere_radius_angstrom) <= 0.0:
            raise ValueError("'sphere_radius_angstrom' must be > 0 when provided.")

        normalized_model = self._normalize_solvation_model(model)
        # Record the water model so both the addSolvent force field (cloned just below)
        # and every post-solvation per-state createSystem load the matching water template.
        self._solvation_water_model = normalized_model

        pre_counts = {state: len(tuple(self._modellers[state].positions)) for state in self.state_names}
        ref_pre_count = pre_counts[reference_state]

        # Preserve which EVB force categories were active before topology changes.
        had_vdw = (
            self._shared_nonbonded_force is not None
            or self._shared_nonbonded_qm_force is not None
            or self._shared_nonbonded_exception_force is not None
            or self._shared_nonbonded_qm_exception_force is not None
            or bool(self._state_specific_nonbonded_forces)
            or bool(self._state_specific_nonbonded_qm_forces)
            or bool(self._state_specific_vdw_exception_forces)
            or bool(self._state_specific_vdw_qm_exception_forces)
            or bool(self._state_specific_electrostatics_exception_forces)
            or bool(self._state_specific_electrostatics_qm_exception_forces)
        )
        previous_mm_cutoff_nm = self._mm_nonbonded_cutoff_nm
        previous_qm_cutoff_nm = self._qm_nonbonded_cutoff_nm
        had_bond = self._shared_bond_force is not None or bool(self._state_specific_bond_forces)
        had_angle = self._shared_angle_force is not None or bool(self._state_specific_angle_forces)
        had_torsion = self._shared_torsion_force is not None or bool(self._state_specific_torsion_forces)

        reference_initial_positions = self._coerce_vec3_tuple(
            self._modellers[reference_state].positions,
            quantity_unit=unit.nanometer,
            label=f"{reference_state} positions",
        )
        reference_modeller = app.Modeller(
            self._modellers[reference_state].topology,
            self._vec3s_to_quantity(reference_initial_positions),
        )
        solvent_force_field = self._clone_force_field_with_generators()

        add_kwargs: Dict[str, object] = {
            "model": normalized_model,
            "neutralize": bool(neutralize),
            "positiveIon": str(positive_ion),
            "negativeIon": str(negative_ion),
        }
        if box_size is not None:
            add_kwargs["boxSize"] = _units.coerce_box_size(box_size)
        elif padding is not None:
            add_kwargs["padding"] = _units.coerce_length_nm(padding, label="padding")
        if box_size is None:
            # boxShape applies to the padding / numAdded paths (a non-orthorhombic box
            # for the same solute-wall distance uses ~29% fewer waters at octahedron).
            add_kwargs["boxShape"] = normalized_shape
        if ionic_strength is not None:
            add_kwargs["ionicStrength"] = _units.coerce_molarity(ionic_strength, label="ionic_strength")
        if num_added is not None:
            add_kwargs["numAdded"] = int(num_added)

        reference_modeller.addSolvent(solvent_force_field, **add_kwargs)

        if sphere_radius_angstrom is not None:
            solvent_names = set(normalize_residue_names(solvent_residue_names))
            radius_nm = float(sphere_radius_angstrom) * ANGSTROM_TO_NM
            reference_positions = self._coerce_vec3_tuple(
                reference_modeller.positions,
                quantity_unit=unit.nanometer,
                label=f"{reference_state} solvated positions",
            )
            if sphere_center is None:
                center_nm = self._mean_position_nm(
                    positions=reference_positions,
                    indices=tuple(range(ref_pre_count)),
                )
            else:
                center_nm = normalize_scaas_center(sphere_center)
            residues_to_remove: list[app.Topology.Residue] = []
            updated_positions = reference_positions
            for residue in reference_modeller.topology.residues():
                residue_name = (residue.name or "").strip().upper()
                if residue_name not in solvent_names:
                    continue
                oxygen_atom = None
                atoms = tuple(residue.atoms())
                for atom in atoms:
                    if atom.element is not None and atom.element.atomic_number == 8:
                        oxygen_atom = atom
                        break
                if oxygen_atom is None and atoms:
                    oxygen_atom = atoms[0]
                if oxygen_atom is None:
                    continue
                pos = updated_positions[oxygen_atom.index]
                dx = float(pos.x) - center_nm[0]
                dy = float(pos.y) - center_nm[1]
                dz = float(pos.z) - center_nm[2]
                distance = math.sqrt(dx * dx + dy * dy + dz * dz)
                if distance > radius_nm:
                    residues_to_remove.append(residue)
            if residues_to_remove:
                reference_modeller.delete(residues_to_remove)

        ref_post_count = len(tuple(reference_modeller.positions))
        added_atom_count = ref_post_count - ref_pre_count
        if added_atom_count <= 0:
            raise RuntimeError(
                "No solvent atoms were added. Adjust 'padding'/'box_size'/'num_added' or check force-field templates."
            )

        solvent_residues = [
            residue
            for residue in reference_modeller.topology.residues()
            if residue.index >= 0 and any(atom.index >= ref_pre_count for atom in residue.atoms())
        ]
        solvent_topology, solvent_atom_indices = self._copy_residue_block_to_topology(
            reference_modeller.topology, solvent_residues
        )
        reference_positions = self._coerce_vec3_tuple(
            reference_modeller.positions,
            quantity_unit=unit.nanometer,
            label=f"{reference_state} solvated positions",
        )
        solvent_positions = tuple(reference_positions[idx] for idx in solvent_atom_indices)
        solvent_positions_quantity = self._vec3s_to_quantity(solvent_positions)
        reference_periodic_vectors = reference_modeller.topology.getPeriodicBoxVectors()

        updated_modellers: Dict[str, app.Modeller] = {reference_state: reference_modeller}
        for state_name in self.state_names:
            if state_name == reference_state:
                continue
            state_initial_positions = self._coerce_vec3_tuple(
                self._modellers[state_name].positions,
                quantity_unit=unit.nanometer,
                label=f"{state_name} positions",
            )
            modeller = app.Modeller(
                self._modellers[state_name].topology,
                self._vec3s_to_quantity(state_initial_positions),
            )
            modeller.add(solvent_topology, solvent_positions_quantity)
            if reference_periodic_vectors is not None:
                modeller.topology.setPeriodicBoxVectors(reference_periodic_vectors)
            updated_modellers[state_name] = modeller

        for state_name in self.state_names:
            before = pre_counts[state_name]
            after = len(tuple(updated_modellers[state_name].positions))
            if after - before != added_atom_count:
                raise RuntimeError(
                    f"State '{state_name}' received {after - before} new atoms, expected {added_atom_count}."
                )

        # Update EVBState payloads and state attributes.
        refreshed_states: Dict[str, EVBState] = {}
        for state_name in self.state_names:
            previous = self._states[state_name]
            modeller = updated_modellers[state_name]
            positions = self._coerce_vec3_tuple(
                modeller.positions,
                quantity_unit=unit.nanometer,
                label=f"{state_name} solvated positions",
            )
            refreshed = EVBState(
                path=previous.path,
                topology=modeller.topology,
                positions=positions,
            )
            refreshed_states[state_name] = refreshed
            setattr(self, state_name, refreshed)
        self._states = refreshed_states

        # Clear atom-reference caches before remapping new solvent atoms.
        self._state_atom_reference_lists = {}
        self._state_atom_identifier_cache = {}
        self._state_canonical_index_maps = {}
        self._state_particle_mass_lists = {}

        self._modellers = {}
        self._topology_bonds = {}
        self._residue_targets = {}
        self._build_residue_targets()
        self._ensure_residue_name_consistency(include_bonds=True)
        self._missing_residues = self._compute_missing_force_field_residues()

        extended_mapping = dict(self._atom_mapping)
        for offset in range(added_atom_count):
            ref_idx = pre_counts[reference_state] + offset
            ref_atom = self._atom_reference_from_index(reference_state, ref_idx)
            extended_mapping.setdefault(ref_atom, ref_atom)
            for state_name in self.state_names:
                state_idx = pre_counts[state_name] + offset
                state_atom = self._atom_reference_from_index(state_name, state_idx)
                extended_mapping[state_atom] = ref_atom
        self._atom_mapping = extended_mapping
        self._atom_canonical_labels = self._build_atom_canonical_labels()

        solute_nb_snapshot: Dict[str, Dict[str, object]] | None = None
        if preserve_solute_parameters:
            solute_nb_snapshot = self._snapshot_per_state_solute_nb_params(
                solute_counts=pre_counts
            )

        self._state_systems = {}
        self._state_systems_cache = None
        self._state_system_signature = None
        self._state_systems_are_prebuilt = False
        self.bonded_force_comparisons = None
        self._invalidate_partition_caches()

        self._reset_evb_force_handles()
        source_state = self._evb_system_source_state or reference_state
        self.build_evb_system(initial_state=source_state)

        if solute_nb_snapshot:
            # Force state-system rebuild now so we can replay overrides onto
            # the fresh forcefield-parameterized NonbondedForces before any
            # scaffold force below reads them.
            self._get_state_systems()
            self._restore_per_state_solute_nb_params(solute_nb_snapshot)

        if rebuild_forces:
            if had_vdw:
                if previous_mm_cutoff_nm is None:
                    self.add_vdw_forces()
                else:
                    kwargs: Dict[str, object] = {"mm_cutoff_distance": previous_mm_cutoff_nm}
                    if previous_qm_cutoff_nm is not None:
                        kwargs["qm_cutoff_distance"] = previous_qm_cutoff_nm
                    self.add_vdw_forces(**kwargs)
            if had_bond:
                self.add_bond_forces()
            if had_angle:
                self.add_angle_forces()
            if had_torsion:
                self.add_torsion_forces()

        return MappingProxyType(
            {
                "reference_state": reference_state,
                "added_atom_count_per_state": int(added_atom_count),
                "added_residue_count_reference": len(solvent_residues),
                "state_atom_counts": MappingProxyType(
                    {
                        state: len(tuple(self._modellers[state].positions))
                        for state in self.state_names
                    }
                ),
            }
        )

    def solvate_scaas_sphere_states(
        self,
        *,
        solvent_radius_angstrom: float,
        reference_state: str | None = None,
        center: unit.Quantity | Sequence[float] | openmm.Vec3 | None = None,
        model: str = "tip3p",
        sphere_buffer_angstrom: float = 2.0,
        neutralize: bool = False,
        ionic_strength: unit.Quantity | float | None = None,
        positive_ion: str = "Na+",
        negative_ion: str = "Cl-",
        solvent_residue_names: Sequence[str] | None = None,
        rebuild_forces: bool = True,
        preserve_solute_parameters: bool = True,
    ) -> Mapping[str, object]:
        """
        Solvate all EVB states using a spherical droplet that matches the SCAAS radius.

        Notes
        -----
        This v1 implementation follows a safety-first policy and rejects ion placement
        options because trimming from a temporary cubic build can distort final
        concentration/neutralization targets.
        """

        radius_a = float(solvent_radius_angstrom)
        if radius_a <= 0.0:
            raise ValueError("'solvent_radius_angstrom' must be > 0.")
        _units.warn_if_scaas_radius_looks_like_nm(radius_a, label="solvent_radius_angstrom")

        buffer_a = float(sphere_buffer_angstrom)
        if buffer_a < 0.0:
            raise ValueError("'sphere_buffer_angstrom' must be >= 0.")

        if neutralize:
            raise ValueError(
                "Sphere solvation v1 does not support 'neutralize=True'. "
                "Use neutralize=False to avoid post-trim ion inconsistencies."
            )
        if ionic_strength is not None:
            raise ValueError(
                "Sphere solvation v1 does not support 'ionic_strength'. "
                "Use ionic_strength=None to avoid post-trim concentration inconsistencies."
            )

        if reference_state is None:
            reference_state = self._evb_system_source_state or self.state_names[0]
        if reference_state not in self._states:
            raise ValueError(
                f"Unknown reference state '{reference_state}'. Available states: {self.state_names!r}."
            )

        reference_positions = self._coerce_vec3_tuple(
            self._modellers[reference_state].positions,
            quantity_unit=unit.nanometer,
            label=f"{reference_state} positions",
        )
        if center is None:
            center_nm = self._mean_position_nm(
                positions=reference_positions,
                indices=tuple(range(len(reference_positions))),
            )
        else:
            _units.warn_if_center_unitless(center, context="solvate_scaas_sphere_states")
            center_nm = normalize_scaas_center(center)

        radius_nm = radius_a * ANGSTROM_TO_NM
        buffer_nm = buffer_a * ANGSTROM_TO_NM
        margin_nm = radius_nm + buffer_nm
        # addSolvent centres the temporary solvation box on the SOLUTE bounding-box
        # centre, but the droplet sphere is carved around `center_nm`. When the solute
        # is larger than / off-centre to the sphere (e.g. an enzyme whose reactive
        # substrate sits off-centre in the protein), a cubic box of side 2*(radius+buffer)
        # leaves the sphere poking outside the water on one face -> a dry cap and an
        # under-solvated droplet. Size each axis so the full sphere around `center_nm`
        # stays inside the box after solute-centring:
        #     half_side_axis >= |center - solute_bbox_centre| + (radius + buffer).
        _sx = [float(p.x) for p in reference_positions]
        _sy = [float(p.y) for p in reference_positions]
        _sz = [float(p.z) for p in reference_positions]
        _bbox_center = (
            0.5 * (min(_sx) + max(_sx)),
            0.5 * (min(_sy) + max(_sy)),
            0.5 * (min(_sz) + max(_sz)),
        )
        _half = tuple(
            max(margin_nm, abs(float(center_nm[k]) - _bbox_center[k]) + margin_nm)
            for k in range(3)
        )
        # Self-check the invariant the sizing establishes (guards against regressions):
        for k in range(3):
            assert _half[k] + 1e-9 >= abs(float(center_nm[k]) - _bbox_center[k]) + margin_nm, (
                "solvation box does not cover the droplet sphere on axis "
                f"{k}: half={_half[k]:.3f} nm, need "
                f"{abs(float(center_nm[k]) - _bbox_center[k]) + margin_nm:.3f} nm"
            )
        box_size = openmm.Vec3(2.0 * _half[0], 2.0 * _half[1], 2.0 * _half[2]) * unit.nanometer

        self._boundary_mode = "scaas_sphere"
        base_summary = self.solvate_states(
            reference_state=reference_state,
            model=model,
            padding=None,
            box_size=box_size,
            ionic_strength=None,
            neutralize=False,
            positive_ion=positive_ion,
            negative_ion=negative_ion,
            num_added=None,
            sphere_radius_angstrom=radius_a,
            sphere_center=center_nm,
            solvent_residue_names=solvent_residue_names,
            rebuild_forces=rebuild_forces,
            preserve_solute_parameters=preserve_solute_parameters,
        )

        cleared_states: list[str] = []
        for state_name in self.state_names:
            topology = self._modellers[state_name].topology
            if topology.getPeriodicBoxVectors() is not None:
                topology.setPeriodicBoxVectors(None)
                cleared_states.append(state_name)
            state_topology = self._states[state_name].topology
            if state_topology.getPeriodicBoxVectors() is not None:
                state_topology.setPeriodicBoxVectors(None)

        payload = dict(base_summary)
        payload.update(
            {
                "solvation_mode": "scaas_sphere",
                "sphere_radius_angstrom": radius_a,
                "sphere_radius_nm": radius_nm,
                "sphere_center_nm": center_nm,
                "sphere_buffer_angstrom": buffer_a,
                "temporary_box_size_nm": (2.0 * _half[0], 2.0 * _half[1], 2.0 * _half[2]),
                "periodic_box_cleared_states": tuple(cleared_states),
            }
        )
        return MappingProxyType(payload)

    _WATER_RESNAMES = frozenset({"HOH", "WAT", "TIP3", "TIP", "H2O", "SOL", "T3P"})
    _DEFAULT_STRUCTURAL_RESNAMES = (
        "CA", "CAL", "ZN", "ZN2", "MN", "MN2", "FE", "FE2", "FE3",
        "CU", "CU1", "NI", "CO", "MG_STRUCT",
    )
    _DEFAULT_ION_RESNAMES = (
        "NA", "SOD", "NA+", "CL", "CLA", "CL-", "K", "POT", "K+",
        "LI", "RB", "CS", "BR", "IOD", "IB", "MG",
    )

    def _build_water_block(self, positions_nm, water_model="tip3p"):
        """A Topology block of one HOH per given O-position, at TIP3P geometry (O-H 0.9572 A,
        H-O-H 104.52 deg), placed so the standard minimize relaxes any overlap. Names O/H1/H2 +
        resname HOH match the tip3p template so the force field parameterizes them on rebuild."""
        self._normalize_solvation_model(water_model)          # validate model name
        top = app.Topology()
        chain = top.addChain()
        d, half = 0.09572, math.radians(104.52 / 2.0)
        coords = []
        for p in positions_nm:
            res = top.addResidue("HOH", chain)
            aO = top.addAtom("O", app.element.oxygen, res)
            aH1 = top.addAtom("H1", app.element.hydrogen, res)
            aH2 = top.addAtom("H2", app.element.hydrogen, res)
            top.addBond(aO, aH1)
            top.addBond(aO, aH2)
            ox, oy, oz = float(p.x), float(p.y), float(p.z)
            coords.append(openmm.Vec3(ox, oy, oz))
            coords.append(openmm.Vec3(ox + d * math.cos(half), oy + d * math.sin(half), oz))
            coords.append(openmm.Vec3(ox + d * math.cos(half), oy - d * math.sin(half), oz))
        return top, self._vec3s_to_quantity(tuple(coords))

    def carve_scaas_sphere_states(
        self,
        *,
        solvent_radius_angstrom: float,
        center: unit.Quantity | Sequence[float] | openmm.Vec3 | None = None,
        reference_state: str | None = None,
        structural_resnames: Sequence[str] = _DEFAULT_STRUCTURAL_RESNAMES,
        ion_resnames: Sequence[str] = _DEFAULT_ION_RESNAMES,
        replace_ions_with_water: bool = True,
        water_model: str = "tip3p",
        rebuild_forces: bool = True,
        preserve_solute_parameters: bool = True,
    ) -> Mapping[str, object]:
        """Carve a SCAAS droplet out of an ALREADY-SOLVATED multi-state system (subtractive sibling
        of :meth:`solvate_scaas_sphere_states`, which is additive).

        Use this to extract a spherical sub-system from a large periodic frame (e.g. a GS-MD snapshot)
        while KEEPING the equilibrated water inside R — the whole point of the GS run. The operation is
        applied IDENTICALLY to every diabatic state (the water/ion block is shared and identically
        ordered across states), so the EVB atom-mapping bijection is preserved.

        Policy:
          * the whole solute + every residue in ``structural_resnames`` (catalytic metals) are ALWAYS
            kept — no solute truncation;
          * water whose O-atom is within R of ``center`` is kept; water outside R is dropped;
          * each NON-structural ion (``ion_resnames``) inside R is REPLACED with a single water at the
            ion's position (``replace_ions_with_water``; then the standard minimize relaxes the overlap),
            so snapshots do not differ by a stray buried counter-ion; ions outside R are dropped.

        ``center`` defaults to the reference state's solute centroid. Restraints on / neutralization of
        everything outside R are separate: SCAAS (:meth:`enable_scaas`) supplies the boundary+restraints;
        :meth:`neutralize_outside_sphere` (opt-in) handles periphery charge if wanted.

        ORDERING — call this BEFORE configuring the reactive core, exactly like
        :meth:`solvate_scaas_sphere_states`. The carve rebuilds the per-state Systems from the force
        field and RENUMBERS residues, so it deliberately does NOT touch QM atoms or reactive forces
        (both are name/index-based and would be invalidated). After carving, (re)configure the reactive
        core on the carved system: ``promote_reactive_constraints_to_bonds`` on each per-state System
        (the reactive X–H is re-constrained by ``constraints=HBonds`` on rebuild), then
        ``set_qm_atoms`` and ``build_reactive_forces`` — the same steps you run on a freshly solvated
        sphere.
        """
        radius_a = float(solvent_radius_angstrom)
        if radius_a <= 0.0:
            raise ValueError("'solvent_radius_angstrom' must be > 0.")
        _units.warn_if_scaas_radius_looks_like_nm(radius_a, label="solvent_radius_angstrom")
        radius_nm = radius_a * ANGSTROM_TO_NM

        if reference_state is None:
            reference_state = self._evb_system_source_state or self.state_names[0]
        if reference_state not in self._states:
            raise ValueError(
                f"Unknown reference state '{reference_state}'. Available states: {self.state_names!r}."
            )

        water_names = set(self._WATER_RESNAMES)
        struct_names = {str(n).strip().upper() for n in structural_resnames}
        ion_names = {str(n).strip().upper() for n in ion_resnames} - struct_names

        ref_positions = self._coerce_vec3_tuple(
            self._modellers[reference_state].positions,
            quantity_unit=unit.nanometer, label=f"{reference_state} positions",
        )
        if center is None:
            solute_idx = tuple(
                a.index for a in self._modellers[reference_state].topology.atoms()
                if (a.residue.name or "").strip().upper() not in (water_names | ion_names)
            )
            if not solute_idx:
                raise RuntimeError("carve: no solute atoms found to centre the sphere on.")
            center_nm = self._mean_position_nm(positions=ref_positions, indices=solute_idx)
        else:
            _units.warn_if_center_unitless(center, context="carve_scaas_sphere_states")
            center_nm = normalize_scaas_center(center)
        cx, cy, cz = float(center_nm[0]), float(center_nm[1]), float(center_nm[2])
        r2 = radius_nm * radius_nm

        def _o_pos(res):
            for a in res.atoms():
                if a.element is not None and a.element.atomic_number == 8:
                    return ref_positions[a.index]
            return ref_positions[next(res.atoms()).index]

        # Decide on the reference state (shared water/ion block => identical residue indices in all states).
        delete_res_indices: set[int] = set()
        replace_ion_positions: list = []
        n_water_dropped = n_ions_removed = 0
        for res in self._modellers[reference_state].topology.residues():
            rn = (res.name or "").strip().upper()
            if rn in water_names:
                p = _o_pos(res)
                if (float(p.x) - cx) ** 2 + (float(p.y) - cy) ** 2 + (float(p.z) - cz) ** 2 > r2:
                    delete_res_indices.add(res.index)
                    n_water_dropped += 1
            elif rn in ion_names:
                a0 = next(res.atoms())
                p = ref_positions[a0.index]
                inside = (float(p.x) - cx) ** 2 + (float(p.y) - cy) ** 2 + (float(p.z) - cz) ** 2 <= r2
                delete_res_indices.add(res.index)
                n_ions_removed += 1
                if inside and replace_ions_with_water:
                    replace_ion_positions.append(openmm.Vec3(float(p.x), float(p.y), float(p.z)))
        n_ions_replaced = len(replace_ion_positions)

        replacement_block = replacement_positions_q = None
        if replace_ion_positions:
            replacement_block, replacement_positions_q = self._build_water_block(
                replace_ion_positions, water_model
            )

        # Capture active force categories + cutoffs so we can rebuild them (as solvate_states does).
        had_vdw = (
            self._shared_nonbonded_force is not None
            or self._shared_nonbonded_qm_force is not None
            or bool(self._state_specific_nonbonded_forces)
            or bool(self._state_specific_nonbonded_qm_forces)
        )
        previous_mm_cutoff_nm = self._mm_nonbonded_cutoff_nm
        previous_qm_cutoff_nm = self._qm_nonbonded_cutoff_nm
        had_bond = self._shared_bond_force is not None or bool(self._state_specific_bond_forces)
        had_angle = self._shared_angle_force is not None or bool(self._state_specific_angle_forces)
        had_torsion = self._shared_torsion_force is not None or bool(self._state_specific_torsion_forces)

        pre_counts = {s: len(tuple(self._modellers[s].positions)) for s in self.state_names}
        # resname per AtomReference on the ORIGINAL topology -> preserve solute (incl reactive) mapping.
        orig_resname: Dict[AtomReference, str] = {}
        for s in self.state_names:
            for a in self._modellers[s].topology.atoms():
                orig_resname[(s, a.residue.chain.id, int(a.residue.id), a.name)] = (
                    a.residue.name or ""
                ).strip().upper()

        updated_modellers: Dict[str, app.Modeller] = {}
        for s in self.state_names:
            positions_q = self._vec3s_to_quantity(
                self._coerce_vec3_tuple(
                    self._modellers[s].positions, quantity_unit=unit.nanometer,
                    label=f"{s} positions",
                )
            )
            modeller = app.Modeller(self._modellers[s].topology, positions_q)
            to_delete = [r for r in modeller.topology.residues() if r.index in delete_res_indices]
            # Guard: the solvent block must be identically ordered across states, so every residue at a
            # shared delete index must be water / non-structural ion in THIS state too. Fail loudly on a
            # misalignment rather than silently deleting a solute or structural residue.
            for r in to_delete:
                rn = (r.name or "").strip().upper()
                if rn not in water_names and rn not in ion_names:
                    raise RuntimeError(
                        f"carve: residue index {r.index} ('{rn}') in state '{s}' is not water/ion; the "
                        "solvent block is not aligned across states -- cannot carve by shared residue index."
                    )
            if to_delete:
                modeller.delete(to_delete)
            if replacement_block is not None:
                modeller.add(replacement_block, replacement_positions_q)
            modeller.topology.setPeriodicBoxVectors(None)
            updated_modellers[s] = modeller

        counts = {s: len(tuple(updated_modellers[s].positions)) for s in self.state_names}
        if len(set(counts.values())) != 1:
            raise RuntimeError(f"carve produced unequal per-state atom counts: {counts}")

        # Rebuild atom mapping: keep solute (non-water/ion) entries verbatim (preserves the reactive
        # core swap); map water by INDEX correspondence to the reference state. The carve applies the
        # identical delete+add to every state, so atom index i is the same water in all states even if
        # Modeller.add assigned a different chain/resSeq -- index correspondence is robust, name identity
        # is not.
        new_mapping: Dict[AtomReference, AtomReference] = {
            k: v for k, v in self._atom_mapping.items()
            if orig_resname.get(k, "") not in (water_names | ion_names)
        }
        ref_atoms = list(updated_modellers[reference_state].topology.atoms())
        for s in self.state_names:
            for i, a in enumerate(updated_modellers[s].topology.atoms()):
                if (a.residue.name or "").strip().upper() in water_names:
                    ra = ref_atoms[i]
                    key = (s, a.residue.chain.id, int(a.residue.id), a.name)
                    new_mapping[key] = (
                        reference_state, ra.residue.chain.id, int(ra.residue.id), ra.name
                    )

        # Refresh EVBState + attributes (mirror solvate_states).
        refreshed_states: Dict[str, EVBState] = {}
        for state_name in self.state_names:
            previous = self._states[state_name]
            modeller = updated_modellers[state_name]
            positions = self._coerce_vec3_tuple(
                modeller.positions, quantity_unit=unit.nanometer,
                label=f"{state_name} carved positions",
            )
            refreshed = EVBState(path=previous.path, topology=modeller.topology, positions=positions)
            refreshed_states[state_name] = refreshed
            setattr(self, state_name, refreshed)
        self._states = refreshed_states

        self._state_atom_reference_lists = {}
        self._state_atom_identifier_cache = {}
        self._state_canonical_index_maps = {}
        self._state_particle_mass_lists = {}
        self._modellers = {}
        self._topology_bonds = {}
        self._residue_targets = {}
        self._build_residue_targets()
        self._ensure_residue_name_consistency(include_bonds=True)
        self._missing_residues = self._compute_missing_force_field_residues()

        self._atom_mapping = new_mapping
        self._atom_canonical_labels = self._build_atom_canonical_labels()

        solute_nb_snapshot = None
        if preserve_solute_parameters:
            surviving_solute = {
                s: sum(
                    1 for a in updated_modellers[s].topology.atoms()
                    if (a.residue.name or "").strip().upper() not in water_names
                )
                for s in self.state_names
            }
            solute_nb_snapshot = self._snapshot_per_state_solute_nb_params(solute_counts=surviving_solute)

        self._state_systems = {}
        self._state_systems_cache = None
        self._state_system_signature = None
        self._state_systems_are_prebuilt = False
        self.bonded_force_comparisons = None
        self._invalidate_partition_caches()
        self._reset_evb_force_handles()

        self._boundary_mode = "scaas_sphere"
        source_state = self._evb_system_source_state or reference_state
        self.build_evb_system(initial_state=source_state)

        if solute_nb_snapshot:
            self._get_state_systems()
            self._restore_per_state_solute_nb_params(solute_nb_snapshot)

        # Category-wise force rebuild, IDENTICAL to solvate_scaas_sphere_states/solvate_states: this
        # method deliberately does NOT touch QM atoms or reactive forces (they are name-based and the
        # rebuild renumbers residues). Configure the reactive core AFTER the carve (see docstring).
        if rebuild_forces:
            if had_vdw:
                if previous_mm_cutoff_nm is None:
                    self.add_vdw_forces()
                else:
                    kw: Dict[str, object] = {"mm_cutoff_distance": previous_mm_cutoff_nm}
                    if previous_qm_cutoff_nm is not None:
                        kw["qm_cutoff_distance"] = previous_qm_cutoff_nm
                    self.add_vdw_forces(**kw)
            if had_bond:
                self.add_bond_forces()
            if had_angle:
                self.add_angle_forces()
            if had_torsion:
                self.add_torsion_forces()

        # Bijection guard: every non-reference state atom maps to a valid reference atom.
        ref_refs = {
            (reference_state, a.residue.chain.id, int(a.residue.id), a.name)
            for a in self._states[reference_state].topology.atoms()
        }
        for s in self.state_names:
            if s == reference_state:
                continue
            for a in self._states[s].topology.atoms():
                key = (s, a.residue.chain.id, int(a.residue.id), a.name)
                tgt = self._atom_mapping.get(key)
                if tgt is None or tgt not in ref_refs:
                    raise RuntimeError(
                        f"carve broke the atom mapping bijection: {key} -> {tgt} not a reference atom."
                    )

        return MappingProxyType(
            {
                "solvation_mode": "scaas_sphere_carved",
                "sphere_radius_angstrom": radius_a,
                "sphere_radius_nm": radius_nm,
                "sphere_center_nm": (cx, cy, cz),
                "n_water_dropped": int(n_water_dropped),
                "n_ions_removed": int(n_ions_removed),
                "n_ions_replaced_with_water": int(n_ions_replaced),
                "state_atom_counts": MappingProxyType(counts),
                "pre_carve_atom_counts": MappingProxyType(pre_counts),
            }
        )

    def neutralize_outside_sphere(
        self,
        *,
        solvent_radius_angstrom: float,
        center: unit.Quantity | Sequence[float] | openmm.Vec3 | None = None,
        reference_state: str | None = None,
        threshold_e: float = 0.1,
        apply: bool = False,
    ) -> Mapping[str, object]:
        """Redistribute the net monopole of each non-solvent residue whose centroid lies OUTSIDE R,
        onto that residue's outside-R atoms, so the droplet periphery is charge-neutral. Applied
        IDENTICALLY to every state system's ``NonbondedForce`` (cancels in the EVB energy gap).

        OPT-IN: ``apply=False`` (default) only REPORTS what would change and mutates nothing — SCAAS
        Q6 surface polarization already screens the net charge, so this is usually unnecessary for
        small spheres. ``apply=True`` performs the redistribution.
        """
        radius_a = float(solvent_radius_angstrom)
        if radius_a <= 0.0:
            raise ValueError("'solvent_radius_angstrom' must be > 0.")
        _units.warn_if_scaas_radius_looks_like_nm(radius_a, label="solvent_radius_angstrom")
        radius_nm = radius_a * ANGSTROM_TO_NM
        if reference_state is None:
            reference_state = self._evb_system_source_state or self.state_names[0]
        if reference_state not in self._states:
            raise ValueError(
                f"Unknown reference state '{reference_state}'. Available states: {self.state_names!r}."
            )

        water_names = set(self._WATER_RESNAMES) | {str(n).strip().upper() for n in self._DEFAULT_ION_RESNAMES}
        struct_names = {str(n).strip().upper() for n in self._DEFAULT_STRUCTURAL_RESNAMES}
        positions = self._coerce_vec3_tuple(
            self._states[reference_state].positions, quantity_unit=unit.nanometer,
            label=f"{reference_state} positions",
        )
        if center is None:
            solute_idx = tuple(
                a.index for a in self._states[reference_state].topology.atoms()
                if (a.residue.name or "").strip().upper() not in water_names
            )
            center_nm = self._mean_position_nm(positions=positions, indices=solute_idx)
        else:
            _units.warn_if_center_unitless(center, context="neutralize_outside_sphere")
            center_nm = normalize_scaas_center(center)
        cx, cy, cz = float(center_nm[0]), float(center_nm[1]), float(center_nm[2])
        r2 = radius_nm * radius_nm

        systems = self._get_state_systems()

        def _nb(system):
            for f in system.getForces():
                if isinstance(f, openmm.NonbondedForce):
                    return f
            raise RuntimeError("state system has no NonbondedForce")

        e = unit.elementary_charge
        total_removed = 0.0
        n_res = 0
        # Selection on the reference state's topology; residue indices are shared across states.
        for res in self._states[reference_state].topology.residues():
            rn = (res.name or "").strip().upper()
            if rn in water_names or rn in struct_names:
                continue                                       # never neutralize water/ions/structural metals
            idxs = [a.index for a in res.atoms()]
            cen = [
                sum(float(getattr(positions[i], ax)) for i in idxs) / len(idxs)
                for ax in ("x", "y", "z")
            ]
            if (cen[0] - cx) ** 2 + (cen[1] - cy) ** 2 + (cen[2] - cz) ** 2 <= r2:
                continue                                       # residue centroid inside R -> keep
            nb0 = _nb(systems[reference_state])
            net = float(sum(nb0.getParticleParameters(i)[0].value_in_unit(e) for i in idxs))
            if abs(net) <= threshold_e:
                continue
            d = [
                ((float(positions[i].x) - cx) ** 2 + (float(positions[i].y) - cy) ** 2
                 + (float(positions[i].z) - cz) ** 2)
                for i in idxs
            ]
            outside = [i for i, di in zip(idxs, d) if di > r2] or list(idxs)
            shift = net / len(outside)
            n_res += 1
            total_removed += net
            if apply:
                for state_name in self.state_names:
                    nb = _nb(systems[state_name])
                    for i in outside:
                        q, sig, eps = nb.getParticleParameters(i)
                        nb.setParticleParameters(i, (q.value_in_unit(e) - shift) * e, sig, eps)

        if apply:
            # Parity check: all states carry the same net charge after redistribution.
            totals = []
            for state_name in self.state_names:
                nb = _nb(systems[state_name])
                totals.append(sum(nb.getParticleParameters(i)[0].value_in_unit(e)
                                  for i in range(nb.getNumParticles())))
            if max(totals) - min(totals) > 1e-6:
                raise RuntimeError(f"neutralize_outside_sphere broke charge parity across states: {totals}")

        return MappingProxyType(
            {
                "applied": bool(apply),
                "sphere_radius_angstrom": radius_a,
                "residues_neutralized": int(n_res),
                "charge_removed_e": float(total_removed),
            }
        )

    def solvate_periodic_box_states(
        self,
        *,
        reference_state: str | None = None,
        model: str = "tip3p",
        padding: unit.Quantity | float | None = 1.0,
        box_size: unit.Quantity | Sequence[float] | None = None,
        box_shape: str = "cube",
        neutralize: bool = True,
        ionic_strength: unit.Quantity | float | None = None,
        positive_ion: str = "Na+",
        negative_ion: str = "Cl-",
        nonbonded_cutoff: unit.Quantity | float = 1.0,
        ewald_error_tolerance: float = 5e-4,
        electrostatics_mode: str = "per_state",
        rebuild_forces: bool = True,
        preserve_solute_parameters: bool = True,
    ) -> Mapping[str, object]:
        """
        Solvate all EVB states in a periodic box for PME (Ewald) electrostatics.

        Unlike :meth:`solvate_scaas_sphere_states`, the periodic box vectors are
        *preserved* and stamped onto the EVB scaffold System so that
        ``CutoffPeriodic``/PME forces engage in the ``Context``. Ion placement
        (``neutralize``/``ionic_strength``) is supported and encouraged: PME on a
        net-charged box applies a uniform neutralizing background whose finite-size
        offset is avoided by neutralizing.

        Parameters
        ----------
        nonbonded_cutoff :
            Real-space cutoff for the periodic nonbonded forces. Defaults to 1.0 nm —
            far shorter than the non-periodic SCAAS default (2.0 nm) because PME treats
            long-range electrostatics in reciprocal space. The box must satisfy the
            minimum-image condition (every edge >= 2 x cutoff); a too-large cutoff for
            the chosen padding raises a clear error.
        box_shape :
            Periodic cell shape passed to OpenMM's ``addSolvent`` (``"cube"`` default,
            ``"dodecahedron"``, or ``"octahedron"``). A non-orthorhombic shape encloses the
            same solute-to-wall distance in ~29% less volume (fewer waters, cheaper PME) and
            matches the truncated-octahedron boxes common in GROMACS EVB protocols. Honoured
            only with ``padding``/``num_added``; combining it with an explicit ``box_size``
            raises.
        electrostatics_mode :
            ``"per_state"`` (default) builds one exact PME force per diabatic state;
            ``"mixed"`` builds a single mixed-charge PME plus a sparse diabatic
            correction (faster, one FFT per step).
        """

        if electrostatics_mode not in ("per_state", "mixed"):
            raise ValueError(
                f"Unknown electrostatics_mode {electrostatics_mode!r}; "
                "expected 'per_state' or 'mixed'."
            )

        cutoff_nm = float(
            _units.coerce_length_nm(nonbonded_cutoff, label="nonbonded_cutoff").value_in_unit(
                unit.nanometer
            )
        )

        previous_boundary_mode = self._boundary_mode
        self._boundary_mode = "periodic"
        self._electrostatics_mode = electrostatics_mode
        self._ewald_error_tolerance = float(ewald_error_tolerance)
        # Drive the force-rebuild inside solvate_states to the periodic cutoff (the
        # rebuild otherwise reuses the previous, droplet-sized cutoff and trips the
        # minimum-image guard).
        self._mm_nonbonded_cutoff_nm = cutoff_nm
        self._qm_nonbonded_cutoff_nm = cutoff_nm

        try:
            base_summary = self.solvate_states(
                reference_state=reference_state,
                model=model,
                padding=padding,
                box_size=box_size,
                box_shape=box_shape,
                ionic_strength=ionic_strength,
                neutralize=neutralize,
                positive_ion=positive_ion,
                negative_ion=negative_ion,
                num_added=None,
                sphere_radius_angstrom=None,
                sphere_center=None,
                solvent_residue_names=None,
                rebuild_forces=rebuild_forces,
                preserve_solute_parameters=preserve_solute_parameters,
            )
        except Exception:
            self._boundary_mode = previous_boundary_mode
            raise

        source_state = self._evb_system_source_state or self.state_names[0]
        box_vectors = self._states[source_state].topology.getPeriodicBoxVectors()
        if box_vectors is None:
            self._boundary_mode = previous_boundary_mode
            raise RuntimeError(
                "Periodic solvation did not produce box vectors; check padding/box_size."
            )
        box_edges_nm = tuple(
            float(box_vectors[i][i].value_in_unit(unit.nanometer)) for i in range(3)
        )

        # Fail early if the requested cutoff violates the minimum-image condition,
        # even when no nonbonded forces were rebuilt during this call.
        try:
            self._validate_periodic_cutoff(cutoff_nm)
        except Exception:
            self._boundary_mode = previous_boundary_mode
            raise

        net_charge = self._total_system_charge_e(source_state)

        payload = dict(base_summary)
        payload.update(
            {
                "solvation_mode": "periodic",
                "electrostatics_mode": electrostatics_mode,
                "ewald_error_tolerance": float(ewald_error_tolerance),
                "box_shape": str(box_shape).strip().lower(),
                "box_edges_nm": box_edges_nm,
                "net_charge_e": net_charge,
            }
        )
        return MappingProxyType(payload)

    def _total_system_charge_e(self, state_name: str) -> float:
        """Net charge (elementary units) of a state's reference NonbondedForce."""
        state_system = self._get_state_systems()[state_name]
        nb = self._find_nonbonded_force(state_system)
        total = 0.0
        for i in range(nb.getNumParticles()):
            charge, _, _ = nb.getParticleParameters(i)
            total += float(charge.value_in_unit(unit.elementary_charge))
        return total

    def _reset_evb_force_handles(self) -> None:
        self._shared_bond_force = None
        self._state_specific_bond_forces = {}
        self._shared_angle_force = None
        self._state_specific_angle_forces = {}
        self._shared_torsion_force = None
        self._state_specific_torsion_forces = {}
        self._shared_nonbonded_force = None
        self._shared_nonbonded_qm_force = None
        self._shared_nonbonded_exception_force = None
        self._shared_nonbonded_qm_exception_force = None
        self._state_specific_nonbonded_forces = {}
        self._state_specific_nonbonded_qm_forces = {}
        self._state_specific_vdw_exception_forces = {}
        self._state_specific_vdw_qm_exception_forces = {}
        self._state_specific_reactive_vdw_forces = {}
        self._endpoint_regular_ff_forces = {}
        self._endpoint_regular_ff_groups = {}
        self._state_specific_electrostatics_exception_forces = {}
        self._state_specific_electrostatics_qm_exception_forces = {}
        self._state_specific_pme_forces = {}
        self._state_pme_g_param_names = {}
        self._mm_nonbonded_cutoff_nm = None
        self._qm_nonbonded_cutoff_nm = None
        # A force rebuild drops the SCAAS forces on the floor. Remember that it happened so
        # attach_scaas can say "you enabled SCAAS too early" instead of silently doing nothing.
        if getattr(self, "_scaas_runtime", None) is not None:
            self._scaas_invalidated_by_rebuild = True
        self._scaas_runtime = None
        self._scaas_backend_controller = None
        self._force_group_ids = {}
        self._force_group_definitions = {}
        self._next_force_group_id = 0

    @staticmethod
    def _copy_residue_block_to_topology(
        topology: app.Topology,
        residues: Sequence[app.Topology.Residue],
    ) -> Tuple[app.Topology, Tuple[int, ...]]:
        fragment = app.Topology()
        periodic_vectors = topology.getPeriodicBoxVectors()
        if periodic_vectors is not None:
            fragment.setPeriodicBoxVectors(periodic_vectors)

        chain_map: Dict[int, app.Topology.Chain] = {}
        atom_map: Dict[int, app.Topology.Atom] = {}
        atom_indices: list[int] = []

        for residue in residues:
            chain = residue.chain
            new_chain = chain_map.get(chain.index)
            if new_chain is None:
                new_chain = fragment.addChain(chain.id)
                chain_map[chain.index] = new_chain
            new_residue = fragment.addResidue(
                residue.name,
                new_chain,
                residue.id,
                residue.insertionCode,
            )
            for atom in residue.atoms():
                new_atom = fragment.addAtom(atom.name, atom.element, new_residue, atom.id)
                atom_map[atom.index] = new_atom
                atom_indices.append(atom.index)

        for atom1, atom2 in topology.bonds():
            if atom1.index not in atom_map or atom2.index not in atom_map:
                continue
            fragment.addBond(atom_map[atom1.index], atom_map[atom2.index])

        return fragment, tuple(atom_indices)

    @staticmethod
    def _mean_position_nm(
        *,
        positions: Tuple[openmm.Vec3, ...],
        indices: Sequence[int],
    ) -> Tuple[float, float, float]:
        if not indices:
            raise ValueError("Cannot compute center from an empty atom-index set.")
        sx = 0.0
        sy = 0.0
        sz = 0.0
        for idx in indices:
            vec = positions[int(idx)]
            sx += float(vec.x)
            sy += float(vec.y)
            sz += float(vec.z)
        inv_n = 1.0 / float(len(indices))
        return (sx * inv_n, sy * inv_n, sz * inv_n)

    def validate_reactive_treatment(self, *, raise_on_error: bool = False) -> None:
        """Warn if reactive (Morse) pairs are left with hard Lennard-Jones interactions.

        Standard EVB replaces the forming/breaking pair's Lennard-Jones with a soft
        exponential (Q6 soft-pair) so the atoms can reach bonding distance without a
        steep ``r**-12`` wall. If reactive Morse bonds were installed (e.g. via
        :meth:`replace_reactive_bonds_with_morse`) but no soft-pair vdW force was added
        (:meth:`add_reactive_softpair_vdw_forces`), the reactive atoms keep hard LJ,
        which destabilises the forming-bond region and can prevent the reaction from
        being sampled at all. This is easy to miss on the manual
        ``auto_build_forces=False`` path, where the soft-pair step must be called
        explicitly after :meth:`add_vdw_forces`.

        Emits a :class:`UserWarning` (or raises ``RuntimeError`` when ``raise_on_error``
        is set) if Morse reactive bonds exist without a corresponding soft-pair. This is
        a heuristic — it checks that the soft-pair step ran, not every per-atom LJ term —
        so it can be ignored if bare LJ on the reactive atoms is genuinely intended.
        """
        has_morse = bool(self._state_specific_reactive_morse_bond_forces)
        has_softpair = bool(self._state_specific_reactive_vdw_forces)
        if has_morse and not has_softpair:
            msg = (
                "EVB reactive Morse bonds are installed but no soft-pair vdW force was "
                "added, so the reactive atoms retain hard Lennard-Jones (r**-12) "
                "interactions. This creates a steep repulsive wall near bonding distance "
                "that destabilises the forming-bond region and can prevent the reaction "
                "from being sampled. Call add_reactive_softpair_vdw_forces() after "
                "add_vdw_forces() (or pass softpair_parameters_by_atom_type to the "
                "constructor). If bare LJ on the reactive atoms is intentional, ignore "
                "this warning."
            )
            if raise_on_error:
                raise RuntimeError(msg)
            warnings.warn(msg, stacklevel=2)


    def _assert_qm_cutoff_covers_system(self, qm_cutoff_nm: float, reactive: Sequence[int]) -> None:
        """Raise if any atom lies further from a reactive atom than the QM cutoff.

        Measured 2026-08-13: E(state2) - E(state1) swings 50-100 kcal/mol, non-monotonically, with
        this cutoff. The diabats differ precisely in the reactive-core charges, so a truncation does
        NOT cancel in the gap. A reference run whose cutoff exceeds its droplet and an enzyme run
        whose cutoff falls inside its own droplet are therefore not comparable -- and alpha, fitted
        per system, silently absorbs the difference. Q's own default is 99 A (no cutoff on reacting
        fragments); this library now defaults to 9.9 nm to match.
        """
        if self._boundary_mode == "periodic":
            return
        reactive = sorted({int(i) for i in reactive})
        positions = getattr(self, "_evb_initial_positions", None)
        if not reactive or positions is None:
            return
        import numpy as _np

        try:
            xyz = _np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
        except AttributeError:
            xyz = _np.asarray([[v.x, v.y, v.z] for v in positions], dtype=float)
        if xyz.ndim != 2 or xyz.shape[0] == 0:
            return
        core = xyz[_np.asarray(reactive, dtype=int)]
        d_max = float(_np.sqrt(((xyz[None, :, :] - core[:, None, :]) ** 2).sum(-1)).max())
        if qm_cutoff_nm < d_max - 1e-9:
            raise ValueError(
                f"qm_cutoff_distance = {qm_cutoff_nm:.2f} nm truncates interactions involving the "
                f"reactive atoms: the system extends {d_max:.2f} nm from them. The EVB gap depends "
                f"strongly and non-monotonically on this cutoff (50-100 kcal/mol measured), because "
                f"the diabats differ in the core charges so the truncation does not cancel in "
                f"E(state2) - E(state1). Use qm_cutoff_distance >= {d_max:.2f} nm (the default 9.9 nm "
                f"matches Q's 99 A 'no cutoff on reacting fragments'), and use the SAME value on "
                f"every arm of a comparison."
            )

    def build_reactive_forces(
        self,
        *,
        reactive_morse_parameters=None,
        soft_pair_parameters_by_atom_type: Mapping[str, tuple[float, float]] | None = None,
        state_dependent_exclusion_patches: bool = True,
        **vdw_kwargs,
    ) -> None:
        """Build the EVB force stack in the canonical order, in one call.

        The correct order for a reactive system is fixed, and two of its steps are easy
        to omit because neither is part of ``auto_build_forces``::

            0. promote_reactive_constraints_to_bonds(state_system, pairs)  # see below
            1. add_vdw_forces
            2. add_bond_forces
            3. add_angle_forces
            4. add_torsion_forces
            5. replace_reactive_bonds_with_morse
            6. add_reactive_softpair_vdw_forces
            7. add_state_dependent_onethree_exclusion_patches

        This method runs steps 1-7. Step 0 cannot: it mutates each per-state
        ``openmm.System`` *before* ``EVBSystem`` exists, so call
        :func:`evbopenmm.forces.promote_reactive_constraints_to_bonds` on every state
        system first when reactive bonds may be constrained (any X-H built under
        ``constraints=HBonds`` is, which silently hides it from reactive-bond detection).

        Steps 5 and 6 are skipped when the system has no forming/breaking pairs.

        Parameters
        ----------
        reactive_morse_parameters :
            Morse parameters for the reactive bonds. Defaults to whatever was passed to
            the constructor. Required when reactive pairs exist.
        soft_pair_parameters_by_atom_type :
            Per-element Q6 soft-pair ``(C, a)`` values. Defaults to the constructor's.
        state_dependent_exclusion_patches :
            Run step 7. Leave enabled unless the reaction provably changes no pair's
            1-2 / 1-3 / 1-4 classification, in which case it is a no-op anyway.
        **vdw_kwargs :
            Forwarded to :meth:`add_vdw_forces` (e.g. ``mm_cutoff_distance``,
            ``qm_cutoff_distance``).

        Notes
        -----
        ``auto_build_forces=True`` runs steps 1-6 only, for backward compatibility: it
        has never applied step 7, and enabling it there would silently change the
        energies of existing setups. Prefer this method for new work.

        Construct with ``auto_build_forces=False`` when using this method: it adds the
        force stack unconditionally, so calling it after the constructor already built
        one would add a second copy.
        """
        morse = (reactive_morse_parameters
                 if reactive_morse_parameters is not None
                 else self._constructor_reactive_morse_parameters)
        softpair = (soft_pair_parameters_by_atom_type
                    if soft_pair_parameters_by_atom_type is not None
                    else self._constructor_softpair_parameters)

        self.add_vdw_forces(**vdw_kwargs)
        self.add_bond_forces()
        self.add_angle_forces()
        self.add_torsion_forces()

        if self.get_forming_breaking_bond_pairs():
            if morse is None:
                raise ValueError(
                    "build_reactive_forces detected forming/breaking reactive pairs but "
                    "no Morse parameters are available. Pass reactive_morse_parameters "
                    "here or to the constructor: a force-field parameterisation only "
                    "produces a harmonic bond, which cannot describe dissociation."
                )
            self.replace_reactive_bonds_with_morse(morse)
            self.add_reactive_softpair_vdw_forces(
                soft_pair_parameters_by_atom_type=softpair,
            )

        if state_dependent_exclusion_patches:
            self.add_state_dependent_onethree_exclusion_patches()

    def validate_diabat_topology(
        self,
        *,
        ignore_residue_names: Iterable[str] | None = ("HOH", "WAT", "TIP3"),
        check_angles: bool = True,
        raise_on_error: bool = False,
    ) -> Mapping[str, tuple]:
        """Audit every diabat's own System against its own bond graph.

        Runs both force-field integrity checks per state:

        * :func:`evbopenmm.forces.audit_dropped_bonded_terms` -- bonds and angles the
          force field silently skipped for want of a matching parameter.
        * :func:`evbopenmm.forces.audit_missing_exclusions` -- 1-2 / 1-3 pairs that are
          not excluded from nonbonded interactions.

        This is the check for hand-assembled diabats. ``prebuilt_state_systems`` are
        trusted verbatim -- :meth:`_validate_prebuilt_state_system` compares only
        particle count and masses -- so a state whose bond list was edited without
        regenerating angles and exclusions is accepted silently. The forming bond then
        makes its neighbours 1-3, they keep interacting at full Coulomb, and in a
        product diabat with large formal charges that is enough to push the reacting
        fragments apart while the energy gap stays finite.

        Detection only: nothing is repaired. Regenerating exclusions on a built System
        is not safely idempotent (OpenMM rejects a duplicate exception), so the fix
        belongs upstream, in whatever produced the state.

        Parameters
        ----------
        ignore_residue_names :
            Residue names to skip. Defaults to bulk water, whose rigid-model constraints
            and exclusions are expected to differ from a plain bond-graph derivation.
        check_angles :
            Passed through to the bonded-term audit.
        raise_on_error :
            Raise ``RuntimeError`` on the first state that fails, instead of warning.

        Returns
        -------
        Mapping[str, tuple]
            ``{state_name: (DroppedTermReport, MissingExclusionReport)}``.
        """
        from .forces import (
            audit_dropped_bonded_terms,
            audit_missing_exclusions,
            bond_graph_from_system,
        )

        results: Dict[str, tuple] = {}
        systems = self._state_systems
        for state_name in self.state_names:
            system = systems.get(state_name) if systems else None
            modeller = self._modellers.get(state_name)
            if system is None or modeller is None:
                continue
            topology = modeller.topology
            # The System, not the Topology, is authoritative for a diabat's
            # connectivity. A Topology read from a PDB without CONECT records carries
            # standard template bonds, so for a reactive state it still holds the bonds
            # this diabat has broken and none of the ones it has formed -- auditing
            # against it flags the broken bond and misses the formed one.
            state_bonds = bond_graph_from_system(system)
            # check_bonds=False: state_bonds came from this System, so a bond check
            # against it compares a set with itself and can never fail. Angles and
            # exclusions, derived from that graph, are the meaningful checks.
            bonded = audit_dropped_bonded_terms(
                system, topology, check_angles=check_angles, check_bonds=False,
                bonds=state_bonds, ignore_residue_names=ignore_residue_names,
                warn=False)
            exclusions = audit_missing_exclusions(
                system, topology, bonds=state_bonds,
                ignore_residue_names=ignore_residue_names, warn=False)
            results[state_name] = (bonded, exclusions)

            if bonded.ok and exclusions.ok:
                continue
            details = "\n".join(
                part for part in (
                    None if bonded.ok else bonded.describe(topology),
                    None if exclusions.ok else exclusions.describe(topology),
                ) if part
            )
            message = (
                f"Diabatic state '{state_name}' has a force field inconsistent with its "
                f"own bond graph. Hand-assembled states are trusted verbatim, so this is "
                f"not caught at construction.\n{details}"
            )
            if raise_on_error:
                raise RuntimeError(message)
            warnings.warn(message, stacklevel=2)
        return results

    def validate_state_dependent_exclusions(self, *, raise_on_error: bool = False) -> None:
        """Warn if reactive bonds exist but state-dependent exclusions were never patched.

        The scaffold's nonbonded exclusion list is built once, from the initial state's
        bond topology. Any pair whose 1-2 / 1-3 / 1-4 classification *differs* between
        diabats is therefore wrong in every state but the initial one: a forming bond
        makes its neighbours 1-3 in the product, and those pairs must be excluded there
        even though they are non-bonded in the reactant.

        :meth:`add_state_dependent_onethree_exclusion_patches` reconciles this per state
        from each diabat's own bond graph, but it is **not** part of
        ``auto_build_forces`` and must be called explicitly. Omitting it leaves the
        forming-bond neighbours interacting through full Coulomb in the state where they
        should be excluded. Observed impact: ~30 kcal/mol of PMF asymmetry for identity
        SN2 (see that method's docstring), and, where the product diabat carries large
        formal charges, a repulsion strong enough to drive the reacting fragments apart
        so the reaction is never sampled.

        Emits a :class:`UserWarning` (or raises ``RuntimeError`` when ``raise_on_error``
        is set) when forming/breaking pairs exist and the patch has not been applied.
        Like :meth:`validate_reactive_treatment` this is a heuristic: it checks that the
        step ran, not that every exclusion is correct.
        """
        if self._onethree_patch_applied:
            return
        try:
            has_reactive = bool(self.get_forming_breaking_bond_pairs())
        except Exception:  # pragma: no cover - detection must never break a run
            return
        if not has_reactive:
            return
        msg = (
            "EVB forming/breaking bonds are present but "
            "add_state_dependent_onethree_exclusion_patches() was never called, so the "
            "scaffold's nonbonded exclusions still reflect only the initial state's "
            "topology. Pairs whose 1-2/1-3/1-4 class differs between diabats are then "
            "mistreated in every other state -- a forming bond's 1-3 neighbours stay "
            "fully interacting in the product diabat instead of being excluded. Call "
            "add_state_dependent_onethree_exclusion_patches() after add_vdw_forces() / "
            "add_bond_forces() / add_angle_forces() / add_torsion_forces(). If the "
            "reaction genuinely does not change any pair's 1-2/1-3/1-4 class, ignore "
            "this warning."
        )
        if raise_on_error:
            raise RuntimeError(msg)
        warnings.warn(msg, stacklevel=2)

    def set_state_lambdas(self, simulation: app.Simulation, lambda_values: Iterable[float] | float) -> Mapping[str, float]:
        """
        Apply lambda parameters to the EVB state-specific forces enforcing a unit-sum constraint.

        Parameters
        ----------
        simulation :
            OpenMM Simulation whose context already uses this EVBSystem's forces.
        lambda_values :
            Iterable containing ``N-1`` lambda values ordered according to :attr:`state_names`.
            A single float is also accepted for the common two-state EVB case.
            The final state's lambda is computed as ``1 - sum(lambda_values)``.

        Returns
        -------
        Mapping[str, float]
            Mapping of ``state_name -> applied lambda`` for convenience.
        """

        # Fire the reactive-treatment check once, right before the first run: by the time
        # lambdas are set the force setup is complete, so this catches a Morse-without-
        # soft-pair build regardless of the order forces were added.
        if not getattr(self, "_reactive_treatment_checked", False):
            self._reactive_treatment_checked = True
            self.validate_reactive_treatment()
            self.validate_state_dependent_exclusions()

        # EVB must be integrated at 1 fs. The reactive Morse / soft-pair coordinate is far stiffer than
        # a normal bond, so a 2 fs step corrupts the barrier; hydrogen-mass repartitioning does NOT help
        # (it lightens only X-H stretches, never the heavy-atom reactive bond). Warn once, right before
        # the first dynamics step of any EVB run.
        if not getattr(self, "_timestep_checked", False):
            self._timestep_checked = True
            try:
                _dt_fs = simulation.integrator.getStepSize().value_in_unit(unit.femtoseconds)
            except Exception:  # pragma: no cover - defensive; never block a run on the check
                _dt_fs = None
            if _dt_fs is not None and _dt_fs > 1.0 + 1e-6:
                warnings.warn(
                    f"EVB is being integrated at {_dt_fs:.2f} fs, but EVB must ALWAYS run at 1 fs. The "
                    "reactive Morse/soft-pair coordinate is too stiff for a larger timestep and the "
                    "barrier will be corrupted; hydrogen-mass repartitioning does not fix this (it only "
                    "lightens X-H stretches, not the heavy-atom reactive bond). Set the integrator step "
                    "size to 1 fs (and do not use HMR).",
                    stacklevel=2,
                )

        state_order = self.state_names
        lambda_sequence = self._coerce_lambda_sequence(lambda_values)
        context = simulation.context
        applied: Dict[str, float] = {}
        for state_name, value in zip(state_order, lambda_sequence):
            parameter_name = self._state_lambda_parameter_names.get(state_name) or f"lambda_{state_name}"
            if self._context_parameter_exists(context, parameter_name):
                try:
                    context.setParameter(parameter_name, value)
                except openmm.OpenMMException as exc:  # pragma: no cover - surface clearer context
                    raise RuntimeError(
                        f"Failed to set lambda parameter '{parameter_name}' for state '{state_name}'."
                    ) from exc
            applied[state_name] = value
        self._sync_pme_g_parameters(context)
        return applied

    def get_state_lambdas(self, simulation: app.Simulation) -> Mapping[str, float]:
        """
        Retrieve the current lambda parameters for every EVB state from the simulation context.

        Parameters
        ----------
        simulation :
            OpenMM Simulation whose context uses this EVBSystem's forces.

        Returns
        -------
        Mapping[str, float]
            Mapping of ``state_name -> lambda`` reporting the context's current parameter values.
        """

        context = simulation.context
        values: Dict[str, float] = {}
        for state_name in self.state_names:
            parameter_name = self._state_lambda_parameter_names.get(state_name) or f"lambda_{state_name}"
            values[state_name] = self._context_get_parameter_or_default(context, parameter_name, default=0.0)
        return values

    def active_state_integration_force_groups(self, active_state: str) -> int:
        """Force-group bitmask that stays nonzero when ``active_state`` is the only
        weighted diabatic state (its lambda = 1, all others 0): the shared groups
        OR'd with that state's per-state groups.

        At a diabatic endpoint every other state's forces are scaled to exactly
        zero (their ``g``/``lambda`` factors vanish), so integrating with only these
        groups is *exact* -- it skips redundant work (chiefly the second per-state
        PME FFT), not an approximation. See :meth:`restrict_integration_to_endpoint`.
        """
        per_state, shared_mask = self._state_force_group_masks()
        if active_state not in per_state:
            raise ValueError(f"Unknown state '{active_state}'; expected one of {self.state_names}.")
        return shared_mask | per_state[active_state]

    def restrict_integration_to_endpoint(
        self,
        simulation: app.Simulation,
        *,
        verify: bool = True,
        tol_kj_mol: float = 1e-2,
        regular_ff: bool = False,
        tol_regular_kj_mol: float = 5.0,
    ) -> int:
        """Skip zero-weighted diabatic states' forces during integration.

        Valid ONLY when the context's lambdas are pinned to a single diabatic
        endpoint (one state at lambda=1, the rest at 0) -- e.g. plain ground-state
        MD on one diabat, or equilibrating an endpoint window. The excluded forces
        contribute exactly zero, so the dynamics are unchanged; the win is not
        evaluating them (chiefly the redundant second PME). Undo with
        :meth:`release_integration_force_groups`. Do NOT use across a mapping-FEP
        lambda schedule -- intermediate lambdas weight every state.

        Note this only affects ``integrator.step`` dynamics; ``minimizeEnergy`` and
        ``getState`` always use the full force set, so diabatic-energy sampling for
        the EVB analysis is unaffected.

        With ``verify`` (default), asserts the masked potential energy equals the
        full potential energy at the current context state within ``tol_kj_mol``,
        proving the lambdas are at an endpoint and no always-on term hides in a
        per-state group. The real failure mode (a non-lambda-scaled per-state force)
        is orders of magnitude larger than mixed-precision rounding, so the loose
        default tolerance still catches it.

        ``regular_ff`` (endpoint regular-force-field integration): integrate the
        active state's reactive bond as a **harmonic matched to the Morse minimum**
        instead of the Morse itself. The Morse's restoring force flattens toward
        dissociation, which goes unstable at a 2 fs timestep; the matched harmonic
        (same ``r0``, curvature ``k = 2*de*alpha^2``) agrees with the Morse near
        ``r0`` -- where the endpoint ensemble lives -- so the sampled distribution
        stays consistent with the EVB diabat while integrating stably. This swaps
        the integrated force groups only; ``getState``/``minimizeEnergy`` and the
        diabatic-energy sampler still use the exact Morse (the harmonic sits in a
        force group outside the energy registry and is gated by ``g_regular_ff``,
        which this method sets to 1). With ``verify`` the harmonic-only and
        Morse-only group energies must agree within ``tol_regular_kj_mol`` at the
        current coordinates (a guard that the bond is near ``r0``).

        Returns the applied force-group bitmask.
        """
        context = simulation.context
        # Reset the regular-FF gate before the endpoint-energy check so a stale
        # value from a prior call cannot perturb the full-vs-masked comparison.
        self._set_regular_ff_gate(context, 0.0)
        lambdas = self.get_state_lambdas(simulation)
        active = [s for s, v in lambdas.items() if abs(float(v) - 1.0) <= 1e-9]
        rest_zero = all(abs(float(v)) <= 1e-9 for s, v in lambdas.items() if s not in active)
        if len(active) != 1 or not rest_zero:
            raise RuntimeError(
                "restrict_integration_to_endpoint requires lambdas pinned to a single diabatic "
                f"endpoint (one state=1, the rest=0); got {dict(lambdas)}."
            )
        active_state = active[0]
        mask = self.active_state_integration_force_groups(active_state)
        if verify:
            e_full = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)
            e_mask = context.getState(getEnergy=True, groups=mask).getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)
            if abs(e_full - e_mask) > tol_kj_mol:
                raise RuntimeError(
                    f"restrict_integration_to_endpoint: masked energy {e_mask:.6f} != full "
                    f"{e_full:.6f} kJ/mol (|diff|={abs(e_full - e_mask):.3e} > {tol_kj_mol}). A "
                    f"zero-weighted state's force is not fully lambda-scaled, or the lambdas are not "
                    f"exactly at the '{active_state}' endpoint."
                )
        if regular_ff:
            if active_state not in self._endpoint_regular_ff_groups:
                raise RuntimeError(
                    f"restrict_integration_to_endpoint(regular_ff=True): the '{active_state}' state "
                    "has no endpoint regular-FF harmonic (no reactive Morse bond bonded in this "
                    "state). regular_ff is only meaningful for states with a reactive bond."
                )
            morse_group = self._force_group_ids[(active_state, "reactive_morse")]
            harmonic_group = self._endpoint_regular_ff_groups[active_state]
            mask = (mask & ~(1 << morse_group)) | (1 << harmonic_group)
            self._set_regular_ff_gate(context, 1.0)
            if verify:
                e_h = context.getState(
                    getEnergy=True, groups=(1 << harmonic_group)).getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole)
                e_m = context.getState(
                    getEnergy=True, groups=(1 << morse_group)).getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole)
                if abs(e_h - e_m) > tol_regular_kj_mol:
                    self._set_regular_ff_gate(context, 0.0)
                    raise RuntimeError(
                        f"restrict_integration_to_endpoint(regular_ff=True): matched-harmonic energy "
                        f"{e_h:.4f} != Morse energy {e_m:.4f} kJ/mol "
                        f"(|diff|={abs(e_h - e_m):.3e} > {tol_regular_kj_mol}). The reactive bond is "
                        f"far from its minimum r0 in the '{active_state}' state; the harmonic "
                        "approximation is only valid near r0."
                    )
        simulation.integrator.setIntegrationForceGroups(mask)
        return mask

    def _set_regular_ff_gate(self, context: openmm.Context, value: float) -> None:
        """Set the endpoint regular-FF global parameter if the context defines it."""
        try:
            context.getParameter(self._regular_ff_global_param)
        except Exception:
            return
        context.setParameter(self._regular_ff_global_param, float(value))

    def release_integration_force_groups(self, simulation: app.Simulation) -> None:
        """Restore evaluation of all force groups during integration (undo
        :meth:`restrict_integration_to_endpoint`), and disable the endpoint
        regular-FF harmonic gate."""
        self._set_regular_ff_gate(simulation.context, 0.0)
        simulation.integrator.setIntegrationForceGroups(0xFFFFFFFF)

    def build_ground_state_md_system(
        self,
        state: str,
        *,
        hydrogen_mass_amu: float = 0.0,
        repartition_water: bool = False,
    ) -> openmm.System:
        """Standalone solvated MM ``System`` for ONE diabatic endpoint -- the
        cheap, 2 fs-stable ground-state sampler (stage 1 of the staged EVB
        sampling protocol).

        At a physical endpoint the bond topology is fixed, so the diabat is plain
        MM: the reactive bond becomes the **harmonic matched to the Morse minimum**
        (``k = 2*de*alpha^2``, same ``r0``) -- the Morse's restoring force flattens
        toward dissociation and goes unstable at 2 fs, the harmonic does not -- while
        the forming/breaking pair keeps its bounded soft-pair. Everything else
        (protein/solvent bonded, LJ, the state's PME charges, the 1-4 overlays, the
        reactive-distance wall, the box, the constraints baked at bundle build) is
        identical to the EVB diabat; the other states' forces are pinned to
        ``lambda=0``. Because this is a freshly **cloned** System, hydrogen-mass
        repartitioning can be applied here (``hydrogen_mass_amu``) -- the lever that
        pushes the timestep to 2 fs+ by slowing the protein/ligand H-angle modes --
        without tripping the EVB loader's mass-validation guard.

        The slow environment ensemble matches the EVB diabat (the FF differs only by
        harmonic-vs-Morse near ``r0``), so a long 2 fs run here converges the slow
        DOF; seed the EVB endpoint from the resulting frames for a short Morse/
        soft-pair refinement. Valid ONLY at a diabatic endpoint -- never intermediate
        lambda; the harmonic+LJ replacement is correct only at fixed topology.

        Parameters
        ----------
        state :
            Diabatic endpoint to build (must be one of ``state_names`` and carry a
            reactive Morse bond, i.e. have an endpoint regular-FF harmonic).
        hydrogen_mass_amu :
            If > 0, repartition each (non-water by default) hydrogen to this mass,
            pulling the difference from its bonded heavy atom (reactive atoms carry
            no H, so the reactive bond is untouched). 0 disables HMR. Masses do not
            affect the configurational ensemble, so HMR is bias-free for the LRA.
        repartition_water :
            Include water hydrogens in HMR (default False -- water is already rigid,
            so it gains nothing and SETTLE handles it either way).
        """
        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("Build the EVB system before 'build_ground_state_md_system'.")
        if state not in self.state_names:
            raise ValueError(
                f"Unknown state '{state}'; expected one of {self.state_names}."
            )
        if state not in self._endpoint_regular_ff_groups:
            raise ValueError(
                f"State '{state}' has no reactive bonded pair / endpoint regular-FF "
                "harmonic. 'build_ground_state_md_system' is only valid at a diabatic "
                "endpoint with a reactive Morse bond; call "
                "'replace_reactive_bonds_with_morse' first."
            )

        # Deep clone so we can rewrite forces/masses without touching the scaffold.
        system = openmm.XmlSerializer.deserialize(
            openmm.XmlSerializer.serialize(self._evb_system_scaffold)
        )

        # Pin global-parameter defaults to the target diabat: target lambda=1 (PME
        # g=1), all other states 0, and activate the matched harmonic (g_regular_ff=1).
        targets: Dict[str, float] = {self._regular_ff_global_param: 1.0}
        for s in self.state_names:
            lam = 1.0 if s == state else 0.0
            targets[self._state_lambda_parameter_names[s]] = lam
            g_name = self._state_pme_g_param_names.get(s)
            if g_name is not None:
                targets[g_name] = math.sqrt(lam)
        for fi in range(system.getNumForces()):
            force = system.getForce(fi)
            if not hasattr(force, "getNumGlobalParameters"):
                continue
            for gi in range(force.getNumGlobalParameters()):
                name = force.getGlobalParameterName(gi)
                if name in targets:
                    force.setGlobalParameterDefaultValue(gi, float(targets[name]))

        # Zero the Morse bonds: the matched harmonic (gated on, lambda=1) now carries
        # the bonded reactive pair, so the Morse must not double-count. (The inactive
        # state's Morse is already off via lambda=0; zeroing it too is harmless.)
        for fi in range(system.getNumForces()):
            force = system.getForce(fi)
            if not isinstance(force, openmm.CustomBondForce):
                continue
            per_bond = [
                force.getPerBondParameterName(i)
                for i in range(force.getNumPerBondParameters())
            ]
            if per_bond == ["r0", "de", "alpha"] and "exp(-alpha" in (
                force.getEnergyFunction().replace(" ", "")
            ):
                for bi in range(force.getNumBonds()):
                    a1, a2, params = force.getBondParameters(bi)
                    force.setBondParameters(bi, a1, a2, [params[0], 0.0, params[2]])

        if hydrogen_mass_amu and hydrogen_mass_amu > 0.0:
            self._apply_hydrogen_mass_repartitioning(
                system, float(hydrogen_mass_amu), repartition_water=repartition_water
            )
        return system

    def apply_hydrogen_mass_repartitioning(
        self,
        hydrogen_mass_amu: float = 3.0,
        *,
        system: Optional[openmm.System] = None,
        repartition_water: bool = False,
    ) -> int:
        """Activate hydrogen-mass repartitioning (HMR) on an EVB ``System`` **in place**,
        so a 2 fs+ timestep is stable -- the standard procedure for the aqueous and enzyme
        FEP/ground-state samplers alike.

        HMR slows the fast X-H stretch/angle modes (which otherwise cap the timestep) by
        moving ``hydrogen_mass_amu - m_H`` Da from each heavy atom onto its bonded hydrogen;
        total mass per H-heavy pair is conserved. It changes only the *dynamics*, not the
        configurational ensemble, so free energies (dG!/dG0) are unchanged -- it is a pure
        speed lever. The reactive Morse coordinate is heavy-atom (reactive atoms carry no H),
        so it is never touched.

        Parameters
        ----------
        hydrogen_mass_amu :
            Target hydrogen mass in Da (typical 3.0-4.0 for a 2 fs step). Must exceed the
            physical H mass (~1.008).
        system :
            System to modify. Defaults to the built EVB scaffold (``_evb_system_scaffold``),
            i.e. the object the FEP simulation runs on. Pass an explicit System to repartition
            a clone instead.
        repartition_water :
            Include water hydrogens (default ``False`` -- water is rigid/SETTLE, gains nothing).

        Returns
        -------
        int
            Number of hydrogens repartitioned.
        """
        if hydrogen_mass_amu is None or float(hydrogen_mass_amu) <= 1.008:
            raise ValueError(
                f"hydrogen_mass_amu must exceed the physical H mass (~1.008 Da); "
                f"got {hydrogen_mass_amu!r} (typical 3.0-4.0)."
            )
        if system is None:
            if self._evb_system_scaffold is None:
                raise ValueError(
                    "Build the EVB system before applying HMR (no scaffold present)."
                )
            system = self._evb_system_scaffold
        if self._evb_system_source_state is None or self._evb_system_source_state not in self._modellers:
            raise ValueError(
                "HMR needs the source-state topology; build/load the EVB system first."
            )
        return self._apply_hydrogen_mass_repartitioning(
            system, float(hydrogen_mass_amu), repartition_water=repartition_water
        )

    def _apply_hydrogen_mass_repartitioning(
        self, system: openmm.System, hydrogen_mass_amu: float, *, repartition_water: bool
    ) -> int:
        """Set every (non-water) hydrogen to ``hydrogen_mass_amu`` Da, pulling the
        added mass from its bonded heavy atom (total mass per H-heavy pair conserved).
        Reactive atoms carry no H and so are untouched. Returns the count repartitioned."""
        topology = self._modellers[self._evb_system_source_state].topology
        water_residues = {"HOH", "WAT", "H2O", "TIP3", "SOL"}
        n_repartitioned = 0
        for bond in topology.bonds():
            a1, a2 = bond.atom1, bond.atom2
            e1 = getattr(a1.element, "symbol", None)
            e2 = getattr(a2.element, "symbol", None)
            if e1 == "H" and e2 != "H":
                h_atom, heavy_atom = a1, a2
            elif e2 == "H" and e1 != "H":
                h_atom, heavy_atom = a2, a1
            else:
                continue
            if (not repartition_water) and heavy_atom.residue.name in water_residues:
                continue
            cur_h = system.getParticleMass(h_atom.index).value_in_unit(unit.dalton)
            cur_heavy = system.getParticleMass(heavy_atom.index).value_in_unit(unit.dalton)
            delta = hydrogen_mass_amu - cur_h
            if delta <= 0.0 or cur_heavy - delta <= 0.0:
                continue
            system.setParticleMass(h_atom.index, hydrogen_mass_amu * unit.dalton)
            system.setParticleMass(heavy_atom.index, (cur_heavy - delta) * unit.dalton)
            n_repartitioned += 1
        return n_repartitioned

    _coerce_length_nm = staticmethod(_units.coerce_length_nm)
    _coerce_molarity = staticmethod(_units.coerce_molarity)
    _coerce_box_size = staticmethod(_units.coerce_box_size)
    _coerce_temperature = staticmethod(_units.coerce_temperature)
    _coerce_timestep = staticmethod(_units.coerce_timestep)
    _coerce_friction = staticmethod(_units.coerce_friction)
    _coerce_force_constant_kj_mol_nm2 = staticmethod(_units.coerce_force_constant_kj_mol_nm2)

    @staticmethod
    def _is_hydrogen_like_atom(atom: app.Topology.Atom) -> bool:
        element = atom.element
        if element is not None and int(element.atomic_number) == 1:
            return True
        name = (atom.name or "").strip().upper()
        return bool(name) and name.startswith("H")

    @staticmethod
    def _default_two_state_lambda_schedule(
        *,
        lambda_start: float,
        lambda_end: float,
        lambda_step: float,
    ) -> Tuple[float, ...]:
        start = float(lambda_start)
        end = float(lambda_end)
        step = float(lambda_step)
        if abs(step) <= 1e-12:
            raise ValueError("'lambda_step' must be non-zero.")
        if step > 0.0 and start > end:
            raise ValueError("Positive 'lambda_step' requires 'lambda_start <= lambda_end'.")
        if step < 0.0 and start < end:
            raise ValueError("Negative 'lambda_step' requires 'lambda_start >= lambda_end'.")

        tolerance = 1e-8
        values: list[float] = []
        current = start
        if step > 0.0:
            while current <= end + tolerance:
                values.append(current)
                current += step
        else:
            while current >= end - tolerance:
                values.append(current)
                current += step

        if not values:
            values = [start]
        if abs(values[-1] - end) > tolerance:
            values.append(end)

        clamped: list[float] = []
        for value in values:
            if value < -tolerance or value > 1.0 + tolerance:
                raise ValueError(
                    f"Two-state schedule generated lambda={value}, outside [0, 1]. "
                    "Adjust lambda_start/lambda_end/lambda_step."
                )
            clipped = min(max(value, 0.0), 1.0)
            clamped.append(clipped)
        return tuple(clamped)

    @staticmethod
    def _default_linear_lambda_schedule(
        *,
        n_states: int,
        lambda_start: float,
        lambda_end: float,
        lambda_step: float,
    ) -> Tuple[Tuple[float, ...], ...]:
        if int(n_states) < 2:
            raise ValueError("Linear lambda schedule requires at least two states.")
        scalar_schedule = EVBSystem._default_two_state_lambda_schedule(
            lambda_start=lambda_start,
            lambda_end=lambda_end,
            lambda_step=lambda_step,
        )
        zero_fillers = (0.0,) * (int(n_states) - 2)
        windows: list[Tuple[float, ...]] = []
        for value in scalar_schedule:
            residual = 1.0 - float(value)
            if residual < 0.0:
                residual = 0.0
            windows.append((float(value),) + zero_fillers + (residual,))
        return tuple(windows)

    def _normalize_full_lambda_sequence(self, values: Sequence[float]) -> Tuple[float, ...]:
        n_states = len(self.state_names)
        if len(values) != n_states:
            raise ValueError(f"Expected {n_states} full-state lambdas but received {len(values)}.")
        tolerance = 1e-8
        normalized = [float(value) for value in values]
        for idx, value in enumerate(normalized):
            if value < -tolerance:
                raise ValueError(f"Lambda for state '{self.state_names[idx]}' is negative ({value}).")
            if value < 0.0:
                normalized[idx] = 0.0
        total = sum(normalized)
        if abs(total - 1.0) > tolerance:
            raise ValueError(
                f"Full-state lambda sequence must sum to 1 (received {total})."
            )
        if abs(total) <= tolerance:
            raise ValueError("Full-state lambda sequence sums to zero.")
        if abs(total - 1.0) > 1e-12:
            normalized = [value / total for value in normalized]
        return tuple(normalized)

    def _resolve_lambda_windows(
        self,
        *,
        lambda_schedule: Sequence[float] | Sequence[Sequence[float]] | Sequence[Mapping[str, float]] | None,
        lambda_start: float,
        lambda_end: float,
        lambda_step: float,
    ) -> Tuple[Tuple[float, ...], ...]:
        n_states = len(self.state_names)
        windows: list[Tuple[float, ...]] = []

        if lambda_schedule is None:
            if n_states == 2:
                for lambda_value in self._default_two_state_lambda_schedule(
                    lambda_start=lambda_start,
                    lambda_end=lambda_end,
                    lambda_step=lambda_step,
                ):
                    windows.append(self._coerce_lambda_sequence(lambda_value))
            else:
                for lambda_vector in self._default_linear_lambda_schedule(
                    n_states=n_states,
                    lambda_start=lambda_start,
                    lambda_end=lambda_end,
                    lambda_step=lambda_step,
                ):
                    windows.append(self._normalize_full_lambda_sequence(lambda_vector))
        else:
            for entry in lambda_schedule:
                if isinstance(entry, Mapping):
                    if all(state in entry for state in self.state_names):
                        full = tuple(float(entry[state]) for state in self.state_names)
                        windows.append(self._normalize_full_lambda_sequence(full))
                    else:
                        missing = [state for state in self.state_names[:-1] if state not in entry]
                        if missing:
                            joined = ", ".join(missing)
                            raise ValueError(
                                "Lambda mapping entries must include either all states or at least "
                                f"states {self.state_names[:-1]}; missing: {joined}."
                            )
                        provided = tuple(float(entry[state]) for state in self.state_names[:-1])
                        windows.append(self._coerce_lambda_sequence(provided))
                elif isinstance(entry, Real):
                    windows.append(self._coerce_lambda_sequence(float(entry)))
                else:
                    sequence = tuple(float(value) for value in entry)
                    if len(sequence) == n_states:
                        windows.append(self._normalize_full_lambda_sequence(sequence))
                    else:
                        windows.append(self._coerce_lambda_sequence(sequence))

        if len(windows) < 2:
            raise ValueError("At least two lambda windows are required for FEP.")
        return tuple(windows)

    def _resolve_reaction_coordinate_coefficients(
        self,
        coefficients: Sequence[float] | None,
    ) -> Tuple[float, ...]:
        n_states = len(self.state_names)
        if coefficients is None:
            defaults = [0.0] * n_states
            defaults[0] = 1.0
            if n_states > 1:
                defaults[1] = -1.0
            return tuple(defaults)

        resolved = tuple(float(value) for value in coefficients)
        if len(resolved) != n_states:
            raise ValueError(
                f"Expected {n_states} reaction-coordinate coefficients, received {len(resolved)}."
            )
        if all(abs(value) <= 1e-14 for value in resolved):
            raise ValueError("Reaction-coordinate coefficients cannot all be zero.")
        return resolved

    def _sample_diabatic_state_energies_kj_mol(self, simulation: app.Simulation) -> Tuple[float, ...]:
        state_masks, shared_mask = self._state_force_group_masks()
        return sample_diabatic_state_energies_via_force_groups_kj_mol(
            simulation=simulation,
            state_names=self.state_names,
            state_group_masks=state_masks,
            shared_group_mask=shared_mask,
            context_lambda_snapshot=self._context_lambda_snapshot,
            apply_context_lambdas=self._apply_context_lambdas,
        )

    def _state_force_group_masks(self) -> Tuple[Dict[str, int], int]:
        """Build (per_state_mask, shared_mask) from the force-group registry.

        ``shared_mask`` covers groups not associated with any specific state
        (shared MM/QM nonbonded, bonded reference forces, SCAAS forces, etc.).
        ``per_state_mask[state]`` covers the per-state force groups for that
        diabatic state. Together, ``shared_mask | per_state_mask[state]``
        names exactly the forces that contribute to the state's diabatic
        energy when its ``lambda_<state>`` is set to 1.0.
        """
        shared_mask = 0
        per_state: Dict[str, int] = {state: 0 for state in self.state_names}
        for group_id, (state, _term) in self._force_group_definitions.items():
            bit = 1 << int(group_id)
            if state is None:
                shared_mask |= bit
            elif state in per_state:
                per_state[state] |= bit
        return per_state, shared_mask

    @staticmethod
    def _effective_energy_from_state_energies(
        lambda_values: Sequence[float],
        state_energies_kj_mol: Sequence[float],
    ) -> float:
        return effective_energy_from_state_energies(lambda_values, state_energies_kj_mol)

    @staticmethod
    def _reaction_coordinate_from_state_energies(
        coefficients: Sequence[float],
        state_energies_kj_mol: Sequence[float],
    ) -> float:
        return reaction_coordinate_from_state_energies(coefficients, state_energies_kj_mol)

    @staticmethod
    def _energy_scale_to_kj(input_energy_unit: str) -> float:
        return energy_scale_to_kj(input_energy_unit)

    @staticmethod
    def _simulation_timestep_ps(
        simulation: app.Simulation,
        *,
        fallback_timestep: unit.Quantity | float | None = None,
    ) -> float:
        return simulation_timestep_ps(
            simulation,
            fallback_timestep=fallback_timestep,
            coerce_timestep=_units.coerce_timestep,
        )

    @staticmethod
    def _steps_from_report_interval_ps(interval_ps: float, timestep_ps: float) -> int:
        return steps_from_report_interval_ps(interval_ps, timestep_ps)

    def _create_langevin_middle_simulation(
        self,
        *,
        topology: app.Topology,
        system: openmm.System,
        temperature: unit.Quantity,
        friction: unit.Quantity,
        timestep: unit.Quantity,
        random_seed: int | None,
        platform: openmm.Platform | None,
        platform_properties: Mapping[str, str] | None,
    ) -> app.Simulation:
        integrator = openmm.LangevinMiddleIntegrator(temperature, friction, timestep)
        if random_seed is not None:
            integrator.setRandomNumberSeed(int(random_seed))
        simulation: app.Simulation
        if platform is None:
            simulation = app.Simulation(topology, system, integrator)
        else:
            simulation = app.Simulation(topology, system, integrator, platform, dict(platform_properties or {}))
        self._attach_scaas_to_simulation(simulation, initialize=False)
        return simulation

    @staticmethod
    def _collect_md_energy_samples(
        simulation: app.Simulation,
        *,
        equilibration_steps: int,
        production_steps: int,
        sample_interval: int,
    ) -> Tuple[MDEnergySampleRecord, ...]:
        return collect_md_energy_samples(
            simulation,
            equilibration_steps=equilibration_steps,
            production_steps=production_steps,
            sample_interval=sample_interval,
        )

    @staticmethod
    def _summarize_md_energy_samples(samples: Sequence[MDEnergySampleRecord]) -> Dict[str, float]:
        return summarize_md_energy_samples(samples)

    @staticmethod
    def _default_equilibration_schedule() -> Tuple[Tuple[float, ...], Tuple[int, ...], Tuple[float, ...]]:
        return default_equilibration_schedule()

    def _coerce_lambda_stage_sequence(
        self,
        lambda_values: Iterable[float] | float | Mapping[str, float],
    ) -> Tuple[float, ...]:
        return coerce_lambda_stage_sequence(
            lambda_values,
            state_names=self.state_names,
            coerce_lambda_sequence=self._coerce_lambda_sequence,
            normalize_full_lambda_sequence=self._normalize_full_lambda_sequence,
        )

    def _apply_lambda_stage_sequence(
        self,
        simulation: app.Simulation,
        lambda_sequence: Sequence[float],
    ) -> Mapping[str, float]:
        return apply_lambda_stage_sequence(
            simulation,
            lambda_sequence,
            state_names=self.state_names,
            apply_context_lambdas=self._apply_context_lambdas,
        )

    @staticmethod
    def _estimate_kinetic_dof(system: openmm.System) -> int:
        return estimate_kinetic_dof(system)

    def equilibrate_schedule(
        self,
        *,
        simulation: app.Simulation | None = None,
        timesteps: Sequence[unit.Quantity | float] | None = None,
        steps: Sequence[int] | None = None,
        temperatures: Sequence[unit.Quantity | float] | None = None,
        equilibration_lambda: Iterable[float] | float | Mapping[str, float] = 0.5,
        lambda_schedule: Sequence[Iterable[float] | float | Mapping[str, float]] | None = None,
        friction: Sequence[unit.Quantity | float] | unit.Quantity | float = 1.0,
        initialize_velocities: str = "first",
        random_seed: int | None = None,
        platform: openmm.Platform | None = None,
        platform_properties: Mapping[str, str] | None = None,
        coordinate_report_interval_ps: float | None = None,
        coordinate_report_file: PathLike | None = None,
        coordinate_report_append: bool = False,
        coordinate_topology_file: PathLike | None = None,
        summary_file: PathLike | None = None,
        coordinate_csv_file: PathLike | None = None,
        coordinate_csv_interval_steps: int | None = None,
        energy_abs_cap_kj_mol: float = 1.0e9,
        monitor_interval_steps: int = 1,
        max_abs_coord_nm_cap: float | None = None,
        solute_restraints: bool = False,
        solute_restraint_atom_indices: Sequence[int] | None = None,
        solute_restraint_k_start: unit.Quantity | float = 2000.0,
        solute_restraint_k_end: unit.Quantity | float = 0.0,
        solute_restraint_k_schedule: Sequence[unit.Quantity | float] | None = None,
        solute_restraint_reference_positions: (
            unit.Quantity | Sequence[openmm.Vec3] | Sequence[Sequence[float]] | None
        ) = None,
        solute_restraint_include_hydrogens: bool = False,
        solute_restraint_solvent_residue_names: Sequence[str] | None = None,
        fail_on_instability: bool = True,
    ) -> Mapping[str, object]:
        """
        Equilibrate an EVB simulation using a staged schedule of dt/steps/temperature.

        Parameters
        ----------
        timesteps, steps, temperatures :
            Per-stage schedules. When any of them are omitted, the missing schedule
            uses the default notebook-inspired protocol:
            ``100 steps @ 0.5 fs, 25 K`` followed by ``100 steps`` windows at
            ``1.0 fs`` with target temperatures ``50..300 K`` in 25 K increments.
        equilibration_lambda :
            Lambda definition used for all stages when ``lambda_schedule`` is not provided.
            Accepts the same forms as mapping windows: scalar (two-state), N-1 partial
            sequence, full N-state sequence, or a mapping keyed by state names.
        lambda_schedule :
            Optional per-stage lambda override (same accepted value forms as
            ``equilibration_lambda``). Must have one entry per stage when provided.
        solute_restraints :
            When ``True``, apply harmonic positional restraints to selected solute atoms
            during equilibration and release them by stage according to
            ``solute_restraint_k_schedule`` (or a linear ``k_start -> k_end`` schedule).
        """

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("EVB system not initialized; call 'build_evb_system' and add forces first.")

        default_timesteps, default_steps, default_temperatures = self._default_equilibration_schedule()
        raw_timesteps = tuple(default_timesteps if timesteps is None else timesteps)
        raw_steps = tuple(default_steps if steps is None else steps)
        raw_temperatures = tuple(default_temperatures if temperatures is None else temperatures)

        n_stages = len(raw_steps)
        if n_stages <= 0:
            raise ValueError("Equilibration schedule must contain at least one stage.")
        if len(raw_timesteps) != n_stages or len(raw_temperatures) != n_stages:
            raise ValueError(
                "Schedule lengths mismatch: 'timesteps', 'steps', and 'temperatures' must all have the same length."
            )

        step_counts = tuple(int(value) for value in raw_steps)
        if any(value <= 0 for value in step_counts):
            raise ValueError("All stage step counts must be > 0.")
        timestep_schedule = tuple(_units.coerce_timestep(value) for value in raw_timesteps)
        temperature_schedule = tuple(_units.coerce_temperature(value) for value in raw_temperatures)

        if isinstance(friction, Real) or isinstance(friction, unit.Quantity):
            stage_frictions = tuple(_units.coerce_friction(friction) for _ in range(n_stages))
        else:
            raw_friction = tuple(friction)
            if len(raw_friction) != n_stages:
                raise ValueError(
                    "When 'friction' is a sequence, it must contain one value per stage "
                    f"(expected {n_stages}, received {len(raw_friction)})."
                )
            stage_frictions = tuple(_units.coerce_friction(value) for value in raw_friction)

        base_lambda = self._coerce_lambda_stage_sequence(equilibration_lambda)
        if lambda_schedule is None:
            stage_lambda_values = tuple(base_lambda for _ in range(n_stages))
        else:
            raw_lambda_schedule = tuple(lambda_schedule)
            if len(raw_lambda_schedule) != n_stages:
                raise ValueError(
                    f"'lambda_schedule' must contain {n_stages} stage entries, received {len(raw_lambda_schedule)}."
                )
            stage_lambda_values = tuple(self._coerce_lambda_stage_sequence(value) for value in raw_lambda_schedule)

        velocity_mode = validate_velocity_mode(initialize_velocities)
        return run_equilibration_schedule(
            simulation=simulation,
            evb_topology=self._modellers[self._evb_system_source_state].topology,
            evb_system_scaffold=self._evb_system_scaffold,
            evb_initial_positions=self.evb_initial_positions,
            state_names=self.state_names,
            step_counts=step_counts,
            timestep_schedule=timestep_schedule,
            temperature_schedule=temperature_schedule,
            stage_frictions=stage_frictions,
            stage_lambda_values=stage_lambda_values,
            velocity_mode=velocity_mode,
            random_seed=random_seed,
            platform=platform,
            platform_properties=platform_properties,
            coordinate_report_interval_ps=coordinate_report_interval_ps,
            coordinate_report_file=coordinate_report_file,
            coordinate_report_append=bool(coordinate_report_append),
            coordinate_topology_file=coordinate_topology_file,
            summary_file=summary_file,
            coordinate_csv_file=coordinate_csv_file,
            coordinate_csv_interval_steps=coordinate_csv_interval_steps,
            energy_abs_cap_kj_mol=float(energy_abs_cap_kj_mol),
            max_abs_coord_nm_cap=None if max_abs_coord_nm_cap is None else float(max_abs_coord_nm_cap),
            monitor_interval_steps=int(monitor_interval_steps),
            solute_restraints=bool(solute_restraints),
            solute_restraint_atom_indices=solute_restraint_atom_indices,
            solute_restraint_k_schedule=solute_restraint_k_schedule,
            solute_restraint_k_start=solute_restraint_k_start,
            solute_restraint_k_end=solute_restraint_k_end,
            solute_restraint_reference_positions=solute_restraint_reference_positions,
            solute_restraint_include_hydrogens=bool(solute_restraint_include_hydrogens),
            solute_restraint_solvent_residue_names=solute_restraint_solvent_residue_names,
            solute_restraint_global_parameter=SOLUTE_RESTRAINT_GLOBAL_PARAMETER,
            fail_on_instability=bool(fail_on_instability),
            create_langevin_middle_simulation=self._create_langevin_middle_simulation,
            attach_scaas_to_simulation=self._attach_scaas_to_simulation,
            context_parameter_exists=self._context_parameter_exists,
            simulation_timestep_ps=lambda simulation: self._simulation_timestep_ps(
                simulation,
                fallback_timestep=timestep_schedule[0],
            ),
            steps_from_report_interval_ps=self._steps_from_report_interval_ps,
            estimate_kinetic_dof=self._estimate_kinetic_dof,
            apply_lambda_stage_sequence=self._apply_lambda_stage_sequence,
            normalize_atom_index_sequence=self._normalize_atom_index_sequence,
            default_solute_restraint_atom_indices=self._default_solute_restraint_atom_indices,
            coerce_vec3_tuple=self._coerce_vec3_tuple,
            coerce_force_constant_kj_mol_nm2=lambda value, label: _units.coerce_force_constant_kj_mol_nm2(
                value,
                label=label,
            ),
            register_force_to_group=self._register_force_to_group,
        )

    def run_adiabatic_md(
        self,
        *,
        coupling_amplitude_kj_mol: float,
        coupling_width_per_kj2_mol2: float = 0.0,
        diabatic_offset_kj_mol: float = 0.0,
        integrator: str = "langevin_middle",
        temperature_kelvin: float = 300.0,
        friction_per_ps: float = 1.0,
        timestep_fs: float = 0.5,
        production_steps: int,
        equilibration_steps: int = 0,
        report_interval_steps: int = 100,
        coord_interval_steps: int | None = None,
        minimise_first: bool = True,
        minimise_max_iterations: int = 200,
        minimisation_state: str | None = None,
        random_seed: int | None = None,
        platform: str | None = None,
        output_directory: PathLike | None = None,
        umbrella_force_constant_per_kj_mol: float = 0.0,
        umbrella_dv_center_kj_mol: float = 0.0,
    ) -> Mapping[str, object]:
        """Propagate MD directly on the EVB adiabatic ground-state surface.

        The mixing parameters ``coupling_amplitude_kj_mol`` (``A`` /
        ``H_{12}``), ``coupling_width_per_kj2_mol2`` (``μ``) and
        ``diabatic_offset_kj_mol`` (``α = Δ_2 − Δ_1``) are inputs from a
        prior FEP fit — this method is the second stage of the workflow
        (FEP-fit → adiabatic-direct-sampling).

        Parameters
        ----------
        coupling_amplitude_kj_mol
            ``A`` in :math:`H_{12}(\\Delta V) = A\\,e^{-\\mu \\Delta V^2}`.
            For constant coupling pass ``coupling_width_per_kj2_mol2=0`` and
            ``A`` = ``H_{12}``.
        coupling_width_per_kj2_mol2
            ``μ`` in (kJ/mol)⁻². ``0`` (default) recovers constant coupling.
        diabatic_offset_kj_mol
            ``α``, folded onto state 2 only.
        integrator
            ``"langevin_middle"`` (default, single force eval/step) or
            ``"nve"`` (velocity-Verlet, two force evals/step).
        temperature_kelvin, friction_per_ps
            Used by Langevin only.
        timestep_fs
            Default **0.5 fs** for adiabatic mode; larger steps validated
            empirically in A4.
        production_steps, equilibration_steps
            Step counts. Equilibration runs first if non-zero.
        report_interval_steps
            CSV-trace stride (one row per ``report_interval_steps``).
        coord_interval_steps
            Optional DCD output stride. ``None`` skips coordinate output.
        minimise_first
            Energy-minimise the EVB scaffold before equilibration. The
            minimisation runs at indicator-mode lambdas for one diabatic
            state (see ``minimisation_state``) — minimising under
            ``λ_1=λ_2=1`` finds an unphysical ``V_1+V_2-V_shared`` minimum
            and is never the right choice.
        minimisation_state
            Name of the state whose indicator lambdas are used during the
            initial minimisation. Defaults to the EVB source state.
        random_seed
            Seed for the Langevin integrator's gaussian random numbers and
            for ``setVelocitiesToTemperature``. Ignored for NVE.
        platform
            OpenMM platform name (``"CUDA"``, ``"OpenCL"``, ``"CPU"``,
            ``"Reference"``). ``None`` lets OpenMM pick.
        output_directory
            If provided, writes ``adiabatic_trace.csv``, ``summary.json``,
            and (when ``coord_interval_steps`` is set) a DCD coordinate
            file to that directory. ``None`` runs in memory only.
        umbrella_force_constant_per_kj_mol
            Optional harmonic umbrella restraint force constant on ``ΔV``,
            in (kJ/mol)⁻¹. Adds ``U_umbrella = 0.5 * K * (ΔV − ΔV_center)²``
            to the adiabatic ground-state surface and the corresponding
            shift to the per-state mixing weights (Hellmann–Feynman). Use
            this for the ΔV-centered umbrella sampling that backs A6's
            ΔG profile reconstruction. Defaults to 0 (no restraint).
        umbrella_dv_center_kj_mol
            Center of the umbrella restraint, in kJ/mol. Ignored when
            ``umbrella_force_constant_per_kj_mol = 0``.

        Returns
        -------
        dict
            Keys: ``trace`` (list of dicts, one per reporter row),
            ``summary`` (means and stds of E_minus / c1sq / dV / etc.),
            ``config`` (echo of input knobs), ``files`` (output file paths
            or ``None``).
        """

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("EVB system not initialized; call build_evb_system and add forces first.")
        if int(production_steps) <= 0:
            raise ValueError("production_steps must be > 0.")
        if int(equilibration_steps) < 0:
            raise ValueError("equilibration_steps must be >= 0.")
        if int(report_interval_steps) <= 0:
            raise ValueError("report_interval_steps must be > 0.")
        if coord_interval_steps is not None and int(coord_interval_steps) <= 0:
            raise ValueError("coord_interval_steps must be > 0 or None.")

        integrator_kind = str(integrator).lower()
        if integrator_kind not in {"langevin_middle", "nve"}:
            raise ValueError(
                f"Unknown integrator {integrator!r}; expected 'langevin_middle' or 'nve'."
            )

        resolved_minimisation_state = (
            minimisation_state
            if minimisation_state is not None
            else self._evb_system_source_state
        )
        if resolved_minimisation_state not in self.state_names:
            raise ValueError(
                f"Unknown minimisation_state {resolved_minimisation_state!r}; "
                f"expected one of {self.state_names!r}."
            )

        masks = self._adiabatic_force_group_masks()
        timestep_ps = float(timestep_fs) * 1e-3

        if integrator_kind == "nve":
            cust_integrator = build_nve_adiabatic_integrator(
                timestep=float(timestep_fs) * unit.femtosecond,
                force_group_masks=masks,
                coupling_amplitude_kj_mol=coupling_amplitude_kj_mol,
                coupling_width_per_kj2_mol2=coupling_width_per_kj2_mol2,
                diabatic_offset_kj_mol=diabatic_offset_kj_mol,
                umbrella_force_constant_per_kj_mol=float(umbrella_force_constant_per_kj_mol),
                umbrella_dv_center_kj_mol=float(umbrella_dv_center_kj_mol),
            )
        else:
            cust_integrator = build_langevin_middle_adiabatic_integrator(
                timestep=float(timestep_fs) * unit.femtosecond,
                temperature=float(temperature_kelvin) * unit.kelvin,
                friction=float(friction_per_ps) / unit.picosecond,
                force_group_masks=masks,
                coupling_amplitude_kj_mol=coupling_amplitude_kj_mol,
                coupling_width_per_kj2_mol2=coupling_width_per_kj2_mol2,
                diabatic_offset_kj_mol=diabatic_offset_kj_mol,
                umbrella_force_constant_per_kj_mol=float(umbrella_force_constant_per_kj_mol),
                umbrella_dv_center_kj_mol=float(umbrella_dv_center_kj_mol),
            )

        platform_obj: openmm.Platform | None
        if platform is None:
            platform_obj = None
        else:
            platform_obj = openmm.Platform.getPlatformByName(str(platform))

        topology = self._modellers[self._evb_system_source_state].topology
        if platform_obj is None:
            simulation = app.Simulation(topology, self._evb_system_scaffold, cust_integrator)
        else:
            simulation = app.Simulation(
                topology, self._evb_system_scaffold, cust_integrator, platform_obj
            )
        simulation.context.setPositions(self.evb_initial_positions)
        # Attach SCAAS runtime updates if `enable_scaas` was called before
        # this run; no-op otherwise. SCAAS forces (sphere restraint,
        # polarisation, shell fix) are registered as state=None and so sit
        # in the universal mask — the adiabatic mixing leaves them at
        # weight 1, which is the correct behaviour for state-independent
        # geometric restraints.
        self._attach_scaas_to_simulation(simulation, initialize=True)

        # Order matters: minimisation must run on a *single* diabatic surface
        # because at λ_1=λ_2=1 (the regime adiabatic mode evaluates in)
        # OpenMM minimises the unphysical sum V_1 + V_2 - V_shared, which
        # parks the geometry in a transition-state-like configuration where
        # both diabatic states' bonded force fields are simultaneously
        # active. Subsequent dynamics from that minimum are pathological —
        # one state immediately wins, releasing tens of kJ/mol of fictitious
        # strain energy as KE and breaking energy conservation.
        if minimise_first:
            self._apply_indicator_lambdas(
                simulation.context, active_state=resolved_minimisation_state
            )
            simulation.minimizeEnergy(maxIterations=int(minimise_max_iterations))

        # Now switch to λ_1 = λ_2 = 1: the bare state-specific terms each
        # contribute their per-group force, and the integrator's w_i weights
        # mix them adiabatically.
        full_indicator = {state: 1.0 for state in self.state_names}
        self._apply_context_lambdas(simulation.context, full_indicator)

        if integrator_kind == "langevin_middle":
            if random_seed is None:
                simulation.context.setVelocitiesToTemperature(
                    float(temperature_kelvin) * unit.kelvin
                )
            else:
                simulation.context.setVelocitiesToTemperature(
                    float(temperature_kelvin) * unit.kelvin, int(random_seed)
                )

        files: Dict[str, str | None] = {
            "trace_csv": None,
            "summary_json": None,
            "coordinates_dcd": None,
        }
        output_dir_path: Path | None
        if output_directory is None:
            output_dir_path = None
            trace_path = None
        else:
            output_dir_path = Path(output_directory)
            output_dir_path.mkdir(parents=True, exist_ok=True)
            trace_path = output_dir_path / "adiabatic_trace.csv"
            files["trace_csv"] = str(trace_path)

        # Even without an output directory we capture the trace in memory by
        # writing to a temp file, parsing it back at the end. This keeps a
        # single code path for both branches.
        if trace_path is None:
            import tempfile

            tmp_handle = tempfile.NamedTemporaryFile(
                mode="w", suffix=".csv", delete=False, encoding="utf-8"
            )
            tmp_handle.close()
            trace_csv_path = Path(tmp_handle.name)
            cleanup_trace = True
        else:
            trace_csv_path = trace_path
            cleanup_trace = False

        if equilibration_steps > 0:
            simulation.step(int(equilibration_steps))

        reporter = AdiabaticReporter(
            trace_csv_path,
            integrator=cust_integrator,
            report_interval_steps=int(report_interval_steps),
            timestep_ps=timestep_ps,
        )
        simulation.reporters.append(reporter)

        coord_reporter: app.DCDReporter | None = None
        if coord_interval_steps is not None and output_dir_path is not None:
            coord_path = output_dir_path / "coordinates.dcd"
            coord_reporter = app.DCDReporter(
                str(coord_path), reportInterval=int(coord_interval_steps)
            )
            simulation.reporters.append(coord_reporter)
            # Drop a topology PDB next to it for downstream tooling.
            topology_pdb = output_dir_path / "coordinates.pdb"
            start_state = simulation.context.getState(getPositions=True)
            with topology_pdb.open("w", encoding="utf-8") as handle:
                app.PDBFile.writeFile(simulation.topology, start_state.getPositions(), handle)
            files["coordinates_dcd"] = str(coord_path)

        try:
            simulation.step(int(production_steps))
        finally:
            reporter.close()
            simulation.reporters.remove(reporter)
            if coord_reporter is not None:
                simulation.reporters.remove(coord_reporter)

        trace_rows = _read_adiabatic_trace_csv(trace_csv_path)
        summary = _summarize_adiabatic_trace(trace_rows)
        config = {
            "coupling_amplitude_kj_mol": float(coupling_amplitude_kj_mol),
            "coupling_width_per_kj2_mol2": float(coupling_width_per_kj2_mol2),
            "diabatic_offset_kj_mol": float(diabatic_offset_kj_mol),
            "integrator": integrator_kind,
            "temperature_kelvin": float(temperature_kelvin),
            "friction_per_ps": float(friction_per_ps),
            "timestep_fs": float(timestep_fs),
            "production_steps": int(production_steps),
            "equilibration_steps": int(equilibration_steps),
            "report_interval_steps": int(report_interval_steps),
            "coord_interval_steps": (
                None if coord_interval_steps is None else int(coord_interval_steps)
            ),
            "minimise_first": bool(minimise_first),
            "minimise_max_iterations": int(minimise_max_iterations),
            "minimisation_state": resolved_minimisation_state,
            "umbrella_force_constant_per_kj_mol": float(umbrella_force_constant_per_kj_mol),
            "umbrella_dv_center_kj_mol": float(umbrella_dv_center_kj_mol),
            "random_seed": None if random_seed is None else int(random_seed),
            "platform": platform,
        }

        if output_dir_path is not None:
            summary_path = output_dir_path / "summary.json"
            payload = {"summary": summary, "config": config, "n_trace_rows": len(trace_rows)}
            summary_path.write_text(
                json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
            )
            files["summary_json"] = str(summary_path)

        if cleanup_trace:
            trace_csv_path.unlink(missing_ok=True)

        return {
            "trace": trace_rows,
            "summary": summary,
            "config": config,
            "files": files,
        }

    def run_endpoint_md_comparison(
        self,
        *,
        states: Sequence[str] | None = None,
        equilibration_steps: int = 0,
        production_steps: int = 5000,
        sample_interval: int = 10,
        temperature: unit.Quantity | float = 300.0,
        timestep: unit.Quantity | float = 1.0,
        friction: unit.Quantity | float = 1.0,
        random_seed: int | None = None,
        platform: openmm.Platform | None = None,
        platform_properties: Mapping[str, str] | None = None,
        initial_source_state: str | None = None,
        initial_positions: unit.Quantity | Sequence[openmm.Vec3] | Sequence[Sequence[float]] | None = None,
        initialize_velocities: bool = True,
    ) -> Dict[str, object]:
        """
        Compare EVB endpoint MD against standalone pure-state force-field MD.

        For each requested state, this runs two simulations under identical protocol
        settings (integrator type, thermostat parameters, random seed strategy,
        equilibration/production lengths, and sampling interval):
        ``(a)`` EVB scaffold with indicator lambdas selecting that state, and
        ``(b)`` standalone OpenMM ``System`` built for that state.
        """

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("EVB system not initialized; call 'build_evb_system' and add forces first.")
        if equilibration_steps < 0:
            raise ValueError("'equilibration_steps' must be >= 0.")
        if production_steps <= 0:
            raise ValueError("'production_steps' must be > 0.")
        if sample_interval <= 0:
            raise ValueError("'sample_interval' must be > 0.")

        selected_states: Tuple[str, ...]
        if states is None:
            selected_states = self.state_names
        else:
            selected_states = tuple(states)
            if not selected_states:
                raise ValueError("Provide at least one state name in 'states'.")
            unknown = sorted(set(selected_states) - set(self.state_names))
            if unknown:
                raise ValueError(
                    f"Unknown state(s) requested for endpoint comparison: {tuple(unknown)!r}. "
                    f"Available states: {self.state_names!r}."
                )
            if len(set(selected_states)) != len(selected_states):
                raise ValueError("'states' contains duplicates; provide each state at most once.")

        temperature_q = _units.coerce_temperature(temperature)
        timestep_q = _units.coerce_timestep(timestep)
        friction_q = _units.coerce_friction(friction)
        evb_source_state = self._evb_system_source_state
        coordinate_source_state = initial_source_state or evb_source_state
        if coordinate_source_state not in self.state_names:
            raise ValueError(
                f"Unknown initial source state '{coordinate_source_state}'. Available states: {self.state_names!r}."
            )

        if initial_positions is None:
            source_positions = tuple(self._states[coordinate_source_state].positions)
        else:
            source_positions = self._coerce_vec3_tuple(
                initial_positions,
                quantity_unit=unit.nanometer,
                expected_count=len(self._state_atom_references(coordinate_source_state)),
                label="initial_positions",
            )

        if coordinate_source_state == evb_source_state:
            evb_start_positions = source_positions
        else:
            evb_start_positions = self._remap_positions_for_state(
                source_state=coordinate_source_state,
                target_state=evb_source_state,
                source_positions=source_positions,
            )

        expected_particles = self._evb_system_scaffold.getNumParticles()
        if len(evb_start_positions) != expected_particles:
            raise ValueError(
                "Initial position count does not match the EVB scaffold. "
                f"Expected {expected_particles} positions, received {len(evb_start_positions)}."
            )

        reference_systems = self._get_state_systems()
        evb_topology = self._modellers[evb_source_state].topology
        state_runs: Dict[str, Dict[str, object]] = {}

        for state_index, state_name in enumerate(selected_states):
            state_seed = None if random_seed is None else int(random_seed) + int(state_index)
            state_runs[state_name] = dict(
                run_endpoint_state_comparison(
                    state_name=state_name,
                    state_seed=state_seed,
                    evb_source_state=evb_source_state,
                    evb_topology=evb_topology,
                    evb_system=self._evb_system_scaffold,
                    evb_start_positions=evb_start_positions,
                    expected_particles=expected_particles,
                    reference_topology=self._modellers[state_name].topology,
                    reference_system=reference_systems[state_name],
                    temperature=temperature_q,
                    friction=friction_q,
                    timestep=timestep_q,
                    platform=platform,
                    platform_properties=platform_properties,
                    initialize_velocities=initialize_velocities,
                    equilibration_steps=equilibration_steps,
                    production_steps=production_steps,
                    sample_interval=sample_interval,
                    create_langevin_middle_simulation=self._create_langevin_middle_simulation,
                    attach_scaas_to_simulation=self._attach_scaas_to_simulation,
                    apply_indicator_lambdas=self._apply_indicator_lambdas,
                    coerce_vec3_tuple=self._coerce_vec3_tuple,
                    vec3s_to_quantity=self._vec3s_to_quantity,
                    collect_md_energy_samples=self._collect_md_energy_samples,
                    get_state_lambdas=self.get_state_lambdas,
                    remap_positions_for_state=self._remap_positions_for_state,
                    summarize_md_energy_samples=self._summarize_md_energy_samples,
                )
            )

        protocol = {
            "equilibration_steps": int(equilibration_steps),
            "production_steps": int(production_steps),
            "sample_interval": int(sample_interval),
            "temperature_kelvin": float(temperature_q.value_in_unit(unit.kelvin)),
            "timestep_femtoseconds": float(timestep_q.value_in_unit(unit.femtoseconds)),
            "friction_per_picosecond": float(friction_q.value_in_unit(unit.picosecond**-1)),
            "random_seed": None if random_seed is None else int(random_seed),
            "initialize_velocities": bool(initialize_velocities),
        }
        initial_conditions = {
            "coordinate_source_state": coordinate_source_state,
            "evb_source_state": evb_source_state,
        }
        return {
            "state_names": selected_states,
            "energy_unit": "kJ/mol",
            "protocol": MappingProxyType(protocol),
            "initial_conditions": MappingProxyType(initial_conditions),
            "runs": MappingProxyType(state_runs),
        }

    def collect_mapping_fep_samples(
        self,
        *,
        simulation: app.Simulation | None = None,
        lambda_schedule: Sequence[float] | Sequence[Sequence[float]] | Sequence[Mapping[str, float]] | None = None,
        lambda_start: float = 1.0,
        lambda_end: float = 0.0,
        lambda_step: float = -0.02,
        equilibration_steps: int = 0,
        production_steps: int = 5000,
        sample_interval: int = 10,
        reaction_coordinate_coefficients: Sequence[float] | None = None,
        temperature: unit.Quantity | float = 300.0,
        timestep: unit.Quantity | float = 1.0,
        friction: unit.Quantity | float = 1.0,
        random_seed: int | None = None,
        platform: openmm.Platform | None = None,
        platform_properties: Mapping[str, str] | None = None,
        coordinate_report_interval_ps: float | None = None,
        coordinate_report_file: PathLike | None = None,
        coordinate_report_append: bool = False,
        coordinate_topology_file: PathLike | None = None,
        pre_equilibration: bool = False,
        pre_equilibration_timesteps: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_steps: Sequence[int] | None = None,
        pre_equilibration_temperatures: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_lambda: Iterable[float] | float | Mapping[str, float] | None = None,
        pre_equilibration_lambda_schedule: Sequence[Iterable[float] | float | Mapping[str, float]] | None = None,
        pre_equilibration_friction: Sequence[unit.Quantity | float] | unit.Quantity | float | None = None,
        pre_equilibration_initialize_velocities: str = "first",
        pre_equilibration_random_seed: int | None = None,
        pre_equilibration_energy_abs_cap_kj_mol: float = 1.0e9,
        pre_equilibration_monitor_interval_steps: int = 1,
        pre_equilibration_max_abs_coord_nm_cap: float | None = None,
        pre_equilibration_solute_restraints: bool = False,
        pre_equilibration_solute_restraint_atom_indices: Sequence[int] | None = None,
        pre_equilibration_solute_restraint_k_start: unit.Quantity | float = 2000.0,
        pre_equilibration_solute_restraint_k_end: unit.Quantity | float = 0.0,
        pre_equilibration_solute_restraint_k_schedule: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_solute_restraint_reference_positions: (
            unit.Quantity | Sequence[openmm.Vec3] | Sequence[Sequence[float]] | None
        ) = None,
        pre_equilibration_solute_restraint_include_hydrogens: bool = False,
        pre_equilibration_solute_restraint_solvent_residue_names: Sequence[str] | None = None,
        coupling_distance_atom_pair: tuple[int, int] | None = None,
    ) -> Tuple[FEPWindowRecord, ...]:
        """
        Run consecutive fixed-lambda windows and sample per-state diabatic energies.

        Notes
        -----
        This implements the sampling side of the Q6-like mapping-potential workflow:
        the MD potential is the existing lambda-weighted EVB mapping potential, while
        per-state energies are reconstructed by temporarily toggling indicator lambdas.
        """

        if equilibration_steps < 0:
            raise ValueError("'equilibration_steps' must be >= 0.")
        if production_steps <= 0:
            raise ValueError("'production_steps' must be > 0.")
        if sample_interval <= 0:
            raise ValueError("'sample_interval' must be > 0.")

        windows = self._resolve_lambda_windows(
            lambda_schedule=lambda_schedule,
            lambda_start=lambda_start,
            lambda_end=lambda_end,
            lambda_step=lambda_step,
        )
        coefficients = self._resolve_reaction_coordinate_coefficients(reaction_coordinate_coefficients)
        temperature_q = _units.coerce_temperature(temperature)
        timestep_q = _units.coerce_timestep(timestep)
        friction_q = _units.coerce_friction(friction)

        if self._evb_system_scaffold is None or self._evb_system_source_state is None:
            raise ValueError("EVB system not initialized; call 'build_evb_system' and add forces first.")

        local_simulation = simulation
        if local_simulation is None:
            integrator = openmm.LangevinMiddleIntegrator(temperature_q, friction_q, timestep_q)
            if random_seed is not None:
                integrator.setRandomNumberSeed(int(random_seed))
            topology = self._modellers[self._evb_system_source_state].topology
            if platform is None:
                local_simulation = app.Simulation(topology, self._evb_system_scaffold, integrator)
            else:
                local_simulation = app.Simulation(
                    topology,
                    self._evb_system_scaffold,
                    integrator,
                    platform,
                    dict(platform_properties or {}),
                )
            local_simulation.context.setPositions(self.evb_initial_positions)
            self._attach_scaas_to_simulation(local_simulation, initialize=True)
        else:
            if local_simulation.system.getNumParticles() != self._evb_system_scaffold.getNumParticles():
                raise ValueError(
                    "Provided simulation uses a system with a different number of particles than this EVB scaffold."
                )
            self._attach_scaas_to_simulation(local_simulation, initialize=True)
            if random_seed is not None and hasattr(local_simulation.integrator, "setRandomNumberSeed"):
                try:
                    local_simulation.integrator.setRandomNumberSeed(int(random_seed))
                except openmm.OpenMMException:
                    pass

        if pre_equilibration:
            resolved_pre_lambda = windows[0] if pre_equilibration_lambda is None else pre_equilibration_lambda
            resolved_pre_friction = friction if pre_equilibration_friction is None else pre_equilibration_friction
            resolved_pre_seed = random_seed if pre_equilibration_random_seed is None else pre_equilibration_random_seed
            self.equilibrate_schedule(
                simulation=local_simulation,
                timesteps=pre_equilibration_timesteps,
                steps=pre_equilibration_steps,
                temperatures=pre_equilibration_temperatures,
                equilibration_lambda=resolved_pre_lambda,
                lambda_schedule=pre_equilibration_lambda_schedule,
                friction=resolved_pre_friction,
                initialize_velocities=pre_equilibration_initialize_velocities,
                random_seed=resolved_pre_seed,
                platform=platform,
                platform_properties=platform_properties,
                energy_abs_cap_kj_mol=pre_equilibration_energy_abs_cap_kj_mol,
                max_abs_coord_nm_cap=pre_equilibration_max_abs_coord_nm_cap,
                monitor_interval_steps=pre_equilibration_monitor_interval_steps,
                solute_restraints=pre_equilibration_solute_restraints,
                solute_restraint_atom_indices=pre_equilibration_solute_restraint_atom_indices,
                solute_restraint_k_start=pre_equilibration_solute_restraint_k_start,
                solute_restraint_k_end=pre_equilibration_solute_restraint_k_end,
                solute_restraint_k_schedule=pre_equilibration_solute_restraint_k_schedule,
                solute_restraint_reference_positions=pre_equilibration_solute_restraint_reference_positions,
                solute_restraint_include_hydrogens=pre_equilibration_solute_restraint_include_hydrogens,
                solute_restraint_solvent_residue_names=pre_equilibration_solute_restraint_solvent_residue_names,
                fail_on_instability=True,
            )
            try:
                local_simulation.integrator.setStepSize(timestep_q)
            except (AttributeError, openmm.OpenMMException):
                pass
            try:
                local_simulation.integrator.setTemperature(temperature_q)
            except (AttributeError, openmm.OpenMMException):
                pass
            try:
                local_simulation.integrator.setFriction(friction_q)
            except (AttributeError, openmm.OpenMMException):
                pass
            if random_seed is not None and hasattr(local_simulation.integrator, "setRandomNumberSeed"):
                try:
                    local_simulation.integrator.setRandomNumberSeed(int(random_seed))
                except openmm.OpenMMException:
                    pass

        added_coordinate_reporter: object | None = None
        (
            added_coordinate_reporter,
            _report_interval_steps,
            _coordinate_path,
            _topology_path,
        ) = setup_coordinate_reporting(
            simulation=local_simulation,
            coordinate_report_interval_ps=coordinate_report_interval_ps,
            coordinate_report_file=coordinate_report_file,
            coordinate_report_append=bool(coordinate_report_append),
            coordinate_topology_file=coordinate_topology_file,
            simulation_timestep_ps=lambda simulation: self._simulation_timestep_ps(
                simulation,
                fallback_timestep=timestep_q,
            ),
            steps_from_report_interval_ps=self._steps_from_report_interval_ps,
        )

        # Optional per-frame reactive-pair distance capture, for a distance-dependent
        # off-diagonal coupling H12(r)=A*exp(-mu*(r-r0)) at analysis time. Orthorhombic
        # minimum image is exact here because the made/broken bond is always close.
        coupling_distance_callable: Callable[[app.Simulation], float] | None = None
        if coupling_distance_atom_pair is not None:
            i_atom = int(coupling_distance_atom_pair[0])
            j_atom = int(coupling_distance_atom_pair[1])

            def coupling_distance_callable(sim: app.Simulation) -> float:  # noqa: F811
                state = sim.context.getState(getPositions=True)
                positions = state.getPositions().value_in_unit(unit.nanometer)
                pi = positions[i_atom]
                pj = positions[j_atom]
                dx, dy, dz = pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]
                box = state.getPeriodicBoxVectors()
                if box is not None:
                    edges = [box[k].value_in_unit(unit.nanometer) for k in range(3)]
                    lengths = (edges[0][0], edges[1][1], edges[2][2])
                    deltas = [dx, dy, dz]
                    for k in range(3):
                        length = lengths[k]
                        if length > 0.0:
                            deltas[k] -= length * round(deltas[k] / length)
                    dx, dy, dz = deltas
                return float((dx * dx + dy * dy + dz * dz) ** 0.5)

        try:
            records = collect_fep_window_records(
                simulation=local_simulation,
                windows=windows,
                state_names=self.state_names,
                equilibration_steps=equilibration_steps,
                production_steps=production_steps,
                sample_interval=sample_interval,
                coefficients=coefficients,
                set_state_lambdas=self.set_state_lambdas,
                sample_diabatic_state_energies_kj_mol=self._sample_diabatic_state_energies_kj_mol,
                effective_energy_from_state_energies=self._effective_energy_from_state_energies,
                reaction_coordinate_from_state_energies=self._reaction_coordinate_from_state_energies,
                sample_coupling_distance_nm=coupling_distance_callable,
            )
        finally:
            remove_coordinate_reporter(local_simulation, added_coordinate_reporter)

        return tuple(records)

    def analyze_mapping_fep_umbrella(
        self,
        window_records: Sequence[FEPWindowRecord],
        *,
        temperature: unit.Quantity | float = 300.0,
        skip_samples: int = 0,
        n_gap_bins: int = 40,
        min_points_per_bin: int = 20,
        input_energy_unit: str = "kcal/mol",
        gap_range: tuple[float, float] | None = None,
        gap_percentile_range: tuple[float, float] | None = None,
    ) -> Dict[str, object]:
        """
        Perform Q6-style postprocessing on pre-sampled lambda windows.

        Parameters
        ----------
        input_energy_unit :
            Unit used by ``window_records`` energies/reaction-coordinate values.
            Defaults to ``"kcal/mol"`` for external post-processing data.
            Use ``"kJ/mol"`` when passing records generated by
            :meth:`collect_mapping_fep_samples`.

        Returns
        -------
        Dict[str, object]
            Summary containing forward/reverse Zwanzig free energies and
            gap-binned diabatic free-energy profiles reported in kcal/mol.
        """
        temperature_q = _units.coerce_temperature(temperature)
        kT = float((unit.MOLAR_GAS_CONSTANT_R * temperature_q).value_in_unit(unit.kilojoules_per_mole))
        return analyze_mapping_fep_umbrella_records(
            window_records,
            state_names=self.state_names,
            kT_kj_mol=kT,
            skip_samples=skip_samples,
            n_gap_bins=n_gap_bins,
            min_points_per_bin=min_points_per_bin,
            input_energy_unit=input_energy_unit,
            gap_range=gap_range,
            gap_percentile_range=gap_percentile_range,

        )

    def _write_mapping_fep_outputs(self, analysis: Mapping[str, object], *, output_directory: PathLike) -> None:
        write_mapping_fep_outputs(analysis, output_directory=output_directory)

    def run_mapping_fep_umbrella(
        self,
        *,
        simulation: app.Simulation | None = None,
        lambda_schedule: Sequence[float] | Sequence[Sequence[float]] | Sequence[Mapping[str, float]] | None = None,
        lambda_start: float = 1.0,
        lambda_end: float = 0.0,
        lambda_step: float = -0.02,
        equilibration_steps: int = 0,
        production_steps: int = 5000,
        sample_interval: int = 10,
        skip_samples: int = 0,
        reaction_coordinate_coefficients: Sequence[float] | None = None,
        n_gap_bins: int = 40,
        min_points_per_bin: int = 20,
        temperature: unit.Quantity | float = 300.0,
        timestep: unit.Quantity | float = 1.0,
        friction: unit.Quantity | float = 1.0,
        random_seed: int | None = None,
        platform: openmm.Platform | None = None,
        platform_properties: Mapping[str, str] | None = None,
        output_directory: PathLike | None = None,
        coordinate_report_interval_ps: float | None = None,
        coordinate_report_file: PathLike | None = None,
        coordinate_report_append: bool = False,
        coordinate_topology_file: PathLike | None = None,
        pre_equilibration: bool = False,
        pre_equilibration_timesteps: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_steps: Sequence[int] | None = None,
        pre_equilibration_temperatures: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_lambda: Iterable[float] | float | Mapping[str, float] | None = None,
        pre_equilibration_lambda_schedule: Sequence[Iterable[float] | float | Mapping[str, float]] | None = None,
        pre_equilibration_friction: Sequence[unit.Quantity | float] | unit.Quantity | float | None = None,
        pre_equilibration_initialize_velocities: str = "first",
        pre_equilibration_random_seed: int | None = None,
        pre_equilibration_energy_abs_cap_kj_mol: float = 1.0e9,
        pre_equilibration_monitor_interval_steps: int = 1,
        pre_equilibration_max_abs_coord_nm_cap: float | None = None,
        pre_equilibration_solute_restraints: bool = False,
        pre_equilibration_solute_restraint_atom_indices: Sequence[int] | None = None,
        pre_equilibration_solute_restraint_k_start: unit.Quantity | float = 2000.0,
        pre_equilibration_solute_restraint_k_end: unit.Quantity | float = 0.0,
        pre_equilibration_solute_restraint_k_schedule: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_solute_restraint_reference_positions: (
            unit.Quantity | Sequence[openmm.Vec3] | Sequence[Sequence[float]] | None
        ) = None,
        pre_equilibration_solute_restraint_include_hydrogens: bool = False,
        pre_equilibration_solute_restraint_solvent_residue_names: Sequence[str] | None = None,
    ) -> Dict[str, object]:
        """
        Execute Q6-style mapping-potential FEP windows and analyze them.
        """
        return run_mapping_window_workflow(
            coordinate_report_interval_ps=coordinate_report_interval_ps,
            coordinate_report_file=coordinate_report_file,
            coordinate_topology_file=coordinate_topology_file,
            output_directory=output_directory,
            collect_records=lambda resolved_coordinate_file, resolved_coordinate_topology_file: self.collect_mapping_fep_samples(
                simulation=simulation,
                lambda_schedule=lambda_schedule,
                lambda_start=lambda_start,
                lambda_end=lambda_end,
                lambda_step=lambda_step,
                equilibration_steps=equilibration_steps,
                production_steps=production_steps,
                sample_interval=sample_interval,
                reaction_coordinate_coefficients=reaction_coordinate_coefficients,
                temperature=temperature,
                timestep=timestep,
                friction=friction,
                random_seed=random_seed,
                platform=platform,
                platform_properties=platform_properties,
                coordinate_report_interval_ps=coordinate_report_interval_ps,
                coordinate_report_file=resolved_coordinate_file,
                coordinate_report_append=coordinate_report_append,
                coordinate_topology_file=resolved_coordinate_topology_file,
                pre_equilibration=pre_equilibration,
                pre_equilibration_timesteps=pre_equilibration_timesteps,
                pre_equilibration_steps=pre_equilibration_steps,
                pre_equilibration_temperatures=pre_equilibration_temperatures,
                pre_equilibration_lambda=pre_equilibration_lambda,
                pre_equilibration_lambda_schedule=pre_equilibration_lambda_schedule,
                pre_equilibration_friction=pre_equilibration_friction,
                pre_equilibration_initialize_velocities=pre_equilibration_initialize_velocities,
                pre_equilibration_random_seed=pre_equilibration_random_seed,
                pre_equilibration_energy_abs_cap_kj_mol=pre_equilibration_energy_abs_cap_kj_mol,
                pre_equilibration_monitor_interval_steps=pre_equilibration_monitor_interval_steps,
                pre_equilibration_max_abs_coord_nm_cap=pre_equilibration_max_abs_coord_nm_cap,
                pre_equilibration_solute_restraints=pre_equilibration_solute_restraints,
                pre_equilibration_solute_restraint_atom_indices=pre_equilibration_solute_restraint_atom_indices,
                pre_equilibration_solute_restraint_k_start=pre_equilibration_solute_restraint_k_start,
                pre_equilibration_solute_restraint_k_end=pre_equilibration_solute_restraint_k_end,
                pre_equilibration_solute_restraint_k_schedule=pre_equilibration_solute_restraint_k_schedule,
                pre_equilibration_solute_restraint_reference_positions=pre_equilibration_solute_restraint_reference_positions,
                pre_equilibration_solute_restraint_include_hydrogens=pre_equilibration_solute_restraint_include_hydrogens,
                pre_equilibration_solute_restraint_solvent_residue_names=pre_equilibration_solute_restraint_solvent_residue_names,
            ),
            analyze_records=lambda records: self.analyze_mapping_fep_umbrella(
                records,
                temperature=temperature,
                skip_samples=skip_samples,
                n_gap_bins=n_gap_bins,
                min_points_per_bin=min_points_per_bin,
                input_energy_unit="kJ/mol",
                gap_range=gap_range,
                gap_percentile_range=gap_percentile_range,
            ),
            write_outputs=(
                lambda analysis, outdir: self._write_mapping_fep_outputs(analysis, output_directory=outdir),
            ),
        )

    def analyze_evb_profile(
        self,
        window_records: Sequence[FEPWindowRecord],
        *,
        diabatic_offsets_kcal_mol: Sequence[float] | None = None,
        coupling_matrix_kcal_mol: Sequence[Sequence[float]] | None = None,
        reaction_coordinate_states: Tuple[int, int] = (0, 1),
        temperature: unit.Quantity | float = 300.0,
        skip_samples: int = 0,
        n_gap_bins: int = 40,
        min_points_per_bin: int = 20,
        input_energy_unit: str = "kcal/mol",
        gap_range: tuple[float, float] | None = None,
        gap_percentile_range: tuple[float, float] | None = None,
        coupling_form: str = "constant",
        coupling_width_matrix_per_kcal2_mol2: Sequence[Sequence[float]] | None = None,
        coupling_decay_per_nm: float = 0.0,
        coupling_reference_distance_nm: float = 0.0,
    ) -> Dict[str, object]:
        """EVB profile analysis for an arbitrary number of diabatic states.

        Thin wrapper over
        :func:`evbopenmm.analysis.analyze_evb_profile_records` that first
        builds the FEP base analysis via :meth:`analyze_mapping_fep_umbrella`.
        For two-state runs pass ``diabatic_offsets_kcal_mol=(0.0, alpha)`` and
        a 2x2 ``coupling_matrix_kcal_mol``.

        ``coupling_form`` selects the off-diagonal element model:
        ``"constant"`` (default) uses ``coupling_matrix_kcal_mol`` directly;
        ``"gaussian"`` treats that matrix as amplitudes A and requires
        ``coupling_width_matrix_per_kcal2_mol2`` for the Gaussian widths mu, so
        ``H_ij(dE) = A_ij * exp(-mu_ij * (E_i - E_j)^2)`` — the classic
        Warshel/Aqvist form.
        """

        base_analysis = self.analyze_mapping_fep_umbrella(
            window_records,
            temperature=temperature,
            skip_samples=skip_samples,
            n_gap_bins=n_gap_bins,
            min_points_per_bin=min_points_per_bin,
            input_energy_unit=input_energy_unit,
            gap_range=gap_range,
            gap_percentile_range=gap_percentile_range,

        )
        return analyze_evb_profile_records(
            window_records,
            state_names=self.state_names,
            base_analysis=base_analysis,
            diabatic_offsets_kcal_mol=diabatic_offsets_kcal_mol,
            coupling_matrix_kcal_mol=coupling_matrix_kcal_mol,
            reaction_coordinate_states=reaction_coordinate_states,
            skip_samples=skip_samples,
            n_gap_bins=n_gap_bins,
            min_points_per_bin=min_points_per_bin,
            input_energy_unit=input_energy_unit,
            gap_range=gap_range,
            gap_percentile_range=gap_percentile_range,

            coupling_form=coupling_form,
            coupling_width_matrix_per_kcal2_mol2=coupling_width_matrix_per_kcal2_mol2,
            coupling_decay_per_nm=coupling_decay_per_nm,
            coupling_reference_distance_nm=coupling_reference_distance_nm,
        )

    def run_mapping_fep_evb_profile(
        self,
        *,
        simulation: app.Simulation | None = None,
        lambda_schedule: Sequence[float] | Sequence[Sequence[float]] | Sequence[Mapping[str, float]] | None = None,
        lambda_start: float = 1.0,
        lambda_end: float = 0.0,
        lambda_step: float = -0.02,
        equilibration_steps: int = 0,
        production_steps: int = 5000,
        sample_interval: int = 10,
        skip_samples: int = 0,
        reaction_coordinate_coefficients: Sequence[float] | None = None,
        n_gap_bins: int = 40,
        min_points_per_bin: int = 20,
        gap_range: tuple[float, float] | None = None,
        gap_percentile_range: tuple[float, float] | None = None,
        alpha_kcal_mol: float = 0.0,
        hij_kcal_mol: float = 0.0,
        temperature: unit.Quantity | float = 300.0,
        timestep: unit.Quantity | float = 1.0,
        friction: unit.Quantity | float = 1.0,
        random_seed: int | None = None,
        platform: openmm.Platform | None = None,
        platform_properties: Mapping[str, str] | None = None,
        output_directory: PathLike | None = None,
        coordinate_report_interval_ps: float | None = None,
        coordinate_report_file: PathLike | None = None,
        coordinate_report_append: bool = False,
        coordinate_topology_file: PathLike | None = None,
        pre_equilibration: bool = False,
        pre_equilibration_timesteps: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_steps: Sequence[int] | None = None,
        pre_equilibration_temperatures: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_lambda: Iterable[float] | float | Mapping[str, float] | None = None,
        pre_equilibration_lambda_schedule: Sequence[Iterable[float] | float | Mapping[str, float]] | None = None,
        pre_equilibration_friction: Sequence[unit.Quantity | float] | unit.Quantity | float | None = None,
        pre_equilibration_initialize_velocities: str = "first",
        pre_equilibration_random_seed: int | None = None,
        pre_equilibration_energy_abs_cap_kj_mol: float = 1.0e9,
        pre_equilibration_monitor_interval_steps: int = 1,
        pre_equilibration_max_abs_coord_nm_cap: float | None = None,
        pre_equilibration_solute_restraints: bool = False,
        pre_equilibration_solute_restraint_atom_indices: Sequence[int] | None = None,
        pre_equilibration_solute_restraint_k_start: unit.Quantity | float = 2000.0,
        pre_equilibration_solute_restraint_k_end: unit.Quantity | float = 0.0,
        pre_equilibration_solute_restraint_k_schedule: Sequence[unit.Quantity | float] | None = None,
        pre_equilibration_solute_restraint_reference_positions: (
            unit.Quantity | Sequence[openmm.Vec3] | Sequence[Sequence[float]] | None
        ) = None,
        pre_equilibration_solute_restraint_include_hydrogens: bool = False,
        pre_equilibration_solute_restraint_solvent_residue_names: Sequence[str] | None = None,
    ) -> Dict[str, object]:
        """
        Execute two-state mapping-potential windows and return EVB profile analysis.
        """
        if len(self.state_names) != 2:
            raise ValueError(
                "run_mapping_fep_evb_profile requires a two-state EVB system. "
                "For N != 2 use the mapping FEP primitives + analyze_evb_profile directly."
            )
        return run_mapping_window_workflow(
            coordinate_report_interval_ps=coordinate_report_interval_ps,
            coordinate_report_file=coordinate_report_file,
            coordinate_topology_file=coordinate_topology_file,
            output_directory=output_directory,
            collect_records=lambda resolved_coordinate_file, resolved_coordinate_topology_file: self.collect_mapping_fep_samples(
                simulation=simulation,
                lambda_schedule=lambda_schedule,
                lambda_start=lambda_start,
                lambda_end=lambda_end,
                lambda_step=lambda_step,
                equilibration_steps=equilibration_steps,
                production_steps=production_steps,
                sample_interval=sample_interval,
                reaction_coordinate_coefficients=reaction_coordinate_coefficients,
                temperature=temperature,
                timestep=timestep,
                friction=friction,
                random_seed=random_seed,
                platform=platform,
                platform_properties=platform_properties,
                coordinate_report_interval_ps=coordinate_report_interval_ps,
                coordinate_report_file=resolved_coordinate_file,
                coordinate_report_append=coordinate_report_append,
                coordinate_topology_file=resolved_coordinate_topology_file,
                pre_equilibration=pre_equilibration,
                pre_equilibration_timesteps=pre_equilibration_timesteps,
                pre_equilibration_steps=pre_equilibration_steps,
                pre_equilibration_temperatures=pre_equilibration_temperatures,
                pre_equilibration_lambda=pre_equilibration_lambda,
                pre_equilibration_lambda_schedule=pre_equilibration_lambda_schedule,
                pre_equilibration_friction=pre_equilibration_friction,
                pre_equilibration_initialize_velocities=pre_equilibration_initialize_velocities,
                pre_equilibration_random_seed=pre_equilibration_random_seed,
                pre_equilibration_energy_abs_cap_kj_mol=pre_equilibration_energy_abs_cap_kj_mol,
                pre_equilibration_monitor_interval_steps=pre_equilibration_monitor_interval_steps,
                pre_equilibration_max_abs_coord_nm_cap=pre_equilibration_max_abs_coord_nm_cap,
                pre_equilibration_solute_restraints=pre_equilibration_solute_restraints,
                pre_equilibration_solute_restraint_atom_indices=pre_equilibration_solute_restraint_atom_indices,
                pre_equilibration_solute_restraint_k_start=pre_equilibration_solute_restraint_k_start,
                pre_equilibration_solute_restraint_k_end=pre_equilibration_solute_restraint_k_end,
                pre_equilibration_solute_restraint_k_schedule=pre_equilibration_solute_restraint_k_schedule,
                pre_equilibration_solute_restraint_reference_positions=pre_equilibration_solute_restraint_reference_positions,
                pre_equilibration_solute_restraint_include_hydrogens=pre_equilibration_solute_restraint_include_hydrogens,
                pre_equilibration_solute_restraint_solvent_residue_names=pre_equilibration_solute_restraint_solvent_residue_names,
            ),
            analyze_records=lambda records: self.analyze_evb_profile(
                records,
                diabatic_offsets_kcal_mol=(0.0, float(alpha_kcal_mol)),
                coupling_matrix_kcal_mol=(
                    (0.0, abs(float(hij_kcal_mol))),
                    (abs(float(hij_kcal_mol)), 0.0),
                ),
                reaction_coordinate_states=(0, 1),
                temperature=temperature,
                skip_samples=skip_samples,
                n_gap_bins=n_gap_bins,
                min_points_per_bin=min_points_per_bin,
                input_energy_unit="kJ/mol",
            ),
            write_outputs=(
                lambda analysis, outdir: self._write_mapping_fep_outputs(analysis, output_directory=outdir),
                lambda analysis, outdir: self._write_evb_profile_outputs(analysis, output_directory=outdir),
            ),
        )

    def aggregate_evb_profiles(
        self,
        analyses: Sequence[Mapping[str, object]],
        *,
        minimum_replicas_per_bin: int = 1,
        gap_tolerance: float = 1e-8,
    ) -> Dict[str, object]:
        """Aggregate replica EVB profile analyses across an arbitrary number of diabatic states."""
        return aggregate_evb_profiles(
            analyses,
            minimum_replicas_per_bin=minimum_replicas_per_bin,
            gap_tolerance=gap_tolerance,
        )

    def _write_evb_profile_outputs(self, analysis: Mapping[str, object], *, output_directory: PathLike) -> None:
        write_evb_profile_outputs(analysis, output_directory=output_directory)

    def _write_evb_replica_outputs(self, aggregation: Mapping[str, object], *, output_directory: PathLike) -> None:
        write_evb_replica_outputs(aggregation, output_directory=output_directory)

    def analyze_evb_profile_replicas(
        self,
        replica_window_records: Sequence[Sequence[FEPWindowRecord]],
        *,
        diabatic_offsets_kcal_mol: Sequence[float] | None = None,
        coupling_matrix_kcal_mol: Sequence[Sequence[float]] | None = None,
        reaction_coordinate_states: Tuple[int, int] = (0, 1),
        temperature: unit.Quantity | float = 300.0,
        skip_samples: int = 0,
        n_gap_bins: int = 40,
        min_points_per_bin: int = 20,
        minimum_replicas_per_bin: int = 1,
        input_energy_unit: str = "kcal/mol",
        gap_range: tuple[float, float] | None = None,
        gap_percentile_range: tuple[float, float] | None = None,
        output_directory: PathLike | None = None,
    ) -> Dict[str, object]:
        """Replica EVB profile analyzer for an arbitrary number of diabatic states.

        Runs :meth:`analyze_evb_profile` on each replica, then averages the
        resulting profiles with :meth:`aggregate_evb_profiles`.
        """
        return analyze_evb_profile_replicas_result(
            replica_window_records,
            analyze_profile=lambda records: self.analyze_evb_profile(
                records,
                diabatic_offsets_kcal_mol=diabatic_offsets_kcal_mol,
                coupling_matrix_kcal_mol=coupling_matrix_kcal_mol,
                reaction_coordinate_states=reaction_coordinate_states,
                temperature=temperature,
                skip_samples=skip_samples,
                n_gap_bins=n_gap_bins,
                min_points_per_bin=min_points_per_bin,
                input_energy_unit=input_energy_unit,
                gap_range=gap_range,
                gap_percentile_range=gap_percentile_range,
            ),
            aggregate_profiles=self.aggregate_evb_profiles,
            minimum_replicas_per_bin=minimum_replicas_per_bin,
            output_directory=output_directory,
            write_mapping_fep_outputs=lambda analysis, outdir: self._write_mapping_fep_outputs(
                analysis,
                output_directory=outdir,
            ),
            write_evb_profile_outputs=lambda analysis, outdir: self._write_evb_profile_outputs(
                analysis,
                output_directory=outdir,
            ),
            write_replica_outputs=lambda aggregation, outdir: self._write_evb_replica_outputs(
                aggregation,
                output_directory=outdir,
            ),
        )

    def plot_evb_profile_replicas(
        self,
        replica_result: Mapping[str, object],
        *,
        profile_block: str = "state_free_energies_kcal_mol",
        curves: Sequence[str] | None = None,
        include_aggregate_mean: bool = True,
        replica_alpha: float = 0.45,
        replica_linewidth: float = 1.1,
        mean_linewidth: float = 2.2,
        figure_size: Tuple[float, float] | None = None,
        title: str | None = None,
        show_legend: bool = True,
    ) -> Dict[str, object]:
        """Plot per-replica EVB profile curves for an arbitrary number of diabatic states."""
        return plot_evb_profile_replicas(
            replica_result,
            profile_block=profile_block,
            curves=curves,
            include_aggregate_mean=include_aggregate_mean,
            replica_alpha=replica_alpha,
            replica_linewidth=replica_linewidth,
            mean_linewidth=mean_linewidth,
            figure_size=figure_size,
            title=title,
            show_legend=show_legend,
        )

    def _build_evb_profile_candidate(
        self,
        window_records: Sequence[FEPWindowRecord],
        *,
        reaction_coordinate_states: Tuple[int, int],
        temperature: unit.Quantity | float,
        skip_samples: int,
        n_gap_bins: int,
        min_points_per_bin: int,
        input_energy_unit: str,
    ) -> Callable[[Mapping[Tuple[str, str], float], Mapping[str, float]], Dict[str, object]]:
        state_names = self.state_names
        n_states = len(state_names)
        state_to_index = {name: idx for idx, name in enumerate(state_names)}

        def _analyze(
            couplings: Mapping[Tuple[str, str], float],
            offsets: Mapping[str, float],
        ) -> Dict[str, object]:
            offset_vec = [0.0] * n_states
            for name, value in offsets.items():
                offset_vec[state_to_index[name]] = float(value)
            coupling_mat = [[0.0] * n_states for _ in range(n_states)]
            for (i, j), value in couplings.items():
                a = state_to_index[i]
                b = state_to_index[j]
                coupling_mat[a][b] = float(value)
                coupling_mat[b][a] = float(value)
            return self.analyze_evb_profile(
                window_records,
                diabatic_offsets_kcal_mol=tuple(offset_vec),
                coupling_matrix_kcal_mol=tuple(tuple(row) for row in coupling_mat),
                reaction_coordinate_states=reaction_coordinate_states,
                temperature=temperature,
                skip_samples=skip_samples,
                n_gap_bins=n_gap_bins,
                min_points_per_bin=min_points_per_bin,
                input_energy_unit=input_energy_unit,
            )

        return _analyze

    def fit_evb_parameters(
        self,
        window_records: Sequence[FEPWindowRecord],
        *,
        coupling_pattern: Sequence[Tuple[str, str]],
        target_dg_dagger_kcal_mol: float | None = None,
        target_dg0_kcal_mol: float | None = None,
        weight_dg_dagger: float = 1.0,
        weight_dg0: float = 1.0,
        diabatic_offset_pattern: Sequence[str] | None = None,
        initial_couplings_kcal_mol: Mapping[Tuple[str, str], float] | None = None,
        initial_offsets_kcal_mol: Mapping[str, float] | None = None,
        coupling_bounds_kcal_mol: Tuple[float, float] | Mapping[Tuple[str, str], Tuple[float, float]] = (-20.0, 20.0),
        offset_bounds_kcal_mol: Tuple[float, float] | Mapping[str, Tuple[float, float]] = (-50.0, 50.0),
        reaction_coordinate_states: Tuple[int, int] | None = None,
        temperature: unit.Quantity | float = 300.0,
        skip_samples: int = 0,
        n_gap_bins: int = 40,
        min_points_per_bin: int = 20,
        input_energy_unit: str = "kcal/mol",
        scipy_method: str = "Powell",
        scipy_tolerance: float | None = None,
        scipy_options: Mapping[str, object] | None = None,
    ) -> Dict[str, object]:
        """N-state EVB parameter fit against ``dG_dagger`` / ``dG0`` targets.

        ``coupling_pattern`` lists the ``(state_i, state_j)`` pairs whose
        ``H_{ij}`` you want to vary; every unlisted off-diagonal stays at 0.
        For a vanilla two-state run, pass
        ``coupling_pattern=[(reactant_name, product_name)]``.

        Offsets for the states in ``diabatic_offset_pattern`` (default: every
        state except ``state_names[0]``) are also varied. The first state is
        always the energy anchor.
        """
        state_names = self.state_names
        if diabatic_offset_pattern is None:
            diabatic_offset_pattern = state_names[1:]
        if reaction_coordinate_states is None:
            reaction_coordinate_states = (0, len(state_names) - 1)

        analyze_candidate = self._build_evb_profile_candidate(
            window_records,
            reaction_coordinate_states=reaction_coordinate_states,
            temperature=temperature,
            skip_samples=skip_samples,
            n_gap_bins=n_gap_bins,
            min_points_per_bin=min_points_per_bin,
            input_energy_unit=input_energy_unit,
        )
        return fit_evb_parameters_nstate(
            state_names=state_names,
            coupling_pattern=coupling_pattern,
            diabatic_offset_pattern=diabatic_offset_pattern,
            analyze_candidate=analyze_candidate,
            target_dg_dagger_kcal_mol=target_dg_dagger_kcal_mol,
            target_dg0_kcal_mol=target_dg0_kcal_mol,
            weight_dg_dagger=weight_dg_dagger,
            weight_dg0=weight_dg0,
            initial_couplings_kcal_mol=initial_couplings_kcal_mol,
            initial_offsets_kcal_mol=initial_offsets_kcal_mol,
            coupling_bounds_kcal_mol=coupling_bounds_kcal_mol,
            offset_bounds_kcal_mol=offset_bounds_kcal_mol,
            scipy_method=scipy_method,
            scipy_tolerance=scipy_tolerance,
            scipy_options=scipy_options,
        )

    def _build_evb_profile_replicas_candidate(
        self,
        replica_window_records: Sequence[Sequence[FEPWindowRecord]],
        *,
        reaction_coordinate_states: Tuple[int, int],
        temperature: unit.Quantity | float,
        skip_samples: int,
        n_gap_bins: int,
        min_points_per_bin: int,
        minimum_replicas_per_bin: int,
        input_energy_unit: str,
    ) -> Callable[[Mapping[Tuple[str, str], float], Mapping[str, float]], Dict[str, object]]:
        state_names = self.state_names
        n_states = len(state_names)
        state_to_index = {name: idx for idx, name in enumerate(state_names)}

        def _analyze(
            couplings: Mapping[Tuple[str, str], float],
            offsets: Mapping[str, float],
        ) -> Dict[str, object]:
            offset_vec = [0.0] * n_states
            for name, value in offsets.items():
                offset_vec[state_to_index[name]] = float(value)
            coupling_mat = [[0.0] * n_states for _ in range(n_states)]
            for (i, j), value in couplings.items():
                a = state_to_index[i]
                b = state_to_index[j]
                coupling_mat[a][b] = float(value)
                coupling_mat[b][a] = float(value)
            return self.analyze_evb_profile_replicas(
                replica_window_records,
                diabatic_offsets_kcal_mol=tuple(offset_vec),
                coupling_matrix_kcal_mol=tuple(tuple(row) for row in coupling_mat),
                reaction_coordinate_states=reaction_coordinate_states,
                temperature=temperature,
                skip_samples=skip_samples,
                n_gap_bins=n_gap_bins,
                min_points_per_bin=min_points_per_bin,
                minimum_replicas_per_bin=minimum_replicas_per_bin,
                input_energy_unit=input_energy_unit,
                output_directory=None,
            )

        return _analyze

    def fit_evb_parameters_replicas(
        self,
        replica_window_records: Sequence[Sequence[FEPWindowRecord]],
        *,
        coupling_pattern: Sequence[Tuple[str, str]],
        target_dg_dagger_kcal_mol: float | None = None,
        target_dg0_kcal_mol: float | None = None,
        weight_dg_dagger: float = 1.0,
        weight_dg0: float = 1.0,
        diabatic_offset_pattern: Sequence[str] | None = None,
        initial_couplings_kcal_mol: Mapping[Tuple[str, str], float] | None = None,
        initial_offsets_kcal_mol: Mapping[str, float] | None = None,
        coupling_bounds_kcal_mol: Tuple[float, float] | Mapping[Tuple[str, str], Tuple[float, float]] = (-20.0, 20.0),
        offset_bounds_kcal_mol: Tuple[float, float] | Mapping[str, Tuple[float, float]] = (-50.0, 50.0),
        reaction_coordinate_states: Tuple[int, int] | None = None,
        temperature: unit.Quantity | float = 300.0,
        skip_samples: int = 0,
        n_gap_bins: int = 40,
        min_points_per_bin: int = 20,
        minimum_replicas_per_bin: int = 1,
        input_energy_unit: str = "kcal/mol",
        scipy_method: str = "Powell",
        scipy_tolerance: float | None = None,
        scipy_options: Mapping[str, object] | None = None,
    ) -> Dict[str, object]:
        """Replica-ensemble variant of :meth:`fit_evb_parameters`.

        Objective uses the aggregate ``dG_dagger`` / ``dG0`` means across
        replicas rather than a single-run value.
        """
        if not replica_window_records:
            raise ValueError("At least one replica window-record sequence is required.")
        state_names = self.state_names
        if diabatic_offset_pattern is None:
            diabatic_offset_pattern = state_names[1:]
        if reaction_coordinate_states is None:
            reaction_coordinate_states = (0, len(state_names) - 1)

        analyze_candidate = self._build_evb_profile_replicas_candidate(
            replica_window_records,
            reaction_coordinate_states=reaction_coordinate_states,
            temperature=temperature,
            skip_samples=skip_samples,
            n_gap_bins=n_gap_bins,
            min_points_per_bin=min_points_per_bin,
            minimum_replicas_per_bin=minimum_replicas_per_bin,
            input_energy_unit=input_energy_unit,
        )
        return fit_evb_parameters_replicas_nstate(
            state_names=state_names,
            replica_count=len(replica_window_records),
            coupling_pattern=coupling_pattern,
            diabatic_offset_pattern=diabatic_offset_pattern,
            analyze_candidate=analyze_candidate,
            target_dg_dagger_kcal_mol=target_dg_dagger_kcal_mol,
            target_dg0_kcal_mol=target_dg0_kcal_mol,
            weight_dg_dagger=weight_dg_dagger,
            weight_dg0=weight_dg0,
            initial_couplings_kcal_mol=initial_couplings_kcal_mol,
            initial_offsets_kcal_mol=initial_offsets_kcal_mol,
            coupling_bounds_kcal_mol=coupling_bounds_kcal_mol,
            offset_bounds_kcal_mol=offset_bounds_kcal_mol,
            scipy_method=scipy_method,
            scipy_tolerance=scipy_tolerance,
            scipy_options=scipy_options,
        )

    def get_energy_breakdown(self, context: openmm.Context) -> Dict[str, object]:
        group_energies = self._collect_group_energies(context)
        return build_evb_energy_breakdown(
            force_group_ids=self._force_group_ids,
            group_energies=group_energies,
            state_names=self.state_names,
        )

    def average_energy_breakdown(
        self,
        context: openmm.Context,
        frames: "Iterable[Sequence | Tuple[Sequence, object]]",
    ) -> Dict[str, object]:
        """Ensemble-average :meth:`get_energy_breakdown` over a set of coordinate frames.

        The per-frame :meth:`get_energy_breakdown` is the building block for an EVB
        energy-component (per-state, per-term) analysis; this aggregates it over an
        ensemble (e.g. trajectory frames) so the caller does not have to loop by hand.
        The ``context``'s current lambda values are used unchanged for every frame, so
        set the desired diabatic state (e.g. ``set_state_lambdas(sim, 1.0)`` for the
        reactant-pure diabat) before calling. See :meth:`diabatic_energy_components`
        for the full two-ensemble decomposition.

        Parameters
        ----------
        context :
            A Context running this EVB system's scaffold.
        frames :
            Iterable of either ``positions`` or ``(positions, box_vectors)``. Positions
            are passed to ``context.setPositions``; box vectors, when given, to
            ``context.setPeriodicBoxVectors`` (required for periodic/PME systems whose
            box fluctuates between frames). Quantities or bare arrays are accepted.

        Returns
        -------
        Dict[str, object]
            The same nested structure as :meth:`get_energy_breakdown`, with every numeric
            energy leaf replaced by its mean over the supplied frames, plus a top-level
            ``"n_frames"`` count.
        """
        frame_list = list(frames)
        if not frame_list:
            raise ValueError("average_energy_breakdown requires at least one frame.")

        accumulator: Dict[str, object] | None = None
        for frame in frame_list:
            positions, box = self._coerce_breakdown_frame(frame)
            if box is not None:
                context.setPeriodicBoxVectors(*box)
            context.setPositions(positions)
            breakdown = self.get_energy_breakdown(context)
            accumulator = (
                breakdown if accumulator is None else accumulate_energy_breakdown(accumulator, breakdown)
            )

        averaged = scale_energy_breakdown(accumulator, 1.0 / len(frame_list))
        averaged["n_frames"] = len(frame_list)
        return averaged

    def diabatic_energy_components(
        self,
        context: openmm.Context,
        frames_per_state: "Mapping[str, Iterable]",
    ) -> Dict[str, object]:
        """Per-(state, term) energy-component decomposition over per-state ensembles.

        For each ``state_name -> frames`` in ``frames_per_state``, set that diabat to its
        pure-state lambda (``lambda_state=1``, others 0) and ensemble-average the energy
        breakdown over its frames. Returns ``{state_name: averaged_breakdown}`` plus a
        ``"delta_terms"`` table giving, for a 2-state system, the per-term difference
        ``<E_term>_(second state) - <E_term>_(first state)`` in kJ/mol -- the building
        block for localising which force-group term drives ``dG0``/``dG_dagger``.

        ``frames_per_state`` keys must be a subset of :attr:`state_names`. Lambda values
        are restored to their prior settings on return.
        """
        unknown = set(frames_per_state) - set(self.state_names)
        if unknown:
            raise ValueError(
                f"frames_per_state references unknown states {sorted(unknown)}; "
                f"valid states are {list(self.state_names)}."
            )
        lam = {name: (self._state_lambda_parameter_names.get(name) or f"lambda_{name}")
               for name in self.state_names}
        saved = {name: self._global_parameter_value(context, lam[name]) for name in self.state_names}
        result: Dict[str, object] = {}
        try:
            for state_name, frames in frames_per_state.items():
                for name in self.state_names:
                    self._set_global_parameter(context, lam[name], 1.0 if name == state_name else 0.0)
                result[state_name] = self.average_energy_breakdown(context, frames)
        finally:
            for name, value in saved.items():
                if value is not None:
                    self._set_global_parameter(context, lam[name], value)

        states = list(frames_per_state)
        if len(states) == 2:
            a, b = states
            terms_a = result[a].get("states", {}).get(a, {})
            terms_b = result[b].get("states", {}).get(b, {})
            delta = {}
            for term in set(terms_a) | set(terms_b):
                ea = energy_value_kj(terms_a.get(term))
                eb = energy_value_kj(terms_b.get(term))
                delta[term] = eb - ea
            result["delta_terms"] = {"from": a, "to": b, "kj_per_mol": delta}
        return result

    @staticmethod
    def _coerce_breakdown_frame(frame):
        """Split a frame into (positions, box-or-None); accept positions or (positions, box)."""
        if isinstance(frame, tuple) and len(frame) == 2 and not unit.is_quantity(frame):
            return frame[0], frame[1]
        return frame, None

    @staticmethod
    def _global_parameter_value(context: openmm.Context, name: str):
        try:
            return context.getParameter(name)
        except Exception:
            return None

    @staticmethod
    def _set_global_parameter(context: openmm.Context, name: str, value: float) -> None:
        try:
            context.setParameter(name, float(value))
        except Exception:
            pass

    def debug_state_energy_differences(
        self,
        context: openmm.Context,
        *,
        system_kwargs: Mapping[str, object] | None = None,
        platform: openmm.Platform | None = None,
        print_summary: bool = True,
    ) -> Dict[str, object]:
        """
        Compare EVB per-state energies against standalone OpenMM systems.

        Parameters
        ----------
        context :
            Context running the EVB system. Its coordinates are reused for all
            reference state evaluations.
        system_kwargs :
            Optional ``ForceField.createSystem`` overrides used when rebuilding the
            pure-state reference systems (e.g., PME settings). Defaults to the cached
            configuration employed when the EVB forces were constructed.
        platform :
            Optional OpenMM platform used when evaluating the reference systems.
        print_summary :
            When ``True`` (default) a human-readable energy table is printed for quick inspection.

        Returns
        -------
        Dict[str, object]
            Mapping with keys:

            ``evb_breakdown`` :
                The result of :meth:`get_energy_breakdown`.
            ``states`` :
                Per-state dictionaries containing the EVB term energies, the reference
                term energies, term-by-term deltas, any \"other\" force contributions,
                and the total discrepancy.
            ``reference_grand_total`` :
                Sum of all pure-state totals.
            ``evb_state_total`` :
                Sum of all EVB state totals (shared forces excluded).
            ``evb_state_total_with_shared`` :
                Sum of EVB pure-state totals including the shared contribution once per state.
            ``shared_term_estimates`` :
                Estimated per-term shared energies inferred from ``reference_terms - evb_terms``
                across states.
            ``shared_term_spread`` :
                Per-term spread (max absolute deviation) of the above shared-energy estimates
                across states, in kJ/mol.
            ``shared_term_estimate_residual`` :
                Difference between ``shared_total`` and the sum of ``shared_term_estimates``.
            ``state_total_delta`` :
                Difference ``reference_grand_total - evb_state_total_with_shared``.

        Notes
        -----
        For each EVB state, the context's lambda parameters are temporarily set to a
        pure-state configuration (``lambda_state=1``, rest ``0``) before collecting
        energy terms. Lambda values are restored to their original settings at the end.
        Reference systems are reused directly: force groups are reassigned temporarily
        during the comparison and restored before returning. Nonbonded parameters are
        likewise reset after measuring van der Waals and electrostatic contributions.
        

        .. warning::
           Comparisons are only meaningful when the EVB scaffold and the
           per-state reference System use **matching nonbonded settings**. The
           scaffold rebuilds nonbonded interactions as cutoff-based custom
           forces, so comparing it against a reference built with PME -- or
           with different MM and QM cutoffs -- measures that mismatch, not the
           diabat force field. On a solvated system the artefact reaches ~1e4
           kcal/mol and appears on *every* state, including unmodified ones,
           which is the tell. Build the reference with the same nonbonded
           method and a single uniform cutoff (or compare in gas phase).
        """

        if self._evb_system_source_state is None:
            raise ValueError("EVB umbrella system not initialized; call 'build_evb_system' first.")

        system_options = dict(system_kwargs or {})
        reference_systems = self.build_state_systems(**system_options)

        original_lambdas = self._context_lambda_snapshot(context)
        evb_terms_by_state: Dict[str, Dict[str, unit.Quantity]] = {}
        shared_total: unit.Quantity | None = None
        for state_name in self.state_names:
            self._apply_indicator_lambdas(context, active_state=state_name)
            state_breakdown = self.get_energy_breakdown(context)
            if shared_total is None:
                shared_total = state_breakdown["shared_total"]
            evb_terms_by_state[state_name] = state_breakdown["states"][state_name]
        self._apply_context_lambdas(context, original_lambdas)
        if shared_total is None:
            shared_total = self._zero_energy()

        source_positions_quantity = context.getState(getPositions=True).getPositions()
        source_positions_nm = source_positions_quantity.value_in_unit(unit.nanometer)
        source_positions = tuple(openmm.Vec3(float(vec[0]), float(vec[1]), float(vec[2])) for vec in source_positions_nm)

        state_results: Dict[str, Dict[str, object]] = {}
        reference_total = self._zero_energy()

        for state_name in self.state_names:
            reference_system = reference_systems[state_name]
            group_metadata, previous_groups = self._assign_reference_force_groups(reference_system)

            integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
            ref_context: openmm.Context | None = None
            try:
                if platform is None:
                    ref_context = openmm.Context(reference_system, integrator)
                else:
                    ref_context = openmm.Context(reference_system, integrator, platform)

                remapped_positions = self._remap_positions_for_state(
                    source_state=self._evb_system_source_state,
                    target_state=state_name,
                    source_positions=source_positions,
                )
                ref_context.setPositions(self._vec3s_to_quantity(remapped_positions))

                group_energies: Dict[int, unit.Quantity] = {}
                for group_id in group_metadata:
                    state = ref_context.getState(getEnergy=True, groups=1 << group_id)
                    group_energies[group_id] = state.getPotentialEnergy()

                reference_terms, other_contributions = self._reduce_reference_terms(
                    ref_context, group_metadata, group_energies
                )
            finally:
                if ref_context is not None:
                    del ref_context
                del integrator
                self._restore_reference_force_groups(reference_system, previous_groups)

            evb_terms = evb_terms_by_state[state_name]
            term_deltas: Dict[str, unit.Quantity] = {}
            for term in STATE_ENERGY_TERMS:
                term_deltas[term] = reference_terms[term] - evb_terms.get(term, self._zero_energy())

            ref_total = sum(reference_terms.values(), self._zero_energy()) + sum(
                other_contributions.values(), self._zero_energy()
            )
            evb_total = evb_terms.get("total", self._zero_energy())
            evb_total_with_shared = evb_total + shared_total
            state_results[state_name] = {
                "reference_terms": reference_terms,
                "reference_other_forces": other_contributions,
                "evb_terms": evb_terms,
                "term_deltas": term_deltas,
                "reference_total": ref_total,
                "evb_total": evb_total,
                "evb_total_with_shared": evb_total_with_shared,
                "total_delta": ref_total - evb_total_with_shared,
            }
            reference_total = reference_total + ref_total

        shared_term_estimates: Dict[str, unit.Quantity] = {}
        shared_term_spread: Dict[str, float] = {}
        n_states = len(self.state_names)
        for term in STATE_ENERGY_TERMS:
            values = [state_results[state_name]["term_deltas"][term] for state_name in self.state_names]
            estimate = sum(values, self._zero_energy()) / n_states
            estimate_kj = estimate.value_in_unit(unit.kilojoules_per_mole)
            spread_kj = max(
                abs(value.value_in_unit(unit.kilojoules_per_mole) - estimate_kj)
                for value in values
            )
            shared_term_estimates[term] = estimate
            shared_term_spread[term] = spread_kj

        for state_name in self.state_names:
            info = state_results[state_name]
            ref_terms = info["reference_terms"]
            evb_terms = info["evb_terms"]
            evb_terms_with_shared_estimate: Dict[str, unit.Quantity] = {}
            corrected_term_deltas: Dict[str, unit.Quantity] = {}
            for term in STATE_ENERGY_TERMS:
                shared_estimate = shared_term_estimates[term]
                evb_plus_shared = evb_terms.get(term, self._zero_energy()) + shared_estimate
                evb_terms_with_shared_estimate[term] = evb_plus_shared
                corrected_term_deltas[term] = ref_terms[term] - evb_plus_shared
            info["evb_terms_with_shared_estimate"] = evb_terms_with_shared_estimate
            info["term_deltas_with_shared_estimate"] = corrected_term_deltas

        shared_term_total = sum(shared_term_estimates.values(), self._zero_energy())
        shared_term_estimate_residual = shared_total - shared_term_total
        evb_state_total = sum(
            (terms.get("total", self._zero_energy()) for terms in evb_terms_by_state.values()),
            self._zero_energy(),
        )
        evb_state_total_with_shared = evb_state_total + shared_total * len(self.state_names)
        grand_total = evb_state_total + shared_total

        debug_data = {
            "evb_breakdown": {
                "states": evb_terms_by_state,
                "shared_total": shared_total,
                "grand_total": grand_total,
            },
            "states": state_results,
            "reference_grand_total": reference_total,
            "evb_state_total": evb_state_total,
            "evb_state_total_with_shared": evb_state_total_with_shared,
            "shared_term_estimates": shared_term_estimates,
            "shared_term_spread": shared_term_spread,
            "shared_term_estimate_residual": shared_term_estimate_residual,
            "state_total_delta": reference_total - evb_state_total_with_shared,
        }
        if print_summary:
            self._print_energy_debug_summary(debug_data)
        return debug_data

    @staticmethod
    def _nonbonded_particle_params(system: openmm.System):
        """(charge_e, sigma_nm, epsilon_kJ) per particle from the first NonbondedForce."""
        nbf = [system.getForce(k) for k in range(system.getNumForces())
               if isinstance(system.getForce(k), openmm.NonbondedForce)]
        if not nbf:
            return None
        force = nbf[0]
        out = []
        for i in range(force.getNumParticles()):
            q, s, e = force.getParticleParameters(i)
            out.append((q.value_in_unit(unit.elementary_charge),
                        s.value_in_unit(unit.nanometer),
                        e.value_in_unit(unit.kilojoules_per_mole)))
        return out

    @staticmethod
    def _compare_nonbonded_consistency(named_systems, allowed_atoms, charge_tol, lj_tol,
                                       index_aligner=None):
        """List per-particle NonbondedForce mismatches across systems, excluding
        ``allowed_atoms``. Compares every system to the first.

        ``index_aligner(ref_name, other_name, ref_index) -> other_index | None`` maps a
        reference-state particle index to the corresponding particle index in another
        system. Diabatic states may legitimately order their atoms differently (e.g. a
        transferring proton sits in a different residue per state in proton-transfer EVB),
        and the EVB scaffold aligns them through the atom mapping -- so the consistency
        check must compare *mapped* counterparts, not the same raw index. Defaults to the
        identity (raw index) for backward compatibility.
        """
        names = list(named_systems)
        params = {n: EVBSystem._nonbonded_particle_params(named_systems[n]) for n in names}
        counts = {n: (len(p) if p is not None else None) for n, p in params.items()}
        ref = names[0]
        natoms = counts[ref]
        mismatches = []
        for n in names[1:]:
            if counts[n] != natoms:
                mismatches.append({"kind": "count", "a": ref, "b": n,
                                   "a_count": natoms, "b_count": counts[n]})
        for i in range(natoms or 0):
            if i in allowed_atoms:
                continue
            base = params[ref][i]
            for n in names[1:]:
                if counts[n] != natoms:
                    continue
                j = index_aligner(ref, n, i) if index_aligner is not None else i
                if j is None:
                    mismatches.append({"kind": "unmapped", "atom": i, "a": ref, "b": n})
                    continue
                p = params[n][j]
                if (abs(p[0] - base[0]) > charge_tol or abs(p[1] - base[1]) > lj_tol
                        or abs(p[2] - base[2]) > lj_tol):
                    mismatches.append({"kind": "particle", "atom": i, "other_atom": j,
                                       "a": ref, "b": n, "a_params": base, "b_params": p})

        return mismatches

    def check_forcefield_consistency(self, reference_system=None, reactive_atom_indices=None,
                                     *, charge_tol=1e-6, lj_tol=1e-6, system_kwargs=None,
                                     raise_on_mismatch=True, print_summary=True):
        """Verify the EVB build's non-reactive force field is identical across diabatic
        states (and, if given, vs an external reference ``openmm.System``), differing
        ONLY at the reactive-region atoms.

        The MD->EVB snapshot handoff is only valid if protein/solvent/ion parameters are
        the same in the equilibration force field and the EVB build. This builds the pure
        per-state Systems (:meth:`build_state_systems`) and compares per-particle
        NonbondedForce parameters (charge, sigma, epsilon) for every NON-reactive atom:
        across states, and optionally against a reference System (e.g. the MD prmtop's
        ``AmberPrmtopFile(...).createSystem()``).

        Parameters
        ----------
        reference_system : openmm.System or None
            External classical System whose non-reactive particles must match the EVB
            source-state build (catches water-model / FF-version / charge divergence
            between the equilibration FF and the EVB FF).
        reactive_atom_indices : Iterable[int] or None
            Atoms allowed to differ. Defaults to every atom in residues named in
            ``reactive_residue_names`` plus any forming/breaking-bond atoms.
        charge_tol, lj_tol : float
            Tolerances (e / nm / kJ-mol) for declaring a parameter unchanged.
        raise_on_mismatch : bool
            Raise ``ValueError`` if any non-reactive parameter differs.

        Returns
        -------
        dict with ``n_atoms``, ``n_reactive_atoms``, ``state_mismatches``,
        ``reference_mismatches`` and ``consistent``.
        """
        state_systems = dict(self.build_state_systems(**(system_kwargs or {})))
        if not state_systems:
            raise ValueError("No per-state systems available; build the EVB system first.")
        names = list(state_systems)

        if reactive_atom_indices is None:
            allowed = set()
            src = self._evb_system_source_state or names[0]
            # Atoms in reactive residues OR per-state (custom OpenFF) residues legitimately
            # differ between diabatic states (reactive bonds / per-state charges).
            allowed_resnames = set(self._reactive_residue_names) | {
                name for name, molecule in self._openff_reference_molecules.items()
                if isinstance(molecule, Mapping)
            }
            allowed_norm = {self._normalize_residue_name(n) for n in allowed_resnames}
            for atom in self._modellers[src].topology.atoms():
                if (atom.residue.name in allowed_resnames
                        or self._normalize_residue_name(atom.residue.name) in allowed_norm):
                    allowed.add(atom.index)
            try:
                for pair in self.get_forming_breaking_bond_pairs(prefer_index_pairs=True):
                    for a in pair:
                        if isinstance(a, int):
                            allowed.add(a)
            except Exception:
                pass
        else:
            allowed = {int(i) for i in reactive_atom_indices}

        def _index_aligner(ref_name, other_name, ref_index):
            # Diabatic states may order atoms differently (e.g. a transferring proton in a
            # different residue per state); align via the EVB atom-mapping's canonical
            # labels. An external reference system shares the ref state's topology order.
            if other_name not in self.state_names:
                return ref_index
            label = self._canonical_atom_label(ref_name, int(ref_index))
            return self._state_canonical_label_index_map(other_name).get(label)

        state_mismatches = self._compare_nonbonded_consistency(
            state_systems, allowed, charge_tol, lj_tol, index_aligner=_index_aligner)

        reference_mismatches = []
        if reference_system is not None:
            ref_name = self._evb_system_source_state or names[0]
            reference_mismatches = self._compare_nonbonded_consistency(
                {ref_name: state_systems[ref_name], "__reference__": reference_system},
                allowed, charge_tol, lj_tol, index_aligner=_index_aligner)

        first_params = self._nonbonded_particle_params(state_systems[names[0]])
        report = {
            "n_atoms": len(first_params) if first_params else 0,
            "n_reactive_atoms": len(allowed),
            "state_mismatches": state_mismatches,
            "reference_mismatches": reference_mismatches,
            "consistent": not state_mismatches and not reference_mismatches,
        }
        if print_summary:
            print(f"[evb.ffcheck] {report['n_atoms']} atoms, {len(allowed)} reactive; "
                  f"non-reactive across-state mismatches={len(state_mismatches)}, "
                  f"vs-reference mismatches={len(reference_mismatches)} -> "
                  f"{'CONSISTENT' if report['consistent'] else 'DIVERGENT'}")
        if raise_on_mismatch and not report["consistent"]:
            first = (state_mismatches + reference_mismatches)[0]
            raise ValueError(
                "EVB force-field consistency check failed: non-reactive parameters differ "
                f"({len(state_mismatches)} across-state, {len(reference_mismatches)} vs "
                f"reference). First mismatch: {first}")
        return report

    def __repr__(self) -> str:
        entries = ", ".join(
            f"{name}={state.path!s} (n_atoms={state.n_atoms})" for name, state in self._states.items()
        )
        return f"EVBSystem({entries})"

    def _register_force_to_group(self, force: openmm.Force, *, term: str, state: str | None) -> None:
        group = self._get_or_create_force_group(term=term, state=state)
        force.setForceGroup(group)

    def _get_or_create_force_group(self, *, term: str, state: str | None) -> int:
        group, next_group_id = get_or_create_force_group(
            force_group_ids=self._force_group_ids,
            force_group_definitions=self._force_group_definitions,
            next_force_group_id=self._next_force_group_id,
            term=term,
            state=state,
        )
        self._next_force_group_id = next_group_id
        return group

    def _collect_group_energies(self, context: openmm.Context) -> Dict[int, unit.Quantity]:
        return collect_group_energies(
            context=context,
            force_group_definitions=self._force_group_definitions,
        )

    def _adiabatic_force_group_masks(self) -> Tuple[int, int, int]:
        """Return (universal, state1, state2) force-group bitmasks for adiabatic mode.

        Used by direct adiabatic ground-state dynamics to evaluate per-state
        diabatic energies and forces in a single force pass via OpenMM's
        per-force-group accessors. The three masks are pairwise disjoint and
        together cover every force currently registered on the EVB scaffold.
        """

        return adiabatic_force_group_masks(
            state_names=self.state_names,
            force_group_definitions=self._force_group_definitions,
        )

    @staticmethod
    def _context_parameter_exists(context: openmm.Context, parameter_name: str) -> bool:
        return context_parameter_exists(context, parameter_name)

    @staticmethod
    def _context_get_parameter_or_default(
        context: openmm.Context,
        parameter_name: str,
        *,
        default: float,
    ) -> float:
        return context_get_parameter_or_default(context, parameter_name, default=default)

    def _context_lambda_snapshot(self, context: openmm.Context) -> Dict[str, float]:
        return context_lambda_snapshot(
            context=context,
            state_names=self.state_names,
            state_lambda_parameter_names=self._state_lambda_parameter_names,
            context_get_parameter_or_default=lambda ctx, name, default: self._context_get_parameter_or_default(
                ctx,
                name,
                default=default,
            ),
        )

    def _apply_context_lambdas(self, context: openmm.Context, mapping: Mapping[str, float]) -> None:
        apply_context_lambdas(
            context=context,
            state_names=self.state_names,
            state_lambda_parameter_names=self._state_lambda_parameter_names,
            mapping=mapping,
            context_parameter_exists=self._context_parameter_exists,
        )
        self._sync_pme_g_parameters(context)

    def _sync_pme_g_parameters(self, context: openmm.Context) -> None:
        """Drive each periodic-PME charge-scale global to ``sqrt(lambda_state)``.

        The per-state PME ``NonbondedForce`` scales its charges by ``g_state`` and
        therefore its Ewald energy by ``g_state^2``; setting ``g_state =
        sqrt(lambda_state)`` makes that energy ``lambda_state * E_elec(q_state)``.
        Reading lambda back from the context (rather than trusting the caller's
        mapping) keeps the two parameters in lock-step through every lambda-setting
        path, including the diabatic force-group sampler that drives lambda to 0/1.
        """
        if not self._state_pme_g_param_names:
            return
        for state_name, g_param in self._state_pme_g_param_names.items():
            if not self._context_parameter_exists(context, g_param):
                continue
            lambda_name = self._state_lambda_parameter_names.get(state_name) or f"lambda_{state_name}"
            lam = self._context_get_parameter_or_default(context, lambda_name, default=0.0)
            context.setParameter(g_param, math.sqrt(max(float(lam), 0.0)))

    def _apply_indicator_lambdas(self, context: openmm.Context, *, active_state: str) -> None:
        apply_indicator_lambdas(
            context=context,
            active_state=active_state,
            state_names=self.state_names,
            apply_context_lambdas=self._apply_context_lambdas,
        )

    def _assign_reference_force_groups(
        self, system: openmm.System
    ) -> Tuple[Dict[int, Dict[str, object]], Tuple[int, ...]]:
        return assign_reference_force_groups(
            system=system,
            categorize_reference_force=self._categorize_reference_force,
        )

    @staticmethod
    def _restore_reference_force_groups(system: openmm.System, previous_groups: Tuple[int, ...]) -> None:
        restore_reference_force_groups(system, previous_groups)

    @staticmethod
    def _categorize_reference_force(force: openmm.Force) -> str:
        return categorize_reference_force(force)

    def _remap_positions_for_state(
        self,
        *,
        source_state: str,
        target_state: str,
        source_positions: Tuple[openmm.Vec3, ...],
    ) -> Tuple[openmm.Vec3, ...]:
        source_atoms = self._state_atom_references(source_state)
        target_atoms = self._state_atom_references(target_state)
        if len(source_positions) != len(source_atoms):
            raise ValueError(
                f"Source state '{source_state}' provides {len(source_positions)} positions but "
                f"{len(source_atoms)} atoms are registered."
            )

        label_to_position: Dict[AtomReference, openmm.Vec3] = {}
        for idx in range(len(source_atoms)):
            label = self._canonical_atom_label(source_state, idx)
            label_to_position[label] = source_positions[idx]

        remapped: list[openmm.Vec3] = []
        for idx in range(len(target_atoms)):
            label = self._canonical_atom_label(target_state, idx)
            try:
                remapped.append(label_to_position[label])
            except KeyError as exc:
                raise ValueError(
                    f"Atom reference {label!r} exists in state '{target_state}' but is unmapped from '{source_state}'."
                ) from exc
        return tuple(remapped)

    @staticmethod
    def _coerce_vec3_tuple(
        values: unit.Quantity | Sequence[openmm.Vec3] | Sequence[Sequence[float]],
        *,
        quantity_unit: unit.Unit,
        expected_count: int | None = None,
        label: str = "values",
    ) -> Tuple[openmm.Vec3, ...]:
        def _coerce_component(component: object) -> float:
            if isinstance(component, unit.Quantity):
                return float(component.value_in_unit(quantity_unit))
            return float(component)

        raw_values: Sequence[object]
        if isinstance(values, unit.Quantity):
            try:
                raw_values = tuple(values.value_in_unit(quantity_unit))
            except (AttributeError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"'{label}' must be convertible to units of {quantity_unit}."
                ) from exc
        else:
            raw_values = tuple(values)

        parsed: list[openmm.Vec3] = []
        for idx, vector in enumerate(raw_values):
            try:
                if isinstance(vector, openmm.Vec3):
                    x, y, z = float(vector.x), float(vector.y), float(vector.z)
                elif isinstance(vector, unit.Quantity):
                    vec = vector.value_in_unit(quantity_unit)
                    if isinstance(vec, openmm.Vec3):
                        x, y, z = float(vec.x), float(vec.y), float(vec.z)
                    else:
                        x = _coerce_component(vec[0])  # type: ignore[index]
                        y = _coerce_component(vec[1])  # type: ignore[index]
                        z = _coerce_component(vec[2])  # type: ignore[index]
                else:
                    x = _coerce_component(vector[0])  # type: ignore[index]
                    y = _coerce_component(vector[1])  # type: ignore[index]
                    z = _coerce_component(vector[2])  # type: ignore[index]
            except (AttributeError, IndexError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"'{label}' entry at index {idx} is not a valid 3-vector."
                ) from exc
            parsed.append(openmm.Vec3(x, y, z))

        if expected_count is not None and len(parsed) != int(expected_count):
            raise ValueError(
                f"'{label}' must contain {int(expected_count)} vectors, received {len(parsed)}."
            )
        return tuple(parsed)

    @staticmethod
    def _vec3s_to_quantity(
        positions: Tuple[openmm.Vec3, ...], *, quantity_unit: unit.Unit = unit.nanometer
    ) -> unit.Quantity:
        coords = [(float(vec.x), float(vec.y), float(vec.z)) for vec in positions]
        return unit.Quantity(coords, quantity_unit)

    def _reduce_reference_terms(
        self,
        context: openmm.Context,
        group_metadata: Mapping[int, Mapping[str, object]],
        group_energies: Mapping[int, unit.Quantity],
    ) -> Tuple[Dict[str, unit.Quantity], Dict[str, unit.Quantity]]:
        return reduce_reference_terms(
            context=context,
            group_metadata=group_metadata,
            group_energies=group_energies,
            split_nonbonded_energy=self._split_nonbonded_energy,
        )

    def _split_nonbonded_energy(
        self,
        *,
        force: openmm.NonbondedForce,
        context: openmm.Context,
        group: int,
    ) -> Tuple[unit.Quantity, unit.Quantity]:
        return split_nonbonded_energy(force=force, context=context, group=group)

    def _print_energy_debug_summary(self, data: Mapping[str, object]) -> None:
        print_energy_debug_summary(data, state_names=self.state_names)

    @staticmethod
    def _zero_energy() -> unit.Quantity:
        return 0.0 * unit.kilojoules_per_mole

    @staticmethod
    def _validate_state_name(name: str) -> None:
        validate_state_name(name)

    @staticmethod
    def _next_auto_name(existing: set[str]) -> str:
        return next_auto_name(existing)

    #: XML files behind the ``"amber14"`` preset.
    AMBER14_FORCE_FIELD_FILES: Tuple[str, ...] = ("amber14-all.xml", "amber14/tip3p.xml")

    #: Absolute paths to the OPLS-AA/M XMLs shipped with the library, behind ``"opls-aam"``.
    #: opls-aam.xml carries a <Script> that installs the geometric-combining CustomNonbondedForce;
    #: opls-aam-patches.xml adds caps/protonation/disulfide + ions (ACE/NME/GLUP/CYX/Ca).
    _OPLS_DATA_DIR = Path(__file__).resolve().parent / "data" / "opls_aam"
    OPLSAAM_FORCE_FIELD_FILES: Tuple[str, ...] = (
        str(_OPLS_DATA_DIR / "opls-aam.xml"),
        str(_OPLS_DATA_DIR / "tip3p.xml"),
        str(_OPLS_DATA_DIR / "opls-aam-patches.xml"),
    )
    #: OPLS-consistent water templates (coulomb14scale=0.5, LJ routed through the geometric
    #: CustomNonbondedForce) keyed by solvation water model. Selected in ``_load_force_field``.
    _OPLS_WATER_MODEL_FILES: "Dict[str, str]" = {
        "tip3p": str(_OPLS_DATA_DIR / "tip3p.xml"),
        "spce": str(_OPLS_DATA_DIR / "spce.xml"),
    }
    #: Preset names that select OPLS-AA/M (and imply geometric LJ combining).
    OPLS_PRESET_NAMES: frozenset = frozenset({"opls-aam", "opls", "oplsaam", "opls-aa/m"})

    @classmethod
    def spec_is_opls(cls, force_field) -> bool:
        """True if the force-field spec selects OPLS-AA/M (=> geometric LJ combining).

        Detects both the preset name ("opls-aam") and the resolved ``opls-aam.xml`` path, so
        auto-geometric survives a spec given as explicit file paths (e.g. a preset expanded by a
        caller, or a prepared bundle that stored the resolved files rather than the preset name).
        """
        if isinstance(force_field, (str, os.PathLike)):
            items = [force_field]
        else:
            items = list(force_field)
        for i in items:
            s = str(i).lower()
            if s in cls.OPLS_PRESET_NAMES or os.path.basename(s) == "opls-aam.xml":
                return True
        return False

    def _opls_host_coulomb14scale(self) -> float | None:
        """Coulomb 1-4 scale of the host force field, but only for an OPLS spec.

        Returned so the SMIRNOFF ligand template can be reconciled to the OPLS 1-4 scale
        (otherwise OpenMM rejects merging the two NonbondedForce tags). ``None`` for non-OPLS
        force fields leaves the standard SMIRNOFF template untouched."""
        if not self.spec_is_opls(self.force_field_name):
            return None
        for gen in getattr(self.force_field, "_forces", ()):
            scale = getattr(gen, "coulomb14scale", None)
            if scale is not None:
                return float(scale)
        return 0.5   # OPLS-AA/M convention if the generator is not introspectable

    def _apply_opls_ligand_overrides(self, system, topology, state_name) -> None:
        """Overwrite each ffld-provided ligand's LJ + bonded + 1-4 with OPLS-AA/M, in place.

        ``self._ligand_ffld_files`` maps residue name -> a single ffld path (all states) or a
        ``{state_name: path}`` mapping. Charges are left untouched (``set_charges=False``) so the
        ligand_charges / SDF machinery stays the single charge axis.

        A per-state mapping that omits a state declares that this ligand species is ABSENT from
        that diabat -- the case of a cross-residue proton transfer, where the substrate is a
        different molecule in each state (keto ligand in the reactant, enediolate in the product,
        each under its own residue name). That omission is honoured, but only after checking the
        residue really is absent from the state's topology: a residue that IS present with no ffld
        entry would otherwise be left silently on its SMIRNOFF parameters while the rest of the
        system is OPLS."""
        if not self._ligand_ffld_files:
            return
        from .parameterizer.opls import (
            apply_opls_ligand, ligand_start_from_topology, parse_ffld,
        )
        for res_name, spec in self._ligand_ffld_files.items():
            if isinstance(spec, Mapping):
                if state_name not in spec:
                    present = any(residue.name.strip() == res_name.strip()
                                  for residue in topology.residues())
                    if present:
                        raise ValueError(
                            f"ligand_ffld_files for residue '{res_name}' is a per-state mapping "
                            f"with no entry for state '{state_name}', but that state's topology "
                            f"DOES contain '{res_name}'. It would keep its SMIRNOFF parameters "
                            f"while the rest of the system is OPLS; provide an ffld for this "
                            f"state (or a single path for all states)."
                        )
                    continue                    # species genuinely absent from this diabat
                path = spec[state_name]
            else:
                path = spec
            ffld = parse_ffld(str(path))
            start = ligand_start_from_topology(topology, res_name)
            # The ffld overwrites the ligand block [start, start+len) BY INDEX, so the topology's
            # ligand-atom order must match the ffld/SDF order. Verify: (1) atom count matches, and
            # (2) every ffld bond apply_opls_ligand could not find in the HarmonicBondForce is
            # actually a System constraint (a constrained X-H) -- otherwise the orders disagree and
            # the ligand would be silently mis-parameterized.
            residue = next(a.residue for a in topology.atoms() if a.index == start)
            n_res = sum(1 for _ in residue.atoms())
            if n_res != len(ffld.atoms):
                raise ValueError(
                    f"ligand_ffld_files: residue '{res_name}' has {n_res} atoms in state "
                    f"'{state_name}' but its ffld '{path}' has {len(ffld.atoms)} -- atom counts "
                    f"must match (same molecule, same order)."
                )
            report = apply_opls_ligand(system, ffld, start, set_charges=False)
            constrained = {
                frozenset(system.getConstraintParameters(c)[:2])
                for c in range(system.getNumConstraints())
            }
            unmatched = [(i, j) for (i, j) in report.bonds_missing
                         if frozenset((start + i, start + j)) not in constrained]
            if unmatched:
                raise ValueError(
                    f"ligand_ffld_files: {len(unmatched)} ffld bond(s) for residue '{res_name}' "
                    f"(state '{state_name}') match no System bond or constraint, e.g. {unmatched[:3]} "
                    f"(ligand-local indices) -- the ffld atom order does not match the topology/SDF "
                    f"order for this ligand. Re-export the ffld from the same structure/order."
                )

    @classmethod
    def resolve_force_field_files(cls, force_field) -> Tuple[str, ...]:
        """Expand a force-field spec into an ordered tuple of OpenMM XML identifiers.

        Accepts the ``"amber14"`` and ``"opls-aam"`` presets, a single XML path or built-in
        name, or a sequence mixing the two -- so extra XMLs layer on top of a preset::

            force_field=("amber14", "my_custom_residues.xml")
            force_field=("opls-aam", "my_reactive_template.xml")

        Later files override earlier ones, per OpenMM's ``ForceField`` semantics.
        """
        if isinstance(force_field, (str, os.PathLike)):
            items = [force_field]
        else:
            items = list(force_field)
        files: list[str] = []
        for item in items:
            text = str(item)
            if text.lower() == "amber14":
                files.extend(cls.AMBER14_FORCE_FIELD_FILES)
            elif text.lower() in cls.OPLS_PRESET_NAMES:
                files.extend(cls.OPLSAAM_FORCE_FIELD_FILES)
            else:
                files.append(text)
        return tuple(files)

    def _apply_water_model_to_files(self, files: Tuple[str, ...]) -> Tuple[str, ...]:
        """Swap the OPLS-shipped ``tip3p.xml`` for the water template matching the most
        recent solvation call (``self._solvation_water_model``).

        Only the library's own OPLS water template is swapped, and only away from the
        default TIP3P, so gas-phase and non-OPLS specs are untouched. Raises if a
        non-default water model is requested but no OPLS water template is present to
        replace (e.g. an SPC/E solvation under a force field that did not ship one)."""
        model = getattr(self, "_solvation_water_model", "tip3p")
        if model == "tip3p":
            return files
        replacement = self._OPLS_WATER_MODEL_FILES.get(model)
        if replacement is None:                      # unknown model -> defer to the name check
            return files
        opls_tip3p = os.path.abspath(self._OPLS_WATER_MODEL_FILES["tip3p"])
        swapped = []
        replaced = False
        for f in files:
            if os.path.abspath(str(f)) == opls_tip3p:
                swapped.append(replacement)
                replaced = True
            else:
                swapped.append(f)
        if not replaced:
            raise ValueError(
                f"Water model '{model}' was requested, but no OPLS water template "
                f"(tip3p.xml) is present in the force-field spec to replace. The "
                f"'{model}' water template is only wired for the 'opls-aam' preset."
            )
        return tuple(swapped)

    def _load_force_field(self, force_field) -> app.ForceField:
        """Build the ``app.ForceField`` backing every per-state ``createSystem``.

        ``force_field`` may be:

        * ``"amber14"`` -- the built-in preset (default, unchanged);
        * a path or built-in name of a single OpenMM force-field XML;
        * a sequence of those, e.g. ``("amber14", "tetrahedral_intermediate.xml")``,
          which is how a reactive state whose topology no standard residue template
          covers (a covalent adduct, a cross-residue forming bond) gets parameterized.

        An already-constructed ``app.ForceField`` is deliberately **not** accepted. The
        spec is re-loaded to build a *fresh* ``ForceField`` per state (see
        ``_clone_force_field_with_generators``), so a shared instance would leak the
        template generators registered for one state into the others; and it could not
        be recorded in a prepared bundle, which stores the spec as JSON.
        """
        if isinstance(force_field, app.ForceField):
            raise TypeError(
                "force_field must be a spec, not a built app.ForceField. Each diabatic "
                "state is given its own freshly-loaded ForceField, so a shared instance "
                "would leak one state's template generators into the others, and it "
                "cannot be serialized into a prepared bundle. Pass the XML files "
                "instead, e.g. force_field=('amber14', 'my_residues.xml')."
            )
        files = self.resolve_force_field_files(force_field)
        files = self._apply_water_model_to_files(files)
        if not files:
            raise ValueError(
                "No force-field files were given. Pass 'amber14', an OpenMM force-field "
                "XML path or name, or a sequence of them."
            )
        try:
            return app.ForceField(*files)
        except Exception as exc:
            raise ValueError(
                f"Could not load force field from {list(files)}: {exc}. Pass 'amber14', "
                "an OpenMM force-field XML path or built-in name, or a sequence of them "
                "(e.g. ('amber14', 'my_residues.xml'))."
            ) from exc

    @staticmethod
    def _load_openff_parameterizer() -> Tuple[object, object]:
        """
        Import OpenFF-based parameterization utilities lazily.

        This keeps EVBSystem importable in environments that run only pre-parameterized
        systems and therefore do not require the OpenFF/RDKit stack.
        """

        try:
            from .parameterizer.openff import instantiate_template_generator, parameterize_with_openff
        except ImportError as exc:  # pragma: no cover - optional dependency path
            raise ImportError(
                "OpenFF parameterization utilities are unavailable. Install optional "
                "dependencies required by 'evbopenmm.parameterizer.openff' "
                "(openff-toolkit, openmmforcefields, rdkit, ambertools)."
            ) from exc
        return parameterize_with_openff, instantiate_template_generator

    def _register_template_generator_specs(self, specs: Iterable[object]) -> None:
        if not specs:
            return
        existing = list(self._template_generator_specs)
        existing_keys = {self._generator_spec_key(spec) for spec in existing}
        updated = False
        for spec in specs:
            key = self._generator_spec_key(spec)
            if key in existing_keys:
                continue
            existing.append(spec)
            existing_keys.add(key)
            updated = True
        if updated:
            self._template_generator_specs = tuple(existing)

    @staticmethod
    def _generator_spec_key(spec: object) -> Tuple[object, ...]:
        if not all(hasattr(spec, attr) for attr in ("kind", "kwargs", "molecules")):
            raise TypeError(
                "Template generator specification must define 'kind', 'kwargs', and 'molecules' attributes."
            )
        kwargs_key = tuple(sorted((str(key), str(value)) for key, value in spec.kwargs.items()))
        molecule_key: list[Tuple[str, str]] = []
        for mol in spec.molecules:
            residue_name = mol.properties.get("residue_name", mol.name)
            try:
                smiles = mol.to_smiles(isomeric=True)
            except (AttributeError, RuntimeError, ValueError):  # pragma: no cover
                smiles = mol.to_json()
            molecule_key.append((residue_name, smiles))
        return (spec.kind, kwargs_key, tuple(molecule_key))

    def _clone_force_field_with_generators(self) -> app.ForceField:
        """
        Construct a fresh ForceField instance and re-register any template generators.

        A new ForceField is required per state to avoid cross-state mutations. However,
        residues parameterized via OpenFF (or other template generators) are registered on
        ``self.force_field``; this helper mirrors those generators onto the fresh copy.
        """

        ff = self._load_force_field(self.force_field_name)
        for generator in self._base_template_generators:
            ff.registerTemplateGenerator(generator)
        if self._template_generator_specs:
            _, instantiate_template_generator = self._load_openff_parameterizer()
        for spec in self._template_generator_specs:
            generator_obj = instantiate_template_generator(spec)
            ff.registerTemplateGenerator(generator_obj.generator)
        return ff

    @staticmethod
    def _deserialize_state_system_xmls(state_system_xmls: Mapping[str, PathLike]) -> Dict[str, openmm.System]:
        systems: Dict[str, openmm.System] = {}
        for state_name, raw_path in state_system_xmls.items():
            path = Path(raw_path).expanduser()
            if not path.is_file():
                raise FileNotFoundError(f"Serialized state-system XML not found for '{state_name}': {path}")
            try:
                system = openmm.XmlSerializer.deserialize(path.read_text(encoding="utf-8"))
            except (openmm.OpenMMException, ValueError) as exc:
                raise ValueError(f"Failed to deserialize state-system XML for '{state_name}': {path}") from exc
            if not isinstance(system, openmm.System):
                raise TypeError(
                    f"Serialized state-system XML for '{state_name}' did not produce an OpenMM System instance."
                )
            systems[str(state_name)] = system
        return systems

    def _load_prebuilt_state_systems(self, systems: Mapping[str, openmm.System]) -> None:
        expected_states = set(self.state_names)
        provided_states = set(systems)
        missing_states = sorted(expected_states - provided_states)
        extra_states = sorted(provided_states - expected_states)
        if missing_states or extra_states:
            details: list[str] = []
            if missing_states:
                details.append("missing -> " + ", ".join(missing_states))
            if extra_states:
                details.append("extra -> " + ", ".join(extra_states))
            raise ValueError("Prepared state systems do not match the configured EVB states: " + "; ".join(details))

        loaded: Dict[str, openmm.System] = {}
        for state_name in self.state_names:
            system = systems[state_name]
            self._validate_prebuilt_state_system(state_name, system)
            loaded[state_name] = system

        self._state_systems = loaded
        self._state_system_signature = tuple()
        self._state_systems_are_prebuilt = True
        self.bonded_force_comparisons = None
        self._invalidate_partition_caches()
        self._prune_cached_state_systems()
        self._state_systems_cache = MappingProxyType(self._state_systems)

    def _validate_prebuilt_state_system(self, state_name: str, system: openmm.System) -> None:
        expected_particles = self._states[state_name].n_atoms
        actual_particles = int(system.getNumParticles())
        if actual_particles != expected_particles:
            raise ValueError(
                f"Prepared system for state '{state_name}' contains {actual_particles} particles, "
                f"expected {expected_particles} from the prepared topology."
            )

        # An OPLS system straight out of createSystem keeps all LJ in the opls-aam.xml Script's
        # geometric CustomNonbondedForce and ZEROES the NonbondedForce LJ. EVBSystem reads per-atom
        # sigma/epsilon from the NonbondedForce, so an unflattened prebuilt OPLS system would give a
        # silent zero-vdW result. Catch it and point at the fix.
        from .forces.opls import _is_geometric_lj_force
        cnb_geom = any(_is_geometric_lj_force(system.getForce(i))
                       for i in range(system.getNumForces()))
        if cnb_geom:
            nb = next((system.getForce(i) for i in range(system.getNumForces())
                       if isinstance(system.getForce(i), openmm.NonbondedForce)), None)
            nb_lj_zero = nb is not None and all(
                nb.getParticleParameters(i)[2].value_in_unit(unit.kilojoule_per_mole) == 0.0
                for i in range(nb.getNumParticles()))
            if nb_lj_zero:
                raise ValueError(
                    f"Prebuilt state '{state_name}' looks like an unflattened OPLS system: its LJ is in "
                    "a geometric CustomNonbondedForce and the NonbondedForce LJ is zeroed, so EVBSystem "
                    "would read zero van der Waals. Call "
                    "evbopenmm.forces.flatten_opls_lj_to_nonbonded(system) before passing it as a "
                    "prebuilt_state_system (or let EVBSystem build the state itself via "
                    "force_field='opls-aam')."
                )

        expected_masses = self._state_particle_masses(state_name)
        for atom_index, expected_mass in enumerate(expected_masses):
            actual_mass = system.getParticleMass(atom_index)
            if not self._particle_masses_match(actual_mass, expected_mass):
                raise ValueError(
                    f"Prepared system for state '{state_name}' has a particle mass mismatch at atom index "
                    f"{atom_index}: system={self._mass_value_dalton(actual_mass):.6f} Da, "
                    f"topology={self._mass_value_dalton(expected_mass):.6f} Da."
                )

        try:
            constraint_map = self._state_constraint_map(state_name=state_name, system=system)
        except (KeyError, ValueError) as exc:
            raise ValueError(
                f"Prepared system for state '{state_name}' defines constraints incompatible with the prepared "
                "topology atom ordering."
            ) from exc
        if len(constraint_map) != int(system.getNumConstraints()):
            raise ValueError(
                f"Prepared system for state '{state_name}' contains duplicate constraint pairs, which makes "
                "the prepared topology/system bundle ambiguous."
            )

    @staticmethod
    def _mass_value_dalton(value: unit.Quantity | float) -> float:
        if isinstance(value, Real):
            return float(value)
        return float(value.value_in_unit(unit.dalton))

    @classmethod
    def _particle_masses_match(
        cls,
        actual_mass: unit.Quantity | float,
        expected_mass: unit.Quantity | float,
        *,
        tolerance_dalton: float = 1e-2,
    ) -> bool:
        return abs(cls._mass_value_dalton(actual_mass) - cls._mass_value_dalton(expected_mass)) <= tolerance_dalton

    @staticmethod
    def _is_trivial_nonbonded_force(force: openmm.NonbondedForce, *, tol: float = 1e-12) -> bool:
        for idx in range(force.getNumParticles()):
            charge, sigma, epsilon = force.getParticleParameters(idx)
            if (
                abs(charge.value_in_unit(unit.elementary_charge)) > tol
                or abs(epsilon.value_in_unit(unit.kilojoules_per_mole)) > tol
            ):
                return False
        for idx in range(force.getNumExceptions()):
            _, _, charge_prod, sigma, epsilon = force.getExceptionParameters(idx)
            if (
                abs(charge_prod.value_in_unit(unit.elementary_charge**2)) > tol
                or abs(epsilon.value_in_unit(unit.kilojoules_per_mole)) > tol
            ):
                return False
        return True

    def _prune_trivial_nonbonded_forces(self, system: openmm.System) -> None:
        remove_indices: list[int] = []
        for idx, force in enumerate(system.getForces()):
            if isinstance(force, openmm.NonbondedForce) and self._is_trivial_nonbonded_force(force):
                remove_indices.append(idx)
        for offset, force_index in enumerate(remove_indices):
            system.removeForce(force_index - offset)

    def _prune_cached_state_systems(self) -> None:
        if not self._state_systems:
            return
        for system in self._state_systems.values():
            self._prune_trivial_nonbonded_forces(system)

    def _snapshot_per_state_solute_nb_params(
        self, *, solute_counts: Mapping[str, int]
    ) -> Dict[str, Dict[str, object]]:
        """Snapshot per-state NonbondedForce solute parameters (particles + exceptions).

        Captures (charge, sigma, epsilon) for particles [0..n_solute) and any
        exception whose both endpoints are solute. Intended to be replayed after
        a state-system rebuild (e.g. during :meth:`solvate_states`) so that any
        caller-applied parameter overrides survive.
        """
        snapshots: Dict[str, Dict[str, object]] = {}
        for state_name in self.state_names:
            sys_obj = self._state_systems.get(state_name)
            if sys_obj is None:
                continue
            n_solute = int(solute_counts.get(state_name, 0))
            if n_solute <= 0:
                continue
            nb: openmm.NonbondedForce | None = None
            for fi in range(sys_obj.getNumForces()):
                force = sys_obj.getForce(fi)
                if isinstance(force, openmm.NonbondedForce):
                    nb = force
                    break
            if nb is None:
                continue
            particle_params: list[tuple[float, float, float]] = []
            for idx in range(min(n_solute, nb.getNumParticles())):
                q, sigma, eps = nb.getParticleParameters(idx)
                particle_params.append((
                    float(q.value_in_unit(unit.elementary_charge)),
                    float(sigma.value_in_unit(unit.nanometer)),
                    float(eps.value_in_unit(unit.kilojoules_per_mole)),
                ))
            exception_params: list[tuple[int, int, float, float, float]] = []
            for ei in range(nb.getNumExceptions()):
                a, b, cp, sigma, eps = nb.getExceptionParameters(ei)
                if int(a) < n_solute and int(b) < n_solute:
                    exception_params.append((
                        int(a), int(b),
                        float(cp.value_in_unit(unit.elementary_charge ** 2)),
                        float(sigma.value_in_unit(unit.nanometer)),
                        float(eps.value_in_unit(unit.kilojoules_per_mole)),
                    ))
            snapshots[state_name] = {
                "n_solute": n_solute,
                "particles": tuple(particle_params),
                "exceptions": tuple(exception_params),
            }
        return snapshots

    def _restore_per_state_solute_nb_params(
        self, snapshots: Mapping[str, Mapping[str, object]]
    ) -> None:
        """Re-apply snapshotted per-state NonbondedForce solute parameters.

        Companion to :meth:`_snapshot_per_state_solute_nb_params`; replays
        the captured (charge, sigma, epsilon) onto solute particles [0..n_solute)
        and matches solute-solute exceptions by (i, j).
        """
        for state_name, snap in snapshots.items():
            sys_obj = self._state_systems.get(state_name)
            if sys_obj is None:
                continue
            nb: openmm.NonbondedForce | None = None
            for fi in range(sys_obj.getNumForces()):
                force = sys_obj.getForce(fi)
                if isinstance(force, openmm.NonbondedForce):
                    nb = force
                    break
            if nb is None:
                continue
            for idx, (q_e, sigma_nm, eps_kj) in enumerate(snap["particles"]):
                if idx >= nb.getNumParticles():
                    break
                nb.setParticleParameters(
                    idx,
                    float(q_e) * unit.elementary_charge,
                    float(sigma_nm) * unit.nanometer,
                    float(eps_kj) * unit.kilojoules_per_mole,
                )
            if not snap["exceptions"]:
                continue
            snap_exc: Dict[tuple[int, int], tuple[float, float, float]] = {
                (min(a, b), max(a, b)): (cp, sigma, eps)
                for (a, b, cp, sigma, eps) in snap["exceptions"]
            }
            for ei in range(nb.getNumExceptions()):
                a, b, cp, sigma, eps = nb.getExceptionParameters(ei)
                key = (min(int(a), int(b)), max(int(a), int(b)))
                target = snap_exc.get(key)
                if target is None:
                    continue
                new_cp, new_sigma, new_eps = target
                nb.setExceptionParameters(
                    ei, int(a), int(b),
                    float(new_cp) * unit.elementary_charge ** 2,
                    float(new_sigma) * unit.nanometer,
                    float(new_eps) * unit.kilojoules_per_mole,
                )

    def _invalidate_partition_caches(self) -> None:
        self._vdw_partition = None
        self._electrostatics_partition = None
        self._vdw_partition_cache_key = None
        self._electrostatics_partition_cache_key = None

    def _compute_missing_force_field_residues(self) -> Dict[str, Tuple[str, ...]]:
        return compute_missing_force_field_residues(
            self._states,
            normalized_force_field_residues=set(self._normalized_force_field_residues),
        )

    def _build_residue_targets(self) -> None:
        index: ResidueTargetIndex = build_residue_targets(self._states)
        self._modellers = dict(index.modellers)
        self._topology_bonds = dict(index.topology_bonds)
        self._residue_targets = {name: list(targets) for name, targets in index.residue_targets.items()}

    def _ensure_residue_name_consistency(self, *, include_bonds: bool = True) -> None:
        # Residues parameterized per-state (per-state SDFs) legitimately differ in
        # their reactive bonds between diabatic states; exempt them from the
        # bond-consistency check (atom identity is still enforced).
        per_state_residues = {
            name
            for name, molecule in self._openff_reference_molecules.items()
            if isinstance(molecule, Mapping)
        } | set(getattr(self, "_reactive_residue_names", frozenset()))
        # Only the CUSTOM (OpenFF-parameterized) / reactive residues need cross-state
        # consistency checking. Every force-field-standard residue (amino acids, water,
        # ions) is parameterized by the force field's own templates and legitimately
        # appears under a shared name with differing atoms -- terminal vs internal
        # backbones, and protonation variants that OpenMM normalizes to one name
        # (HID/HIE/HIP -> "HIS"). Exempt all of those so protein EVB systems build.
        custom_residues = {self._normalize_residue_name(n) for n in self._openff_reference_molecules}
        custom_residues |= {self._normalize_residue_name(n) for n in per_state_residues}
        skip_residues = {name for name in self._residue_targets if name not in custom_residues}
        ensure_residue_name_consistency(
            residue_targets=self._residue_targets,
            residue_signature=lambda target, include: self._residue_signature(target, include_bonds=include),
            format_signature_difference=self._format_signature_difference,
            include_bonds=include_bonds,
            skip_residues=skip_residues,
            bonds_optional_residues=per_state_residues,
        )

    def _residue_signature(
        self, target: ResidueTarget, *, include_bonds: bool
    ) -> Tuple[Tuple[Tuple[str, str | None], ...], Tuple[Tuple[str, str], ...]]:
        return residue_signature(
            target=target,
            topology_bonds=self._topology_bonds,
            include_bonds=include_bonds,
        )

    def _residue_bond_pairs(self, target: ResidueTarget) -> Tuple[CanonicalPairKey, ...]:
        return residue_bond_pairs(
            target=target,
            topology_bonds=self._topology_bonds,
            atom_reference_from_topology_atom=self._atom_reference_from_topology_atom,
            canonical_pair_key=self._canonical_pair_key,
        )

    @staticmethod
    def _format_signature_difference(
        reference: Tuple[Tuple[Tuple[str, str | None], ...], Tuple[Tuple[str, str], ...]],
        candidate: Tuple[Tuple[Tuple[str, str | None], ...], Tuple[Tuple[str, str], ...]],
    ) -> str:
        return format_signature_difference(reference, candidate)

    def parameterize_missing_residues(
        self,
        *,
        ligand_charges: Mapping[str, float],
        ligand_sdf_files: Mapping[str, object] | None = None,
        ligand_smiles: Mapping[str, str] | None = None,
        smirnoff_forcefield: str = "openff-2.2.0.offxml",
        charge_method: str = "am1bcc",
        output_directory: PathLike | None = None,
    ) -> Dict[str, object]:
        """
        Use OpenFF/SMIRNOFF to derive parameters for residues missing from the base force field.

        Parameters
        ----------
        ligand_charges :
            Mapping of residue name -> expected total charge (in units of e).
        ligand_sdf_files :
            Optional mapping of residue name -> SDF/MOL definition. Each entry can be either a path
            or a mapping with ``path``, ``index`` (for multi-entry SDFs), and ``charge_method``.
        ligand_smiles :
            Optional mapping of residue name -> SMILES string, used when an SDF is not provided.
        smirnoff_forcefield :
            OFFXML file that defines the small-molecule parameters (defaults to Sage 2.2).
        charge_method :
            Partial-charge generation method understood by the OpenFF toolkit.
        output_directory :
            Folder where intermediate OpenFF artifacts (cached molecules) will be stored.

        Returns
        -------
        Dict[str, object]
            Mapping of residue name -> OpenFF Molecule objects that were registered.

        Notes
        -----
        The atom ordering of the provided ligand definitions (SDF/MOL/SMILES) is assumed
        to match the ordering of atoms for the corresponding residues in every PDB state.
        This is required so that reference bond connectivity can be reconstructed reliably.
        """
        missing_names = sorted({res for residues in self._missing_residues.values() for res in residues})
        if not missing_names:
            return {}

        # Under an OPLS host FF every SMIRNOFF-templated ligand MUST have an OPLS ffld: the
        # ligand's 1-4 Coulomb scale is reconciled to the OPLS 0.5 (so the templates merge) and
        # apply_opls_ligand then restores full OPLS LJ/bonded/1-4. Without an ffld the ligand
        # would keep SMIRNOFF LJ under OPLS 1-4 scaling -- a silent force-field mix. Fail loudly.
        if self.spec_is_opls(self.force_field_name):
            ffld_keys = {str(k).upper() for k in self._ligand_ffld_files}
            no_ffld = [r for r in missing_names if r.upper() not in ffld_keys]
            if no_ffld:
                raise ValueError(
                    "Under an OPLS force field every SMIRNOFF-templated ligand needs an OPLS ffld "
                    "via ligand_ffld_files (a SMIRNOFF ligand in an OPLS system is a silent FF mix). "
                    f"Missing ffld for residue(s): {', '.join(sorted(no_ffld))}."
                )

        target_mapping: Dict[str, list[ResidueTarget]] = {}
        normalized_missing = {name.upper() for name in missing_names}
        for normalized_name, contexts in self._residue_targets.items():
            if normalized_name in normalized_missing:
                target_mapping[normalized_name] = contexts

        output_path = None if output_directory is None else Path(output_directory)

        parameterize_with_openff, _ = self._load_openff_parameterizer()
        molecules, generator_specs = parameterize_with_openff(
            self.force_field,
            target_mapping,
            residue_names=missing_names,
            ligand_charges=ligand_charges,
            ligand_sdf_files=ligand_sdf_files,
            ligand_smiles=ligand_smiles,
            smirnoff_forcefield=smirnoff_forcefield,
            charge_method=charge_method,
            output_directory=output_path,
            coulomb14scale=self._opls_host_coulomb14scale(),
        )
        self._register_template_generator_specs(generator_specs)

        updated_states: Set[str] = set()
        for residue_name, molecule in molecules.items():
            display = residue_name.strip() or residue_name
            normalized = self._normalize_residue_name(display)
            self._force_field_residue_names.add(display)
            self._normalized_force_field_residues.add(normalized)
            self._openff_reference_molecules[normalized] = molecule
            if isinstance(molecule, Mapping):
                # Per-state ligand: apply each diabatic state's own molecular
                # graph (bonds) to that state's residue only, so the reactant and
                # product topologies differ as the chemistry requires.
                for state_name, state_molecule in molecule.items():
                    modified = self._apply_reference_bonds_from_openff(
                        normalized, state_molecule, only_state=state_name
                    )
                    updated_states.update(modified)
            else:
                modified = self._apply_reference_bonds_from_openff(normalized, molecule)
                updated_states.update(modified)

        self._missing_residues = self._compute_missing_force_field_residues()
        if self._missing_residues:
            details = "; ".join(
                f"{state}: {', '.join(residues)}" for state, residues in self._missing_residues.items()
            )
            raise ValueError(
                "Force field coverage remains incomplete after OpenFF parameterization: " + details
            )
        # Force-field templates/topology bond graphs changed; rebuild state systems on next request.
        self._state_systems = {}
        self._state_systems_cache = None
        self._state_system_signature = None
        self._state_systems_are_prebuilt = False
        self.bonded_force_comparisons = None
        self._invalidate_partition_caches()
        return molecules

    def _apply_reference_bonds_from_openff(
        self, normalized_residue_name: str, molecule: object, only_state: str | None = None
    ) -> Set[str]:
        """
        Ensure that residues associated with ``normalized_residue_name`` have all bonds present.

        This uses the bonding pattern from an OpenFF Molecule as the reference. The residue atoms
        in each EVB state must appear in the same order as the atoms in the reference molecule.

        ``only_state`` restricts the application to a single diabatic state, used when a residue
        has per-state molecular graphs (different bonds in reactant vs product).

        Returns
        -------
        Set[str]
            Names of EVB states whose modeller topologies were modified.
        """

        contexts = self._residue_targets.get(normalized_residue_name)
        if not contexts:
            return set()
        if only_state is not None:
            contexts = [ctx for ctx in contexts if ctx.state_name == only_state]
            if not contexts:
                return set()

        n_reference_atoms = molecule.n_atoms
        reference_bonds = tuple((bond.atom1_index, bond.atom2_index) for bond in molecule.bonds)
        modified_states: Set[str] = set()

        for target in contexts:
            residue_atoms = tuple(target.residue.atoms())
            if len(residue_atoms) != n_reference_atoms:
                display_name = target.residue.name.strip() or normalized_residue_name
                raise ValueError(
                    f"Residue '{display_name}' in state '{target.state_name}' contains {len(residue_atoms)} atoms, "
                    f"but the OpenFF reference expects {n_reference_atoms}. Ensure the ligand definition uses the "
                    "same atom ordering as the PDB residues."
                )

            for idx, (pdb_atom, reference_atom) in enumerate(zip(residue_atoms, molecule.atoms)):
                pdb_element = pdb_atom.element.symbol if pdb_atom.element else None

                reference_element: str | None = None
                if hasattr(reference_atom, "element") and reference_atom.element is not None:
                    reference_element = reference_atom.element.symbol
                elif hasattr(reference_atom, "symbol"):
                    reference_element = reference_atom.symbol
                else:
                    atomic_number = getattr(reference_atom, "atomic_number", None)
                    if atomic_number:
                        element_obj = app.element.Element.getByAtomicNumber(int(atomic_number))
                        reference_element = element_obj.symbol if element_obj is not None else None

                if pdb_element is None or reference_element is None:
                    continue
                if pdb_element.upper() != reference_element.upper():
                    display_name = target.residue.name.strip() or normalized_residue_name
                    pdb_atom_name = pdb_atom.name.strip() or f"atom {idx + 1}"
                    ref_atom_name = reference_atom.name or f"atom {idx + 1}"
                    raise ValueError(
                        f"Residue '{display_name}' in state '{target.state_name}' has atom '{pdb_atom_name}' "
                        f"with element '{pdb_element}', but the SDF reference atom '{ref_atom_name}' at the same "
                        f"position is '{reference_element}'. Ensure the ligand SDF/PDB atom ordering matches."
                    )

            atom_index_map = {atom: idx for idx, atom in enumerate(residue_atoms)}
            modeller = target.modeller
            existing_pairs: Set[Tuple[int, int]] = set()
            for bond in modeller.topology.bonds():
                if bond[0].residue is target.residue and bond[1].residue is target.residue:
                    pair = tuple(sorted((atom_index_map[bond[0]], atom_index_map[bond[1]])))
                    existing_pairs.add(pair)

            added = False
            for ref_i, ref_j in reference_bonds:
                pair = tuple(sorted((ref_i, ref_j)))
                if pair in existing_pairs:
                    continue
                modeller.topology.addBond(residue_atoms[ref_i], residue_atoms[ref_j])
                existing_pairs.add(pair)
                added = True

            if added:
                modified_states.add(target.state_name)

        for state_name in modified_states:
            modeller = self._modellers[state_name]
            self._topology_bonds[state_name] = tuple(tuple(bond) for bond in modeller.topology.bonds())

        return modified_states

    def _get_state_systems(self) -> Mapping[str, openmm.System]:
        """
        Return the cached per-state OpenMM systems, building them once if needed.
        """

        if self._state_systems_cache is None:
            self._state_systems_cache = self.build_state_systems()
        return self._state_systems_cache

    def build_state_systems(self, **system_kwargs) -> Mapping[str, openmm.System]:
        """
        Construct OpenMM System objects for each EVB state using the configured force field.

        The EVB method currently operates without periodic boundary conditions (gas phase);
        this helper intentionally mirrors that limitation until spherical boundaries are
        implemented.

        Parameters
        ----------
        **system_kwargs :
            Keyword arguments forwarded to :meth:`openmm.app.ForceField.createSystem`. These options
            are applied uniformly to every state. Examples include ``nonbondedMethod=app.PME`` or
            ``constraints=app.HBonds``.

        Returns
        -------
        Mapping[str, openmm.System]
            Read-only mapping of state name -> OpenMM System.

        Notes
        -----
        Systems are cached per unique set of ``system_kwargs``. Rebuilding with the same arguments
        returns the cached objects, while different arguments trigger reconstruction.
        """

        signature = tuple(sorted(system_kwargs.items(), key=lambda entry: entry[0]))
        if self._state_systems_are_prebuilt:
            if signature != (self._state_system_signature or tuple()):
                raise ValueError(
                    "This EVBSystem was initialized from serialized state systems; rebuilding with "
                    "different createSystem options is not supported. Re-run prepare to change the "
                    "pure-state force-field setup."
                )
            self._prune_cached_state_systems()
            if self._state_systems_cache is None:
                self._state_systems_cache = MappingProxyType(self._state_systems)
            return self._state_systems_cache

        if self._missing_residues:
            pending = "; ".join(
                f"{state}: {', '.join(residues)}" for state, residues in sorted(self._missing_residues.items())
            )
            raise ValueError(
                "Cannot build OpenMM systems until all residues are parameterized. Missing -> " + pending
            )

        if self._state_systems and signature == self._state_system_signature:
            self._prune_cached_state_systems()
            if self._state_systems_cache is None:
                self._state_systems_cache = MappingProxyType(self._state_systems)
            return self._state_systems_cache

        built: Dict[str, openmm.System] = {}
        opls = self.spec_is_opls(self.force_field_name)
        for state_name, modeller in self._modellers.items():
            # fresh ForceField per state to avoid cross-state side effects
            ff = self._clone_force_field_with_generators()
            system = ff.createSystem(modeller.topology, **system_kwargs)
            if opls:
                # OPLS-AA/M needs two post-createSystem steps Amber doesn't (see forces/opls.py):
                # enumerate the periodic impropers OpenMM's matcher can't apply, and move LJ out of
                # the Script's geometric CustomNonbondedForce into the NonbondedForce so the vdW
                # partition builder can read per-atom sigma/epsilon (geometric mixing is re-applied
                # via lj_combining='geometric', auto-selected for the OPLS preset).
                from .forces import add_opls_impropers, flatten_opls_lj_to_nonbonded
                add_opls_impropers(system, modeller.topology, ff,
                                   self.resolve_force_field_files(self.force_field_name))
                flatten_opls_lj_to_nonbonded(system)
                # Full-OPLS ligand: a novel ligand comes in via SMIRNOFF (to satisfy the
                # topology), then its LJ + bonded + 1-4 are overwritten with OPLS-AA/M here so
                # the whole system is genuinely OPLS (charges left to the ligand_charges/SDF
                # machinery, keeping the charge set an orthogonal axis).
                self._apply_opls_ligand_overrides(system, modeller.topology, state_name)
            self._prune_trivial_nonbonded_forces(system)

            built[state_name] = system

        self._state_systems = built
        self._state_system_signature = signature
        self._state_systems_are_prebuilt = False
        self.bonded_force_comparisons = None
        self._invalidate_partition_caches()
        self._prune_cached_state_systems()
        self._state_systems_cache = MappingProxyType(self._state_systems)
        return self._state_systems_cache

    def _find_nonbonded_force(self, system: openmm.System) -> openmm.NonbondedForce:
        return find_nonbonded_force(system)

    def _collect_vdw_terms_by_state(
        self, systems: Mapping[str, openmm.System]
    ) -> Tuple[
        Dict[str, Dict[AtomReference, Tuple[float, float]]],
        Dict[str, Dict[CanonicalPairKey, Tuple[float, float]]],
        Dict[str, Tuple[CanonicalPairKey, ...]],
    ]:
        return collect_vdw_terms_by_state(
            systems,
            find_nonbonded_force=self._find_nonbonded_force,
            extract_vdw_particle_parameters=self._extract_vdw_particle_parameters,
            extract_vdw_exception_parameters=self._extract_vdw_exception_parameters,
        )

    def _extract_vdw_particle_parameters(
        self, state_name: str, force: openmm.NonbondedForce
    ) -> Dict[AtomReference, Tuple[float, float]]:
        return extract_vdw_particle_parameters(
            state_name,
            force,
            canonical_atom_label=self._canonical_atom_label,
        )

    def _extract_vdw_exception_parameters(
        self, state_name: str, force: openmm.NonbondedForce
    ) -> Tuple[Dict[CanonicalPairKey, Tuple[float, float]], Tuple[CanonicalPairKey, ...]]:
        return extract_vdw_exception_parameters(
            state_name,
            force,
            canonical_atom_label=self._canonical_atom_label,
            canonical_pair_key=self._canonical_pair_key,
        )

    def _collect_electrostatics_terms_by_state(
        self, systems: Mapping[str, openmm.System]
    ) -> Tuple[
        Dict[str, Dict[AtomReference, float]],
        Dict[str, Dict[CanonicalPairKey, float]],
        Dict[str, Tuple[CanonicalPairKey, ...]],
    ]:
        return collect_electrostatics_terms_by_state(
            systems,
            find_nonbonded_force=self._find_nonbonded_force,
            extract_electrostatics_particle_parameters=self._extract_electrostatics_particle_parameters,
            extract_electrostatics_exception_parameters=self._extract_electrostatics_exception_parameters,
        )

    def _extract_electrostatics_particle_parameters(
        self, state_name: str, force: openmm.NonbondedForce
    ) -> Dict[AtomReference, float]:
        return extract_electrostatics_particle_parameters(
            state_name,
            force,
            canonical_atom_label=self._canonical_atom_label,
        )

    def _extract_electrostatics_exception_parameters(
        self, state_name: str, force: openmm.NonbondedForce
    ) -> Tuple[Dict[CanonicalPairKey, float], Tuple[CanonicalPairKey, ...]]:
        return extract_electrostatics_exception_parameters(
            state_name,
            force,
            canonical_atom_label=self._canonical_atom_label,
            canonical_pair_key=self._canonical_pair_key,
        )

    def _partition_vdw_atoms(
        self,
        particle_parameters: Mapping[str, Mapping[AtomReference, Tuple[float, float]]],
        *,
        tolerance: float,
    ) -> Tuple[Dict[AtomReference, Tuple[float, float]], Dict[str, Dict[AtomReference, Tuple[float, float]]]]:
        return partition_vdw_atoms(
            particle_parameters,
            tolerance=tolerance,
            parameter_sequences_close=self._parameter_sequences_close,
        )

    def _partition_electrostatics_atoms(
        self,
        particle_parameters: Mapping[str, Mapping[AtomReference, float]],
        *,
        tolerance: float,
    ) -> Tuple[Dict[AtomReference, float], Dict[str, Dict[AtomReference, float]]]:
        return partition_electrostatics_atoms(
            particle_parameters,
            tolerance=tolerance,
            parameter_sequences_close=self._parameter_sequences_close,
        )

    def _pair_is_universal(self, pair: CanonicalPairKey, universal_atoms: Mapping[AtomReference, object]) -> bool:
        return pair_is_universal(pair, universal_atoms)

    def _partition_vdw_exceptions(
        self,
        exception_parameters: Mapping[str, Mapping[CanonicalPairKey, Tuple[float, float]]],
        *,
        tolerance: float,
        universal_atoms: Mapping[AtomReference, Tuple[float, float]],
    ) -> Tuple[Dict[CanonicalPairKey, Tuple[float, float]], Dict[str, Dict[CanonicalPairKey, Tuple[float, float]]]]:
        return partition_vdw_exceptions(
            exception_parameters,
            tolerance=tolerance,
            universal_atoms=universal_atoms,
            parameter_sequences_close=self._parameter_sequences_close,
            pair_is_universal=self._pair_is_universal,
        )

    def _partition_electrostatics_exceptions(
        self,
        exception_parameters: Mapping[str, Mapping[CanonicalPairKey, float]],
        *,
        tolerance: float,
        universal_atoms: Mapping[AtomReference, float],
    ) -> Tuple[Dict[CanonicalPairKey, float], Dict[str, Dict[CanonicalPairKey, float]]]:
        return partition_electrostatics_exceptions(
            exception_parameters,
            tolerance=tolerance,
            universal_atoms=universal_atoms,
            parameter_sequences_close=self._parameter_sequences_close,
            pair_is_universal=self._pair_is_universal,
        )

    @staticmethod
    def _partition_zero_exclusions(
        zero_pairs: Mapping[str, Tuple[CanonicalPairKey, ...]]
    ) -> Tuple[Tuple[CanonicalPairKey, ...], Dict[str, Tuple[CanonicalPairKey, ...]]]:
        return partition_zero_exclusions(zero_pairs)

    def _demote_atoms_for_exception_pairs(
        self,
        universal_atoms: Dict[AtomReference, object],
        per_state_atoms: Dict[str, Dict[AtomReference, object]],
        particle_parameters: Mapping[str, Mapping[AtomReference, object]],
        state_specific_exceptions: Mapping[str, Mapping[CanonicalPairKey, object]],
        state_specific_zero_exclusions: Mapping[str, Tuple[CanonicalPairKey, ...]],
    ) -> None:
        demote_atoms_for_exception_pairs(
            universal_atoms,
            per_state_atoms,
            particle_parameters,
            state_specific_exceptions,
            state_specific_zero_exclusions,
        )

    def _vdw_partition_terms(self, *, tolerance: float) -> VanDerWaalsPartition:
        systems = self._get_state_systems()
        cache_key = (self._state_system_signature, float(tolerance))
        if self._vdw_partition is not None and self._vdw_partition_cache_key == cache_key:
            return self._vdw_partition

        particle_parameters, exception_parameters, zero_exclusions = self._collect_vdw_terms_by_state(systems)
        universal_atoms, per_state_atoms = self._partition_vdw_atoms(particle_parameters, tolerance=tolerance)
        universal_exceptions, per_state_exceptions = self._partition_vdw_exceptions(
            exception_parameters,
            tolerance=tolerance,
            universal_atoms=universal_atoms,
        )
        universal_zero_exclusions, state_specific_zero_exclusions = self._partition_zero_exclusions(zero_exclusions)
        self._demote_atoms_for_exception_pairs(
            universal_atoms,
            per_state_atoms,
            particle_parameters,
            per_state_exceptions,
            state_specific_zero_exclusions,
        )
        partition = VanDerWaalsPartition(
            universal_atoms=universal_atoms,
            state_specific_atoms=per_state_atoms,
            universal_exceptions=universal_exceptions,
            state_specific_exceptions=per_state_exceptions,
            universal_zero_exclusions=universal_zero_exclusions,
            state_specific_zero_exclusions=state_specific_zero_exclusions,
        )
        self._vdw_partition = partition
        self._vdw_partition_cache_key = cache_key
        return partition

    def _electrostatics_partition_terms(self, *, tolerance: float) -> ElectrostaticsPartition:
        systems = self._get_state_systems()
        cache_key = (self._state_system_signature, float(tolerance))
        if self._electrostatics_partition is not None and self._electrostatics_partition_cache_key == cache_key:
            return self._electrostatics_partition

        particle_charges, exception_parameters, zero_exclusions = self._collect_electrostatics_terms_by_state(systems)
        universal_atoms, per_state_atoms = self._partition_electrostatics_atoms(particle_charges, tolerance=tolerance)
        universal_exceptions, per_state_exceptions = self._partition_electrostatics_exceptions(
            exception_parameters,
            tolerance=tolerance,
            universal_atoms=universal_atoms,
        )
        universal_zero_exclusions, state_specific_zero_exclusions = self._partition_zero_exclusions(zero_exclusions)
        self._demote_atoms_for_exception_pairs(
            universal_atoms,
            per_state_atoms,
            particle_charges,
            per_state_exceptions,
            state_specific_zero_exclusions,
        )
        partition = ElectrostaticsPartition(
            universal_atoms=universal_atoms,
            state_specific_atoms=per_state_atoms,
            universal_exceptions=universal_exceptions,
            state_specific_exceptions=per_state_exceptions,
            universal_zero_exclusions=universal_zero_exclusions,
            state_specific_zero_exclusions=state_specific_zero_exclusions,
        )
        self._electrostatics_partition = partition
        self._electrostatics_partition_cache_key = cache_key
        return partition

    def compare_bonded_terms(self, *, absolute_tolerance: float = 1e-6) -> Mapping[str, BondedForceComparison]:
        """
        Compare bonded force-field terms across all EVB states.

        Parameters
        ----------
        absolute_tolerance :
            Numerical tolerance (in native units) used when determining whether parameters
            are identical across states.

        Returns
        -------
        Mapping[str, BondedForceComparison]
            Mapping of bonded force label -> comparison summary. Each summary separates
            atom tuples into three categories:

            * ``universal_identical`` – present in every state with parameters matching
              within ``absolute_tolerance``.
            * ``universal_conflicting`` – present in every state but at least one state
              differs in its parameters.
            * ``partial`` – present in only a subset of states (each entry lists which
              states include that bonded term and their parameters).

        Notes
        -----
        Atom identities are normalized via the user-supplied mapping (when provided).
        This ensures that bonded interactions are compared between states even when atom
        indices differ across the corresponding PDB files.
        """

        systems = self._get_state_systems()
        terms_by_force = self._collect_bonded_terms_by_state(systems)
        comparisons: Dict[str, BondedForceComparison] = {}
        for force_name, state_terms in terms_by_force.items():
            comparisons[force_name] = self._summarize_bonded_terms(
                force_name,
                state_terms,
                tolerance=absolute_tolerance,
            )
        return comparisons

    def _collect_bonded_terms_by_state(
        self, systems: Mapping[str, openmm.System]
    ) -> Dict[str, Dict[str, Dict[CanonicalAtomKey, BondedTermEntry]]]:
        return collect_bonded_terms_by_state(
            systems=systems,
            canonical_bond_key_from_indices=self._canonical_bond_key_from_indices,
            canonical_angle_key_from_indices=self._canonical_angle_key_from_indices,
            canonical_torsion_key_from_indices=self._canonical_torsion_key_from_indices,
        )

    def _summarize_bonded_terms(
        self,
        force_name: str,
        terms_by_state: Mapping[str, Dict[CanonicalAtomKey, BondedTermEntry]],
        *,
        tolerance: float,
    ) -> BondedForceComparison:
        return summarize_bonded_terms(
            force_name=force_name,
            terms_by_state=terms_by_state,
            tolerance=tolerance,
        )

    @staticmethod
    def _parameter_sequences_close(left, right, tolerance: float) -> bool:
        return parameter_sequences_close(left, right, tolerance)
