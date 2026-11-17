from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from types import MappingProxyType
from typing import Dict, Tuple

import openmm
from openmm import unit

from ..models import BondedForceComparison, BondedTermEntry
from ..types import AtomReference, CanonicalAtomKey, CanonicalPairKey


def build_shared_bond_force(
    *,
    shared_entries: Mapping[CanonicalAtomKey, BondedTermEntry],
    initial_state: str,
    particle_indices_from_canonical_key: Callable[[str, CanonicalAtomKey], tuple[int, ...]],
) -> openmm.HarmonicBondForce:
    force = openmm.HarmonicBondForce()
    for canonical_key, entry in shared_entries.items():
        idx_a, idx_b = particle_indices_from_canonical_key(initial_state, canonical_key)
        length, k = entry.parameters
        force.addBond(
            idx_a,
            idx_b,
            length * unit.nanometer,
            k * unit.kilojoules_per_mole / unit.nanometer**2,
        )
    return force


def build_state_bond_force(
    *,
    lambda_param: str,
    entries: Mapping[CanonicalAtomKey, BondedTermEntry],
    initial_state: str,
    particle_indices_from_canonical_key: Callable[[str, CanonicalAtomKey], tuple[int, ...]],
) -> openmm.CustomBondForce:
    force = openmm.CustomBondForce(f"{lambda_param} * 0.5 * k * (r - r0)^2")
    force.addGlobalParameter(lambda_param, 0.0)
    force.addPerBondParameter("r0")
    force.addPerBondParameter("k")
    for canonical_key, entry in entries.items():
        idx_a, idx_b = particle_indices_from_canonical_key(initial_state, canonical_key)
        length, k_value = entry.parameters
        force.addBond(idx_a, idx_b, (length, k_value))
    return force


def build_shared_angle_force(
    *,
    shared_entries: Mapping[CanonicalAtomKey, BondedTermEntry],
    initial_state: str,
    particle_indices_from_canonical_key: Callable[[str, CanonicalAtomKey], tuple[int, ...]],
) -> openmm.HarmonicAngleForce:
    force = openmm.HarmonicAngleForce()
    for canonical_key, entry in shared_entries.items():
        idx_a, idx_b, idx_c = particle_indices_from_canonical_key(initial_state, canonical_key)
        theta0, k = entry.parameters
        force.addAngle(
            idx_a,
            idx_b,
            idx_c,
            theta0 * unit.radian,
            k * unit.kilojoules_per_mole / unit.radian**2,
        )
    return force


def build_state_angle_force(
    *,
    lambda_param: str,
    entries: Mapping[CanonicalAtomKey, BondedTermEntry],
    initial_state: str,
    particle_indices_from_canonical_key: Callable[[str, CanonicalAtomKey], tuple[int, ...]],
) -> openmm.CustomAngleForce:
    force = openmm.CustomAngleForce(f"{lambda_param} * 0.5 * k * (theta - theta0)^2")
    force.addGlobalParameter(lambda_param, 0.0)
    force.addPerAngleParameter("theta0")
    force.addPerAngleParameter("k")
    for canonical_key, entry in entries.items():
        idx_a, idx_b, idx_c = particle_indices_from_canonical_key(initial_state, canonical_key)
        theta0_value, k_value = entry.parameters
        force.addAngle(idx_a, idx_b, idx_c, (theta0_value, k_value))
    return force


def build_shared_torsion_force(
    *,
    shared_entries: Mapping[CanonicalAtomKey, BondedTermEntry],
    initial_state: str,
    particle_indices_from_canonical_key: Callable[[str, CanonicalAtomKey], tuple[int, ...]],
    torsion_components: Callable[[Tuple[object, ...]], tuple[tuple[int, float, float], ...]],
) -> openmm.PeriodicTorsionForce:
    force = openmm.PeriodicTorsionForce()
    for canonical_key, entry in shared_entries.items():
        idx_a, idx_b, idx_c, idx_d = particle_indices_from_canonical_key(initial_state, canonical_key)
        for periodicity, phase, k in torsion_components(entry.parameters):
            force.addTorsion(
                idx_a,
                idx_b,
                idx_c,
                idx_d,
                periodicity,
                phase * unit.radian,
                k * unit.kilojoules_per_mole,
            )
    return force


def build_state_torsion_force(
    *,
    lambda_param: str,
    entries: Mapping[CanonicalAtomKey, BondedTermEntry],
    initial_state: str,
    particle_indices_from_canonical_key: Callable[[str, CanonicalAtomKey], tuple[int, ...]],
    torsion_components: Callable[[Tuple[object, ...]], tuple[tuple[int, float, float], ...]],
) -> openmm.CustomTorsionForce:
    force = openmm.CustomTorsionForce(f"{lambda_param} * k * (1 + cos(periodicity*theta - phase))")
    force.addGlobalParameter(lambda_param, 0.0)
    force.addPerTorsionParameter("periodicity")
    force.addPerTorsionParameter("phase")
    force.addPerTorsionParameter("k")
    for canonical_key, entry in entries.items():
        idx_a, idx_b, idx_c, idx_d = particle_indices_from_canonical_key(initial_state, canonical_key)
        for periodicity, phase_value, k_value in torsion_components(entry.parameters):
            force.addTorsion(idx_a, idx_b, idx_c, idx_d, (periodicity, phase_value, k_value))
    return force


def extract_harmonic_bond_terms(
    *,
    state_name: str,
    force: openmm.HarmonicBondForce,
    canonical_bond_key_from_indices: Callable[[str, int, int], CanonicalAtomKey],
) -> dict[CanonicalAtomKey, BondedTermEntry]:
    terms: dict[CanonicalAtomKey, BondedTermEntry] = {}
    for idx in range(force.getNumBonds()):
        atom1, atom2, length, k = force.getBondParameters(idx)
        key = canonical_bond_key_from_indices(state_name, atom1, atom2)
        params = (
            length.value_in_unit(unit.nanometer),
            k.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2),
        )
        terms[key] = BondedTermEntry(parameters=params, original_indices=(atom1, atom2))
    return terms


def extract_harmonic_angle_terms(
    *,
    state_name: str,
    force: openmm.HarmonicAngleForce,
    canonical_angle_key_from_indices: Callable[[str, int, int, int], CanonicalAtomKey],
) -> dict[CanonicalAtomKey, BondedTermEntry]:
    terms: dict[CanonicalAtomKey, BondedTermEntry] = {}
    for idx in range(force.getNumAngles()):
        atom1, atom2, atom3, angle, k = force.getAngleParameters(idx)
        key = canonical_angle_key_from_indices(state_name, atom1, atom2, atom3)
        params = (
            angle.value_in_unit(unit.radian),
            k.value_in_unit(unit.kilojoules_per_mole / unit.radian**2),
        )
        terms[key] = BondedTermEntry(parameters=params, original_indices=(atom1, atom2, atom3))
    return terms


def extract_periodic_torsion_terms(
    *,
    state_name: str,
    force: openmm.PeriodicTorsionForce,
    canonical_torsion_key_from_indices: Callable[[str, int, int, int, int], CanonicalAtomKey],
) -> dict[CanonicalAtomKey, BondedTermEntry]:
    grouped: dict[CanonicalAtomKey, list[tuple[int, float, float]]] = {}
    originals: dict[CanonicalAtomKey, tuple[int, int, int, int]] = {}
    for idx in range(force.getNumTorsions()):
        atom1, atom2, atom3, atom4, periodicity, phase, k = force.getTorsionParameters(idx)
        key = canonical_torsion_key_from_indices(state_name, atom1, atom2, atom3, atom4)
        params = (
            int(periodicity),
            phase.value_in_unit(unit.radian),
            k.value_in_unit(unit.kilojoules_per_mole),
        )
        grouped.setdefault(key, []).append(params)
        originals.setdefault(key, (atom1, atom2, atom3, atom4))

    terms: dict[CanonicalAtomKey, BondedTermEntry] = {}
    for key, entries in grouped.items():
        terms[key] = BondedTermEntry(parameters=tuple(sorted(entries)), original_indices=originals[key])
    return terms


def collect_bonded_terms_by_state(
    *,
    systems: Mapping[str, openmm.System],
    canonical_bond_key_from_indices: Callable[[str, int, int], CanonicalAtomKey],
    canonical_angle_key_from_indices: Callable[[str, int, int, int], CanonicalAtomKey],
    canonical_torsion_key_from_indices: Callable[[str, int, int, int, int], CanonicalAtomKey],
) -> dict[str, dict[str, dict[CanonicalAtomKey, BondedTermEntry]]]:
    force_labels = ("bonds", "angles", "torsions")
    per_force: dict[str, dict[str, dict[CanonicalAtomKey, BondedTermEntry]]] = {
        label: {state: {} for state in systems.keys()} for label in force_labels
    }

    for state_name, system in systems.items():
        for force in system.getForces():
            if isinstance(force, openmm.HarmonicBondForce):
                per_force["bonds"][state_name].update(
                    extract_harmonic_bond_terms(
                        state_name=state_name,
                        force=force,
                        canonical_bond_key_from_indices=canonical_bond_key_from_indices,
                    )
                )
            elif isinstance(force, openmm.HarmonicAngleForce):
                per_force["angles"][state_name].update(
                    extract_harmonic_angle_terms(
                        state_name=state_name,
                        force=force,
                        canonical_angle_key_from_indices=canonical_angle_key_from_indices,
                    )
                )
            elif isinstance(force, openmm.PeriodicTorsionForce):
                per_force["torsions"][state_name].update(
                    extract_periodic_torsion_terms(
                        state_name=state_name,
                        force=force,
                        canonical_torsion_key_from_indices=canonical_torsion_key_from_indices,
                    )
                )
    return per_force


def parameter_sequences_close(left: object, right: object, tolerance: float) -> bool:
    if isinstance(left, tuple) and isinstance(right, tuple):
        if len(left) != len(right):
            return False
        return all(parameter_sequences_close(l_entry, r_entry, tolerance) for l_entry, r_entry in zip(left, right))
    try:
        left_value = float(left)
        right_value = float(right)
    except (TypeError, ValueError):
        return False
    return abs(left_value - right_value) <= tolerance


def summarize_bonded_terms(
    *,
    force_name: str,
    terms_by_state: Mapping[str, dict[CanonicalAtomKey, BondedTermEntry]],
    tolerance: float,
) -> BondedForceComparison:
    all_keys: set[CanonicalAtomKey] = set()
    for terms in terms_by_state.values():
        all_keys.update(terms.keys())

    universal_identical: dict[CanonicalAtomKey, BondedTermEntry] = {}
    universal_conflicting: dict[CanonicalAtomKey, dict[str, BondedTermEntry]] = {}
    partial: dict[CanonicalAtomKey, dict[str, BondedTermEntry]] = {}

    state_names = tuple(terms_by_state.keys())
    for key in sorted(all_keys):
        state_params = {state: entries[key] for state, entries in terms_by_state.items() if key in entries}
        if len(state_params) == len(state_names):
            iterator = iter(state_params.values())
            first_entry = next(iterator)
            if all(
                parameter_sequences_close(first_entry.parameters, other_entry.parameters, tolerance)
                for other_entry in iterator
            ):
                universal_identical[key] = first_entry
            else:
                universal_conflicting[key] = state_params
        else:
            partial[key] = state_params

    return BondedForceComparison(
        force_name=force_name,
        universal_identical=universal_identical,
        universal_conflicting=universal_conflicting,
        partial=partial,
    )


def extract_harmonic_bond_terms_by_index(system: openmm.System) -> dict[tuple[int, int], BondedTermEntry]:
    terms: dict[tuple[int, int], BondedTermEntry] = {}
    for force in system.getForces():
        if not isinstance(force, openmm.HarmonicBondForce):
            continue
        for idx in range(force.getNumBonds()):
            atom1, atom2, length, k = force.getBondParameters(idx)
            pair = (atom1, atom2) if atom1 <= atom2 else (atom2, atom1)
            params = (
                length.value_in_unit(unit.nanometer),
                k.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2),
            )
            terms[pair] = BondedTermEntry(parameters=params, original_indices=(atom1, atom2))
    return terms


def states_share_directional_index_layout(
    *,
    reference_atoms: Sequence[AtomReference],
    comparison_atom_sets: Sequence[Sequence[AtomReference]],
) -> bool:
    reference_labels = tuple((entry[1], entry[2], entry[3]) for entry in reference_atoms)
    for candidate_atoms in comparison_atom_sets:
        if len(candidate_atoms) != len(reference_labels):
            return False
        candidate_labels = tuple((entry[1], entry[2], entry[3]) for entry in candidate_atoms)
        if candidate_labels != reference_labels:
            return False
    return True


def collect_forming_breaking_bond_pairs_by_index(
    *,
    reference_state: str,
    comparison_states: Sequence[str],
    state_atom_references: Mapping[str, Sequence[AtomReference]],
    state_bond_terms_by_index: Mapping[str, Mapping[tuple[int, int], BondedTermEntry]],
    atom_reference_from_index: Callable[[str, int], AtomReference],
    canonical_pair_key: Callable[[AtomReference, AtomReference], CanonicalPairKey],
) -> Mapping[CanonicalPairKey, Mapping[str, BondedTermEntry]]:
    state_group = (reference_state, *comparison_states)
    if not states_share_directional_index_layout(
        reference_atoms=state_atom_references[reference_state],
        comparison_atom_sets=[state_atom_references[state_name] for state_name in comparison_states],
    ):
        return MappingProxyType({})

    reference_atom_count = len(state_atom_references[reference_state])
    for state_name in comparison_states:
        if len(state_atom_references[state_name]) != reference_atom_count:
            return MappingProxyType({})

    all_pairs: set[tuple[int, int]] = set()
    for state_name in state_group:
        all_pairs.update(state_bond_terms_by_index[state_name].keys())

    reactive_pairs: dict[CanonicalPairKey, Mapping[str, BondedTermEntry]] = {}
    for pair in sorted(all_pairs):
        entries = {
            state_name: state_bond_terms_by_index[state_name][pair]
            for state_name in state_group
            if pair in state_bond_terms_by_index[state_name]
        }
        if not entries or len(entries) == len(state_group):
            continue
        atom1_ref = atom_reference_from_index(reference_state, pair[0])
        atom2_ref = atom_reference_from_index(reference_state, pair[1])
        key = canonical_pair_key(atom1_ref, atom2_ref)
        reactive_pairs[key] = MappingProxyType(dict(entries))
    return MappingProxyType(reactive_pairs)


def group_bonded_terms_by_pair(
    *,
    bonds_by_state: Mapping[str, Mapping[CanonicalAtomKey, BondedTermEntry]],
    canonical_pair_key: Callable[[AtomReference, AtomReference], CanonicalPairKey],
) -> Mapping[CanonicalPairKey, Mapping[str, BondedTermEntry]]:
    grouped: dict[CanonicalPairKey, Dict[str, BondedTermEntry]] = {}
    for state_name, state_terms in bonds_by_state.items():
        for canonical_key, entry in state_terms.items():
            if len(canonical_key) != 2:
                continue
            pair = canonical_pair_key(canonical_key[0], canonical_key[1])
            grouped.setdefault(pair, {})[state_name] = entry
    return MappingProxyType({pair: MappingProxyType(dict(entries)) for pair, entries in grouped.items()})


__all__ = [
    "build_shared_bond_force",
    "build_state_bond_force",
    "build_shared_angle_force",
    "build_state_angle_force",
    "build_shared_torsion_force",
    "build_state_torsion_force",
    "extract_harmonic_bond_terms",
    "extract_harmonic_angle_terms",
    "extract_periodic_torsion_terms",
    "collect_bonded_terms_by_state",
    "parameter_sequences_close",
    "summarize_bonded_terms",
    "extract_harmonic_bond_terms_by_index",
    "states_share_directional_index_layout",
    "collect_forming_breaking_bond_pairs_by_index",
    "group_bonded_terms_by_pair",
]
