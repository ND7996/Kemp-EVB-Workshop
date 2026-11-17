from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Dict, Set, Tuple

import openmm
from openmm import unit

from ..types import AtomReference, CanonicalPairKey


def find_nonbonded_force(system: openmm.System) -> openmm.NonbondedForce:
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            return force
    raise ValueError("No NonbondedForce present in state system.")


def collect_vdw_terms_by_state(
    systems: Mapping[str, openmm.System],
    *,
    find_nonbonded_force: Callable[[openmm.System], openmm.NonbondedForce],
    extract_vdw_particle_parameters: Callable[[str, openmm.NonbondedForce], Dict[AtomReference, Tuple[float, float]]],
    extract_vdw_exception_parameters: Callable[
        [str, openmm.NonbondedForce],
        Tuple[Dict[CanonicalPairKey, Tuple[float, float]], Tuple[CanonicalPairKey, ...]],
    ],
) -> Tuple[
    Dict[str, Dict[AtomReference, Tuple[float, float]]],
    Dict[str, Dict[CanonicalPairKey, Tuple[float, float]]],
    Dict[str, Tuple[CanonicalPairKey, ...]],
]:
    particle_parameters: Dict[str, Dict[AtomReference, Tuple[float, float]]] = {}
    exception_parameters: Dict[str, Dict[CanonicalPairKey, Tuple[float, float]]] = {}
    zero_exceptions: Dict[str, Tuple[CanonicalPairKey, ...]] = {}

    for state_name, system in systems.items():
        nb_force = find_nonbonded_force(system)
        particle_parameters[state_name] = extract_vdw_particle_parameters(state_name, nb_force)
        exception_params, exclusions = extract_vdw_exception_parameters(state_name, nb_force)
        exception_parameters[state_name] = exception_params
        zero_exceptions[state_name] = exclusions

    return particle_parameters, exception_parameters, zero_exceptions


def extract_vdw_particle_parameters(
    state_name: str,
    force: openmm.NonbondedForce,
    *,
    canonical_atom_label: Callable[[str, int], AtomReference],
) -> Dict[AtomReference, Tuple[float, float]]:
    params: Dict[AtomReference, Tuple[float, float]] = {}
    for idx in range(force.getNumParticles()):
        _, sigma, epsilon = force.getParticleParameters(idx)
        label = canonical_atom_label(state_name, idx)
        params[label] = (
            sigma.value_in_unit(unit.nanometer),
            epsilon.value_in_unit(unit.kilojoules_per_mole),
        )
    return params


def extract_vdw_exception_parameters(
    state_name: str,
    force: openmm.NonbondedForce,
    *,
    canonical_atom_label: Callable[[str, int], AtomReference],
    canonical_pair_key: Callable[[AtomReference, AtomReference], CanonicalPairKey],
) -> Tuple[Dict[CanonicalPairKey, Tuple[float, float]], Tuple[CanonicalPairKey, ...]]:
    parameters: Dict[CanonicalPairKey, Tuple[float, float]] = {}
    exclusions: Set[CanonicalPairKey] = set()
    exclusion_tolerance = 1e-12
    for idx in range(force.getNumExceptions()):
        atom1, atom2, charge_prod, sigma, epsilon = force.getExceptionParameters(idx)
        label1 = canonical_atom_label(state_name, atom1)
        label2 = canonical_atom_label(state_name, atom2)
        pair_key = canonical_pair_key(label1, label2)
        epsilon_value = epsilon.value_in_unit(unit.kilojoules_per_mole)
        charge_value = charge_prod.value_in_unit(unit.elementary_charge**2)
        if abs(epsilon_value) <= exclusion_tolerance and abs(charge_value) <= exclusion_tolerance:
            exclusions.add(pair_key)
            continue
        parameters[pair_key] = (
            sigma.value_in_unit(unit.nanometer),
            epsilon_value,
        )
    return parameters, tuple(sorted(exclusions))


def collect_electrostatics_terms_by_state(
    systems: Mapping[str, openmm.System],
    *,
    find_nonbonded_force: Callable[[openmm.System], openmm.NonbondedForce],
    extract_electrostatics_particle_parameters: Callable[[str, openmm.NonbondedForce], Dict[AtomReference, float]],
    extract_electrostatics_exception_parameters: Callable[
        [str, openmm.NonbondedForce],
        Tuple[Dict[CanonicalPairKey, float], Tuple[CanonicalPairKey, ...]],
    ],
) -> Tuple[
    Dict[str, Dict[AtomReference, float]],
    Dict[str, Dict[CanonicalPairKey, float]],
    Dict[str, Tuple[CanonicalPairKey, ...]],
]:
    particle_charges: Dict[str, Dict[AtomReference, float]] = {}
    exception_parameters: Dict[str, Dict[CanonicalPairKey, float]] = {}
    zero_exceptions: Dict[str, Tuple[CanonicalPairKey, ...]] = {}

    for state_name, system in systems.items():
        nb_force = find_nonbonded_force(system)
        particle_charges[state_name] = extract_electrostatics_particle_parameters(state_name, nb_force)
        exception_params, exclusions = extract_electrostatics_exception_parameters(state_name, nb_force)
        exception_parameters[state_name] = exception_params
        zero_exceptions[state_name] = exclusions
    return particle_charges, exception_parameters, zero_exceptions


def extract_electrostatics_particle_parameters(
    state_name: str,
    force: openmm.NonbondedForce,
    *,
    canonical_atom_label: Callable[[str, int], AtomReference],
) -> Dict[AtomReference, float]:
    charges: Dict[AtomReference, float] = {}
    for idx in range(force.getNumParticles()):
        charge, _, _ = force.getParticleParameters(idx)
        label = canonical_atom_label(state_name, idx)
        charges[label] = charge.value_in_unit(unit.elementary_charge)
    return charges


def extract_electrostatics_exception_parameters(
    state_name: str,
    force: openmm.NonbondedForce,
    *,
    canonical_atom_label: Callable[[str, int], AtomReference],
    canonical_pair_key: Callable[[AtomReference, AtomReference], CanonicalPairKey],
) -> Tuple[Dict[CanonicalPairKey, float], Tuple[CanonicalPairKey, ...]]:
    parameters: Dict[CanonicalPairKey, float] = {}
    exclusions: Set[CanonicalPairKey] = set()
    exclusion_tolerance = 1e-12
    for idx in range(force.getNumExceptions()):
        atom1, atom2, charge_prod, sigma, epsilon = force.getExceptionParameters(idx)
        label1 = canonical_atom_label(state_name, atom1)
        label2 = canonical_atom_label(state_name, atom2)
        pair_key = canonical_pair_key(label1, label2)
        epsilon_value = epsilon.value_in_unit(unit.kilojoules_per_mole)
        charge_value = charge_prod.value_in_unit(unit.elementary_charge**2)
        if abs(epsilon_value) <= exclusion_tolerance and abs(charge_value) <= exclusion_tolerance:
            exclusions.add(pair_key)
            continue
        if abs(charge_value) <= exclusion_tolerance:
            continue
        parameters[pair_key] = charge_value
    return parameters, tuple(sorted(exclusions))


def partition_vdw_atoms(
    particle_parameters: Mapping[str, Mapping[AtomReference, Tuple[float, float]]],
    *,
    tolerance: float,
    parameter_sequences_close: Callable[[object, object, float], bool],
) -> Tuple[Dict[AtomReference, Tuple[float, float]], Dict[str, Dict[AtomReference, Tuple[float, float]]]]:
    universal: Dict[AtomReference, Tuple[float, float]] = {}
    per_state: Dict[str, Dict[AtomReference, Tuple[float, float]]] = {}

    state_names = tuple(particle_parameters.keys())
    all_atoms: Set[AtomReference] = set()
    for params in particle_parameters.values():
        all_atoms.update(params.keys())

    for canonical_atom in sorted(all_atoms):
        entries = {
            state: params[canonical_atom]
            for state, params in particle_parameters.items()
            if canonical_atom in params
        }
        complete = len(entries) == len(state_names)
        if complete and entries:
            iterator = iter(entries.values())
            first = next(iterator)
            if all(parameter_sequences_close(first, other, tolerance) for other in iterator):
                universal[canonical_atom] = first
                continue
        for state_name, params in entries.items():
            per_state.setdefault(state_name, {})[canonical_atom] = params
    return universal, per_state


def partition_electrostatics_atoms(
    particle_parameters: Mapping[str, Mapping[AtomReference, float]],
    *,
    tolerance: float,
    parameter_sequences_close: Callable[[object, object, float], bool],
) -> Tuple[Dict[AtomReference, float], Dict[str, Dict[AtomReference, float]]]:
    universal: Dict[AtomReference, float] = {}
    per_state: Dict[str, Dict[AtomReference, float]] = {}

    state_names = tuple(particle_parameters.keys())
    all_atoms: Set[AtomReference] = set()
    for params in particle_parameters.values():
        all_atoms.update(params.keys())

    for canonical_atom in sorted(all_atoms):
        entries = {
            state: params[canonical_atom]
            for state, params in particle_parameters.items()
            if canonical_atom in params
        }
        complete = len(entries) == len(state_names)
        if complete and entries:
            iterator = iter(entries.values())
            first = next(iterator)
            if all(parameter_sequences_close(first, other, tolerance) for other in iterator):
                universal[canonical_atom] = first
                continue
        for state_name, value in entries.items():
            per_state.setdefault(state_name, {})[canonical_atom] = value
    return universal, per_state


def pair_is_universal(pair: CanonicalPairKey, universal_atoms: Mapping[AtomReference, object]) -> bool:
    return pair[0] in universal_atoms and pair[1] in universal_atoms


def partition_vdw_exceptions(
    exception_parameters: Mapping[str, Mapping[CanonicalPairKey, Tuple[float, float]]],
    *,
    tolerance: float,
    universal_atoms: Mapping[AtomReference, Tuple[float, float]],
    parameter_sequences_close: Callable[[object, object, float], bool],
    pair_is_universal: Callable[[CanonicalPairKey, Mapping[AtomReference, object]], bool],
) -> Tuple[Dict[CanonicalPairKey, Tuple[float, float]], Dict[str, Dict[CanonicalPairKey, Tuple[float, float]]]]:
    universal: Dict[CanonicalPairKey, Tuple[float, float]] = {}
    per_state: Dict[str, Dict[CanonicalPairKey, Tuple[float, float]]] = {}

    state_names = tuple(exception_parameters.keys())
    all_pairs: Set[CanonicalPairKey] = set()
    for params in exception_parameters.values():
        all_pairs.update(params.keys())

    for pair in sorted(all_pairs):
        entries = {
            state: params[pair]
            for state, params in exception_parameters.items()
            if pair in params
        }
        if not entries:
            continue
        complete = len(entries) == len(state_names)
        add_to_universal = False
        reference_entry: Tuple[float, float] | None = None
        if complete:
            iterator = iter(entries.values())
            reference_entry = next(iterator)
            if all(parameter_sequences_close(reference_entry, other, tolerance) for other in iterator):
                add_to_universal = True
        if add_to_universal and reference_entry is not None and pair_is_universal(pair, universal_atoms):
            universal[pair] = reference_entry
        else:
            for state_name, params in entries.items():
                per_state.setdefault(state_name, {})[pair] = params
    return universal, per_state


def partition_electrostatics_exceptions(
    exception_parameters: Mapping[str, Mapping[CanonicalPairKey, float]],
    *,
    tolerance: float,
    universal_atoms: Mapping[AtomReference, float],
    parameter_sequences_close: Callable[[object, object, float], bool],
    pair_is_universal: Callable[[CanonicalPairKey, Mapping[AtomReference, object]], bool],
) -> Tuple[Dict[CanonicalPairKey, float], Dict[str, Dict[CanonicalPairKey, float]]]:
    universal: Dict[CanonicalPairKey, float] = {}
    per_state: Dict[str, Dict[CanonicalPairKey, float]] = {}

    state_names = tuple(exception_parameters.keys())
    all_pairs: Set[CanonicalPairKey] = set()
    for params in exception_parameters.values():
        all_pairs.update(params.keys())

    for pair in sorted(all_pairs):
        entries = {
            state: params[pair]
            for state, params in exception_parameters.items()
            if pair in params
        }
        if not entries:
            continue
        complete = len(entries) == len(state_names)
        add_to_universal = False
        reference_entry: float | None = None
        if complete:
            iterator = iter(entries.values())
            reference_entry = next(iterator)
            if all(parameter_sequences_close(reference_entry, other, tolerance) for other in iterator):
                add_to_universal = True
        if add_to_universal and reference_entry is not None and pair_is_universal(pair, universal_atoms):
            universal[pair] = reference_entry
        else:
            for state_name, value in entries.items():
                per_state.setdefault(state_name, {})[pair] = value
    return universal, per_state


def partition_zero_exclusions(
    zero_pairs: Mapping[str, Tuple[CanonicalPairKey, ...]]
) -> Tuple[Tuple[CanonicalPairKey, ...], Dict[str, Tuple[CanonicalPairKey, ...]]]:
    if not zero_pairs:
        return tuple(), {}
    zero_sets = {state: set(pairs) for state, pairs in zero_pairs.items()}
    iterator = iter(zero_sets.values())
    try:
        universal = set(next(iterator))
    except StopIteration:
        return tuple(), {}
    for remaining in iterator:
        universal &= remaining
    universal_sorted = tuple(sorted(universal))
    per_state: Dict[str, Tuple[CanonicalPairKey, ...]] = {}
    for state, pairs in zero_sets.items():
        difference = pairs - universal
        if difference:
            per_state[state] = tuple(sorted(difference))
    return universal_sorted, per_state


def demote_atoms_for_exception_pairs(
    universal_atoms: Dict[AtomReference, object],
    per_state_atoms: Dict[str, Dict[AtomReference, object]],
    particle_parameters: Mapping[str, Mapping[AtomReference, object]],
    state_specific_exceptions: Mapping[str, Mapping[CanonicalPairKey, object]],
    state_specific_zero_exclusions: Mapping[str, Tuple[CanonicalPairKey, ...]],
) -> None:
    atoms_to_demote: Set[AtomReference] = set()
    for per_state in state_specific_exceptions.values():
        for pair in per_state.keys():
            atoms_to_demote.update(pair)
    for pairs in state_specific_zero_exclusions.values():
        for pair in pairs:
            atoms_to_demote.update(pair)

    if not atoms_to_demote:
        return

    for atom in atoms_to_demote:
        if atom not in universal_atoms:
            continue
        for state_name, parameters in particle_parameters.items():
            if atom not in parameters:
                raise ValueError(f"Atom {atom!r} missing from particle parameters for state '{state_name}'.")
            per_state_atoms.setdefault(state_name, {})[atom] = parameters[atom]
        universal_atoms.pop(atom, None)


__all__ = [
    "find_nonbonded_force",
    "collect_vdw_terms_by_state",
    "extract_vdw_particle_parameters",
    "extract_vdw_exception_parameters",
    "collect_electrostatics_terms_by_state",
    "extract_electrostatics_particle_parameters",
    "extract_electrostatics_exception_parameters",
    "partition_vdw_atoms",
    "partition_electrostatics_atoms",
    "pair_is_universal",
    "partition_vdw_exceptions",
    "partition_electrostatics_exceptions",
    "partition_zero_exclusions",
    "demote_atoms_for_exception_pairs",
]
