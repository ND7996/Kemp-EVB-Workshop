from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from numbers import Real
from typing import Dict

import openmm
from openmm import unit

from ..types import CanonicalPairKey


def build_evb_system_scaffold(
    *,
    initial_state: str,
    state_particle_masses: Callable[[str], Sequence[unit.Quantity]],
    add_shared_constraints_to_evb_system: Callable[[openmm.System, str], int],
) -> openmm.System:
    system = openmm.System()
    for mass in state_particle_masses(initial_state):
        system.addParticle(mass)
    add_shared_constraints_to_evb_system(system, initial_state)
    return system


def constraint_distance_nm(value: unit.Quantity | float) -> float:
    if isinstance(value, Real):
        return float(value)
    return float(value.value_in_unit(unit.nanometer))


def state_constraint_map(
    *,
    state_name: str,
    system: openmm.System,
    canonical_pair_key_from_indices: Callable[[str, int, int], CanonicalPairKey],
    constraint_distance_nm: Callable[[unit.Quantity | float], float],
) -> Dict[CanonicalPairKey, float]:
    mapping: Dict[CanonicalPairKey, float] = {}
    for idx in range(system.getNumConstraints()):
        atom1, atom2, distance = system.getConstraintParameters(idx)
        pair = canonical_pair_key_from_indices(state_name, int(atom1), int(atom2))
        distance_nm = constraint_distance_nm(distance)
        previous = mapping.get(pair)
        if previous is not None and abs(previous - distance_nm) > 1e-6:
            raise ValueError(
                f"State '{state_name}' defines inconsistent constraint distances for pair {pair!r}."
            )
        mapping[pair] = distance_nm
    return mapping


def add_shared_constraints_to_evb_system(
    *,
    system: openmm.System,
    initial_state: str,
    state_names: Sequence[str],
    reference_systems: Mapping[str, openmm.System],
    state_constraint_map: Callable[[str, openmm.System], Dict[CanonicalPairKey, float]],
    particle_indices_from_canonical_pair: Callable[[str, CanonicalPairKey], tuple[int, int]],
    tolerance_nm: float = 1e-6,
) -> int:
    if not reference_systems:
        return 0

    per_state_constraints: Dict[str, Dict[CanonicalPairKey, float]] = {}
    for state_name in state_names:
        state_system = reference_systems[state_name]
        per_state_constraints[state_name] = state_constraint_map(state_name, state_system)

    common_pairs = set(per_state_constraints[state_names[0]].keys())
    for state_name in state_names[1:]:
        common_pairs.intersection_update(per_state_constraints[state_name].keys())

    added = 0
    for pair in sorted(common_pairs):
        distances = [per_state_constraints[state_name][pair] for state_name in state_names]
        if max(distances) - min(distances) > float(tolerance_nm):
            continue
        atom1, atom2 = particle_indices_from_canonical_pair(initial_state, pair)
        system.addConstraint(int(atom1), int(atom2), float(sum(distances) / len(distances)) * unit.nanometer)
        added += 1
    return added


__all__ = [
    "build_evb_system_scaffold",
    "constraint_distance_nm",
    "state_constraint_map",
    "add_shared_constraints_to_evb_system",
]
