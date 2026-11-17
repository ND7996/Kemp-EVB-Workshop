from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Dict, Iterable, Set, Tuple, TypeVar

import openmm

from ..models import ElectrostaticsPartition, VanDerWaalsPartition
from ..types import AtomReference, CanonicalPairKey

_T = TypeVar("_T")


def split_shared_exception_terms(
    *,
    partition: VanDerWaalsPartition,
    electrostatics_partition: ElectrostaticsPartition,
    all_exception_parameters: Mapping[CanonicalPairKey, Tuple[float, float]],
    qm_labels: set[AtomReference],
    pair_in_qm_channel: Callable[[CanonicalPairKey], bool],
) -> tuple[
    list[tuple[CanonicalPairKey, float, float, float]],
    list[tuple[CanonicalPairKey, float, float, float]],
]:
    combined_exception_pairs = set(partition.universal_exceptions.keys())
    combined_exception_pairs.update(electrostatics_partition.universal_exceptions.keys())
    mm_terms: list[tuple[CanonicalPairKey, float, float, float]] = []
    qm_terms: list[tuple[CanonicalPairKey, float, float, float]] = []
    for pair in sorted(combined_exception_pairs):
        if pair in partition.universal_exceptions:
            sigma, epsilon = partition.universal_exceptions[pair]
        elif pair in all_exception_parameters:
            sigma, epsilon = all_exception_parameters[pair]
        else:
            sigma, epsilon = 1.0, 0.0
        charge_prod = electrostatics_partition.universal_exceptions.get(pair, 0.0)
        if abs(epsilon) <= 1e-12 and abs(charge_prod) <= 1e-12:
            continue
        if pair_in_qm_channel(pair):
            qm_terms.append((pair, sigma, epsilon, charge_prod))
        else:
            mm_terms.append((pair, sigma, epsilon, charge_prod))
    return mm_terms, qm_terms


def build_shared_exception_force(
    *,
    terms: Sequence[tuple[CanonicalPairKey, float, float, float]],
    one_4pi_eps0: float,
    initial_state: str,
    particle_indices_from_canonical_pair: Callable[[str, CanonicalPairKey], tuple[int, int]],
    periodic: bool = False,
) -> openmm.CustomBondForce | None:
    if not terms:
        return None
    force = openmm.CustomBondForce(
        "4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*charge_prod/r"
    )
    force.setUsesPeriodicBoundaryConditions(periodic)
    force.addGlobalParameter("ONE_4PI_EPS0", one_4pi_eps0)
    force.addPerBondParameter("sigma")
    force.addPerBondParameter("epsilon")
    force.addPerBondParameter("charge_prod")
    for pair, sigma, epsilon, charge_prod in terms:
        idx_a, idx_b = particle_indices_from_canonical_pair(initial_state, pair)
        force.addBond(idx_a, idx_b, (sigma, epsilon, charge_prod))
    return force


def collect_state_exclusion_pairs(
    *,
    partition: VanDerWaalsPartition,
    electrostatics_partition: ElectrostaticsPartition,
    state_name: str,
    exception_parameters: Mapping[CanonicalPairKey, object],
    electro_exception_parameters: Mapping[CanonicalPairKey, object],
) -> Set[CanonicalPairKey]:
    state_exclusion_pairs: Set[CanonicalPairKey] = set(partition.universal_zero_exclusions)
    state_exclusion_pairs.update(electrostatics_partition.universal_zero_exclusions)
    state_exclusion_pairs.update(partition.universal_exceptions.keys())
    state_exclusion_pairs.update(electrostatics_partition.universal_exceptions.keys())
    state_exclusion_pairs.update(partition.state_specific_zero_exclusions.get(state_name, ()))
    state_exclusion_pairs.update(electrostatics_partition.state_specific_zero_exclusions.get(state_name, ()))
    state_exclusion_pairs.update(exception_parameters.keys())
    state_exclusion_pairs.update(electro_exception_parameters.keys())
    return state_exclusion_pairs


def populate_state_combined_particles(
    *,
    force: openmm.CustomNonbondedForce,
    force_qm: openmm.CustomNonbondedForce | None,
    reference_labels: Sequence[AtomReference],
    vdw_atom_parameters: Mapping[AtomReference, Tuple[float, float]],
    vdw_universal_atoms: Mapping[AtomReference, Tuple[float, float]],
    electro_atom_parameters: Mapping[AtomReference, float],
    electro_universal_atoms: Mapping[AtomReference, float],
    qm_labels: set[AtomReference],
) -> Set[int]:
    """Populate a state's combined vdW+electrostatics CustomNonbondedForce.

    Each particle receives ``(sigma, epsilon, charge, active_vdw,
    active_elec, qm)``. Returns the union of vdW-active and
    electrostatics-active indices, which is the right subset to pass to
    ``restrict_to_subset_interaction_group`` because a pair contributes
    iff at least one atom is active in **either** term.
    """
    active_indices: Set[int] = set()
    for index, label in enumerate(reference_labels):
        if label in vdw_atom_parameters:
            sigma, epsilon = vdw_atom_parameters[label]
            active_vdw = 1.0
            active_indices.add(index)
        elif label in vdw_universal_atoms:
            sigma, epsilon = vdw_universal_atoms[label]
            active_vdw = 0.0
        else:
            sigma = 0.0
            epsilon = 0.0
            active_vdw = 0.0
        if label in electro_atom_parameters:
            charge = electro_atom_parameters[label]
            active_elec = 1.0
            active_indices.add(index)
        elif label in electro_universal_atoms:
            charge = electro_universal_atoms[label]
            active_elec = 0.0
        else:
            charge = 0.0
            active_elec = 0.0
        qm_flag = 1.0 if label in qm_labels else 0.0
        force.addParticle((sigma, epsilon, charge, active_vdw, active_elec, qm_flag))
        if force_qm is not None:
            force_qm.addParticle((sigma, epsilon, charge, active_vdw, active_elec, qm_flag))
    return active_indices


def build_vdw_bond_terms(
    *,
    exception_parameters: Mapping[CanonicalPairKey, Tuple[float, float]],
    restored_pairs: Sequence[CanonicalPairKey],
    atom_parameters: Mapping[AtomReference, Tuple[float, float]],
    universal_atoms: Mapping[AtomReference, Tuple[float, float]],
) -> list[tuple[CanonicalPairKey, Tuple[float, float]]]:
    terms = list(sorted(exception_parameters.items()))
    for pair in restored_pairs:
        atom_a, atom_b = pair
        sigma_a, epsilon_a = atom_parameters.get(atom_a, universal_atoms.get(atom_a, (0.0, 0.0)))
        sigma_b, epsilon_b = atom_parameters.get(atom_b, universal_atoms.get(atom_b, (0.0, 0.0)))
        epsilon_mix = (epsilon_a * epsilon_b) ** 0.5
        if abs(epsilon_mix) <= 1e-12:
            continue
        sigma_mix = 0.5 * (sigma_a + sigma_b)
        terms.append((pair, (sigma_mix, epsilon_mix)))
    return terms


def build_electro_bond_terms(
    *,
    electro_exception_parameters: Mapping[CanonicalPairKey, float],
    restored_pairs: Sequence[CanonicalPairKey],
    atom_parameters: Mapping[AtomReference, float],
    universal_atoms: Mapping[AtomReference, float],
) -> list[tuple[CanonicalPairKey, float]]:
    terms = list(sorted(electro_exception_parameters.items()))
    for pair in restored_pairs:
        atom_a, atom_b = pair
        charge_a = atom_parameters.get(atom_a, universal_atoms.get(atom_a, 0.0))
        charge_b = atom_parameters.get(atom_b, universal_atoms.get(atom_b, 0.0))
        charge_prod = charge_a * charge_b
        if abs(charge_prod) <= 1e-12:
            continue
        terms.append((pair, charge_prod))
    return terms


def split_terms_by_qm_channel(
    terms: Iterable[tuple[CanonicalPairKey, _T]],
    *,
    pair_in_qm_channel: Callable[[CanonicalPairKey], bool],
) -> tuple[list[tuple[CanonicalPairKey, _T]], list[tuple[CanonicalPairKey, _T]]]:
    mm_terms: list[tuple[CanonicalPairKey, _T]] = []
    qm_terms: list[tuple[CanonicalPairKey, _T]] = []
    for pair, values in terms:
        if pair_in_qm_channel(pair):
            qm_terms.append((pair, values))
        else:
            mm_terms.append((pair, values))
    return mm_terms, qm_terms


def build_lambda_scaled_vdw_exception_force(
    *,
    lambda_param: str,
    terms: Sequence[tuple[CanonicalPairKey, Tuple[float, float]]],
    initial_state: str,
    particle_indices_from_canonical_pair: Callable[[str, CanonicalPairKey], tuple[int, int]],
    periodic: bool = False,
) -> openmm.CustomBondForce | None:
    if not terms:
        return None
    force = openmm.CustomBondForce(f"{lambda_param} * 4*epsilon * ( (sigma/r)^12 - (sigma/r)^6 )")
    force.setUsesPeriodicBoundaryConditions(periodic)
    force.addGlobalParameter(lambda_param, 0.0)
    force.addPerBondParameter("sigma")
    force.addPerBondParameter("epsilon")
    for pair, (sigma, epsilon) in terms:
        idx_a, idx_b = particle_indices_from_canonical_pair(initial_state, pair)
        force.addBond(idx_a, idx_b, (sigma, epsilon))
    return force


def build_lambda_scaled_electro_exception_force(
    *,
    lambda_param: str,
    terms: Sequence[tuple[CanonicalPairKey, float]],
    initial_state: str,
    particle_indices_from_canonical_pair: Callable[[str, CanonicalPairKey], tuple[int, int]],
    one_4pi_eps0: float,
    periodic: bool = False,
) -> openmm.CustomBondForce | None:
    if not terms:
        return None
    force = openmm.CustomBondForce(f"{lambda_param} * ONE_4PI_EPS0 * charge_prod / r")
    force.setUsesPeriodicBoundaryConditions(periodic)
    force.addGlobalParameter(lambda_param, 0.0)
    force.addGlobalParameter("ONE_4PI_EPS0", one_4pi_eps0)
    force.addPerBondParameter("charge_prod")
    for pair, charge_prod in terms:
        idx_a, idx_b = particle_indices_from_canonical_pair(initial_state, pair)
        force.addBond(idx_a, idx_b, (charge_prod,))
    return force


def restrict_to_subset_interaction_group(
    force: openmm.CustomNonbondedForce,
    subset_indices: Set[int],
) -> None:
    """Limit pair evaluation to pairs involving at least one atom in subset.

    Used wherever a CustomNonbondedForce's expression contains a mask of the
    form ``(1 - (1 - flag1)*(1 - flag2))`` (zero unless one of the two atoms
    is in the named subset). Without an interaction group OpenMM evaluates
    the expression for every pair within the cutoff and discards the zero
    result. An interaction group of ``{subset} × {all}`` makes the neighbor
    list itself skip pairs that lie entirely outside the subset, cutting
    pair work by roughly ``N_atoms / |subset|``.

    Applies to both the per-state vdW/electrostatics forces (``active``
    subset = state-specific atoms) and the shared QM channel (``qm`` subset
    = reactive atoms).
    """
    if not subset_indices:
        return
    all_indices = set(range(force.getNumParticles()))
    force.addInteractionGroup(subset_indices, all_indices)


__all__ = [
    "split_shared_exception_terms",
    "build_shared_exception_force",
    "collect_state_exclusion_pairs",
    "populate_state_combined_particles",
    "restrict_to_subset_interaction_group",
    "build_vdw_bond_terms",
    "build_electro_bond_terms",
    "split_terms_by_qm_channel",
    "build_lambda_scaled_vdw_exception_force",
    "build_lambda_scaled_electro_exception_force",
]
