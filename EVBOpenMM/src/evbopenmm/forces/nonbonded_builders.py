from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Dict, Iterable, Set, Tuple

import openmm
from openmm import unit

from ..models import ElectrostaticsPartition, VanDerWaalsPartition
from ..types import CanonicalPairKey


def collect_exception_parameters_across_states(
    *,
    reference_systems: Mapping[str, openmm.System],
    find_nonbonded_force: Callable[[openmm.System], openmm.NonbondedForce],
    extract_vdw_exception_parameters: Callable[
        [str, openmm.NonbondedForce],
        Tuple[Dict[CanonicalPairKey, Tuple[float, float]], Tuple[CanonicalPairKey, ...]],
    ],
) -> Dict[CanonicalPairKey, Tuple[float, float]]:
    all_exception_parameters: Dict[CanonicalPairKey, Tuple[float, float]] = {}
    for state_name, state_system in reference_systems.items():
        nb_force = find_nonbonded_force(state_system)
        params, _ = extract_vdw_exception_parameters(state_name, nb_force)
        for pair, values in params.items():
            all_exception_parameters.setdefault(pair, values)
    return all_exception_parameters


def collect_global_exclusion_pairs(
    *,
    partition: VanDerWaalsPartition,
    electrostatics_partition: ElectrostaticsPartition,
) -> Tuple[CanonicalPairKey, ...]:
    global_exclusion_pairs: Set[CanonicalPairKey] = set(partition.universal_zero_exclusions)
    global_exclusion_pairs.update(electrostatics_partition.universal_zero_exclusions)
    for zero_set in partition.state_specific_zero_exclusions.values():
        global_exclusion_pairs.update(zero_set)
    for zero_set in electrostatics_partition.state_specific_zero_exclusions.values():
        global_exclusion_pairs.update(zero_set)
    global_exclusion_pairs.update(partition.universal_exceptions.keys())
    global_exclusion_pairs.update(electrostatics_partition.universal_exceptions.keys())
    for exception_map in partition.state_specific_exceptions.values():
        global_exclusion_pairs.update(exception_map.keys())
    for exception_map in electrostatics_partition.state_specific_exceptions.values():
        global_exclusion_pairs.update(exception_map.keys())
    return tuple(sorted(global_exclusion_pairs))


def pair_in_qm_channel(pair: CanonicalPairKey, *, qm_labels: set[object]) -> bool:
    return pair[0] in qm_labels or pair[1] in qm_labels


def lj_sigma_combining_expr(lj_combining: str = "lorentz-berthelot") -> str:
    """Return the OpenMM CustomNonbondedForce sub-expression for the combined sigma.

    ``"lorentz-berthelot"`` (default, Amber/CHARMM) -> arithmetic ``0.5*(sigma1+sigma2)``;
    ``"geometric"`` (OPLS-AA/M) -> ``sqrt(sigma1*sigma2)``.  Epsilon is geometric in both.
    """
    rule = lj_combining.lower()
    if rule in ("geometric", "opls"):
        return "sqrt(sigma1*sigma2)"
    if rule in ("lorentz-berthelot", "lb", "arithmetic", "amber", "charmm"):
        return "(0.5*(sigma1 + sigma2))"
    raise ValueError(f"Unknown lj_combining {lj_combining!r} "
                     "(use 'lorentz-berthelot' or 'geometric').")


def build_channel_force(
    *,
    expression: str,
    cutoff_distance: unit.Quantity,
    one_4pi_eps0: float,
    periodic: bool = False,
) -> openmm.CustomNonbondedForce:
    force = openmm.CustomNonbondedForce(expression)
    force.setNonbondedMethod(
        openmm.CustomNonbondedForce.CutoffPeriodic
        if periodic
        else openmm.CustomNonbondedForce.CutoffNonPeriodic
    )
    force.setCutoffDistance(cutoff_distance)
    force.setUseLongRangeCorrection(False)
    force.addGlobalParameter("ONE_4PI_EPS0", one_4pi_eps0)
    force.addPerParticleParameter("sigma")
    force.addPerParticleParameter("epsilon")
    force.addPerParticleParameter("charge")
    force.addPerParticleParameter("qm")
    return force


def build_state_combined_force(
    *,
    lambda_param: str,
    cutoff_distance: unit.Quantity,
    qm_channel: bool | None,
    one_4pi_eps0: float,
    periodic: bool = False,
    include_coulomb: bool = True,
    lj_combining: str = "lorentz-berthelot",
) -> openmm.CustomNonbondedForce:
    """Per-state vdW + electrostatics overlay in a single CustomNonbondedForce.

    Replaces the historical pair of separate forces (vdW + electrostatics)
    with one kernel that evaluates both per-pair contributions inside the
    same neighbor-list pass. Two independent per-particle ``active_*``
    flags let each term mask itself: a pair contributes vdW only when at
    least one atom is in the vdW-active set, and Coulomb only when at
    least one atom is in the electrostatics-active set, so systems where
    the two state-specific subsets differ stay correct.

    ``qm_channel`` controls the QM-channel mask:
    - ``False``: MM-channel — multiply by ``(1 - qm1)*(1 - qm2)`` (pairs
      where neither atom is reactive).
    - ``True``:  QM-channel — multiply by ``(1 - (1 - qm1)*(1 - qm2))``
      (pairs where at least one atom is reactive).
    - ``None``:  no channel mask — evaluate every pair within cutoff. Use
      this when MM and QM cutoffs are unified so the two complementary
      channels can collapse into a single force.
    """
    if qm_channel is None:
        channel_prefix = ""
    else:
        channel = (
            "(1 - (1 - qm1)*(1 - qm2))" if qm_channel else "((1 - qm1)*(1 - qm2))"
        )
        channel_prefix = f"{channel} * "
    _sig = lj_sigma_combining_expr(lj_combining)
    vdw_term = (
        "(1 - (1 - active_vdw1)*(1 - active_vdw2)) * "
        "4*sqrt(epsilon1*epsilon2) * "
        f"(({_sig}/r)^12 - ({_sig}/r)^6)"
    )
    coulomb_term = (
        "(1 - (1 - active_elec1)*(1 - active_elec2)) * ONE_4PI_EPS0 * "
        "charge1 * charge2 / r"
    )
    if include_coulomb:
        body = f"{vdw_term} + {coulomb_term}"
    else:
        # Periodic PME mode: electrostatics is handled by dedicated PME
        # NonbondedForce objects, so this custom force carries van der Waals only.
        body = vdw_term
    expression = f"{lambda_param} * {channel_prefix}({body})"
    force = openmm.CustomNonbondedForce(expression)
    force.setNonbondedMethod(
        openmm.CustomNonbondedForce.CutoffPeriodic
        if periodic
        else openmm.CustomNonbondedForce.CutoffNonPeriodic
    )
    force.setCutoffDistance(cutoff_distance)
    force.setUseLongRangeCorrection(False)
    force.addGlobalParameter(lambda_param, 0.0)
    force.addGlobalParameter("ONE_4PI_EPS0", one_4pi_eps0)
    force.addPerParticleParameter("sigma")
    force.addPerParticleParameter("epsilon")
    # The per-particle parameter schema is identical whether or not the Coulomb
    # term is present so callers (e.g. populate_state_combined_particles) can add
    # the same 6-tuple. In periodic PME mode ``charge`` is simply unreferenced.
    force.addPerParticleParameter("charge")
    force.addPerParticleParameter("active_vdw")
    force.addPerParticleParameter("active_elec")
    force.addPerParticleParameter("qm")
    return force


__all__ = [
    "collect_exception_parameters_across_states",
    "collect_global_exclusion_pairs",
    "pair_in_qm_channel",
    "lj_sigma_combining_expr",
    "build_channel_force",
    "build_state_combined_force",
]
