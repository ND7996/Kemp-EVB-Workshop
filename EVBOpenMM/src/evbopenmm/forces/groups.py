from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Dict, Tuple

import openmm
from openmm import unit

from ..constants import REFERENCE_FORCE_CATEGORIES


def get_or_create_force_group(
    *,
    force_group_ids: Dict[tuple[str | None, str], int],
    force_group_definitions: Dict[int, tuple[str | None, str]],
    next_force_group_id: int,
    term: str,
    state: str | None,
) -> tuple[int, int]:
    key = (state, term)
    if key in force_group_ids:
        return force_group_ids[key], next_force_group_id
    if next_force_group_id >= 32:
        raise RuntimeError(
            "EVBSystem requires more than 32 OpenMM force groups. Reduce the number of tracked terms."
        )
    group = next_force_group_id
    force_group_ids[key] = group
    force_group_definitions[group] = (state, term)
    return group, next_force_group_id + 1


def collect_group_energies(
    *,
    context: openmm.Context,
    force_group_definitions: Mapping[int, tuple[str | None, str]],
) -> Dict[int, unit.Quantity]:
    energies: Dict[int, unit.Quantity] = {}
    for group in force_group_definitions:
        state = context.getState(getEnergy=True, groups=1 << group)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole) * unit.kilojoules_per_mole
        energies[group] = energy
    return energies


def context_parameter_exists(context: openmm.Context, parameter_name: str) -> bool:
    try:
        context.getParameter(parameter_name)
        return True
    except openmm.OpenMMException:
        return False


def context_get_parameter_or_default(
    context: openmm.Context,
    parameter_name: str,
    *,
    default: float,
) -> float:
    try:
        return float(context.getParameter(parameter_name))
    except openmm.OpenMMException:
        return float(default)


def context_lambda_snapshot(
    *,
    context: openmm.Context,
    state_names: Sequence[str],
    state_lambda_parameter_names: Mapping[str, str],
    context_get_parameter_or_default: Callable[[openmm.Context, str, float], float],
) -> Dict[str, float]:
    values: Dict[str, float] = {}
    for state_name in state_names:
        parameter_name = state_lambda_parameter_names.get(state_name) or f"lambda_{state_name}"
        values[state_name] = context_get_parameter_or_default(context, parameter_name, 0.0)
    return values


def apply_context_lambdas(
    *,
    context: openmm.Context,
    state_names: Sequence[str],
    state_lambda_parameter_names: Mapping[str, str],
    mapping: Mapping[str, float],
    context_parameter_exists: Callable[[openmm.Context, str], bool],
) -> None:
    for state_name in state_names:
        parameter_name = state_lambda_parameter_names.get(state_name) or f"lambda_{state_name}"
        value = float(mapping.get(state_name, 0.0))
        if context_parameter_exists(context, parameter_name):
            context.setParameter(parameter_name, value)


def apply_indicator_lambdas(
    *,
    context: openmm.Context,
    active_state: str,
    state_names: Sequence[str],
    apply_context_lambdas: Callable[[openmm.Context, Mapping[str, float]], None],
) -> None:
    if active_state not in state_names:
        raise ValueError(f"Unknown state '{active_state}' supplied for lambda toggling.")
    indicator = {state: (1.0 if state == active_state else 0.0) for state in state_names}
    apply_context_lambdas(context, indicator)


def assign_reference_force_groups(
    *,
    system: openmm.System,
    categorize_reference_force: Callable[[openmm.Force], str],
) -> Tuple[Dict[int, Dict[str, object]], Tuple[int, ...]]:
    metadata: Dict[int, Dict[str, object]] = {}
    previous_groups: list[int] = []
    for group_id in range(system.getNumForces()):
        if group_id >= 32:
            raise RuntimeError(
                "Reference system exposes more than 32 forces; unable to assign unique force groups for debugging."
            )
        force = system.getForce(group_id)
        previous_groups.append(force.getForceGroup())
        force.setForceGroup(group_id)
        metadata[group_id] = {
            "category": categorize_reference_force(force),
            "label": f"{force.__class__.__name__}[{group_id}]",
            "force": force,
        }
    return metadata, tuple(previous_groups)


def restore_reference_force_groups(system: openmm.System, previous_groups: Tuple[int, ...]) -> None:
    for idx, original_group in enumerate(previous_groups):
        force = system.getForce(idx)
        force.setForceGroup(original_group)


def categorize_reference_force(force: openmm.Force) -> str:
    for cls, category in REFERENCE_FORCE_CATEGORIES:
        if isinstance(force, cls):
            return category
    return "other"


def adiabatic_force_group_masks(
    *,
    state_names: Sequence[str],
    force_group_definitions: Mapping[int, Tuple[str | None, str]],
) -> Tuple[int, int, int]:
    """Classify registered force groups for direct adiabatic ground-state dynamics.

    Returns ``(universal_mask, state1_mask, state2_mask)`` as OpenMM force-group
    bitmasks. Forces registered with ``state=None`` (e.g. ``shared`` bonded
    terms, ``scaas_*`` restraints) belong to the universal mask; forces
    registered with ``state=<name>`` belong to that state's mask.

    The three masks are pairwise disjoint by construction — every entry in
    ``force_group_definitions`` is keyed by a unique ``(state, term)`` pair via
    :func:`get_or_create_force_group`, so each group bit appears in exactly one
    of the returned masks.
    """

    state_list = tuple(state_names)
    if len(state_list) != 2:
        raise ValueError(
            "Adiabatic ground-state dynamics currently require a two-state EVB "
            f"scaffold; received {len(state_list)} state(s): {state_list!r}."
        )
    state1, state2 = state_list

    masks: Dict[str | None, int] = {None: 0, state1: 0, state2: 0}
    for group_id, (state, term) in force_group_definitions.items():
        bit = 1 << int(group_id)
        if state not in masks:
            raise ValueError(
                f"Unexpected state {state!r} in force-group definition "
                f"(group={group_id}, term={term!r}); expected None or one of "
                f"{state_list!r}."
            )
        masks[state] |= bit
    return masks[None], masks[state1], masks[state2]


__all__ = [
    "get_or_create_force_group",
    "collect_group_energies",
    "context_parameter_exists",
    "context_get_parameter_or_default",
    "context_lambda_snapshot",
    "apply_context_lambdas",
    "apply_indicator_lambdas",
    "assign_reference_force_groups",
    "restore_reference_force_groups",
    "categorize_reference_force",
    "adiabatic_force_group_masks",
]
