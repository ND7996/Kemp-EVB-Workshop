from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence

import openmm
from openmm import unit

from ..constants import STATE_ENERGY_TERMS


def build_evb_energy_breakdown(
    *,
    force_group_ids: Mapping[tuple[str | None, str], int],
    group_energies: Mapping[int, unit.Quantity],
    state_names: Sequence[str],
) -> dict[str, object]:
    if not force_group_ids:
        raise RuntimeError("No force-group assignments detected. Add EVB forces before requesting energies.")

    zero_energy = 0.0 * unit.kilojoules_per_mole
    shared_total = zero_energy
    for (state_name, _term), shared_group in force_group_ids.items():
        if state_name is not None:
            continue
        shared_total = shared_total + group_energies.get(shared_group, zero_energy)

    states_breakdown: dict[str, dict[str, unit.Quantity]] = {}
    cumulative_state_energy = zero_energy
    for state_name in state_names:
        term_map: dict[str, unit.Quantity] = {}
        state_total = zero_energy
        for term in STATE_ENERGY_TERMS:
            if term == "nonbonded":
                group_terms = ("nonbonded", "nonbonded_qm")
            else:
                group_terms = (term,)

            energy = zero_energy
            for group_term in group_terms:
                group = force_group_ids.get((state_name, group_term))
                if group is None:
                    continue
                energy = energy + group_energies.get(group, zero_energy)
            term_map[term] = energy
            state_total = state_total + energy
        term_map["total"] = state_total
        states_breakdown[state_name] = term_map
        cumulative_state_energy = cumulative_state_energy + state_total

    grand_total = cumulative_state_energy + shared_total
    return {
        "states": states_breakdown,
        "shared_total": shared_total,
        "grand_total": grand_total,
    }


def energy_value_kj(value: object) -> float:
    """kJ/mol float from an OpenMM energy Quantity (or a bare number); 0.0 for None."""
    if value is None:
        return 0.0
    if unit.is_quantity(value):
        return float(value.value_in_unit(unit.kilojoule_per_mole))
    return float(value)


def accumulate_energy_breakdown(acc: object, new: object) -> object:
    """Recursively add two ``build_evb_energy_breakdown`` structures (numeric/Quantity leaves)."""
    if isinstance(acc, Mapping):
        return {key: accumulate_energy_breakdown(acc[key], new[key]) for key in acc}
    if unit.is_quantity(acc) or isinstance(acc, (int, float)):
        return acc + new
    return acc  # non-numeric leaf (labels etc.) -- keep as is


def scale_energy_breakdown(breakdown: object, factor: float) -> object:
    """Recursively scale the numeric/Quantity leaves of a breakdown structure."""
    if isinstance(breakdown, Mapping):
        return {key: scale_energy_breakdown(value, factor) for key, value in breakdown.items()}
    if unit.is_quantity(breakdown) or isinstance(breakdown, (int, float)):
        return breakdown * factor
    return breakdown


def split_nonbonded_energy(
    *,
    force: openmm.NonbondedForce,
    context: openmm.Context,
    group: int,
) -> tuple[unit.Quantity, unit.Quantity]:
    particle_parameters = [force.getParticleParameters(idx) for idx in range(force.getNumParticles())]
    exception_parameters = [force.getExceptionParameters(idx) for idx in range(force.getNumExceptions())]

    zero_charge = 0.0 * unit.elementary_charge
    zero_charge_prod = 0.0 * unit.elementary_charge**2
    zero_epsilon = 0.0 * unit.kilojoules_per_mole

    for idx, (_, sigma, epsilon) in enumerate(particle_parameters):
        force.setParticleParameters(idx, zero_charge, sigma, epsilon)
    for idx, (i, j, _, sigma, epsilon) in enumerate(exception_parameters):
        force.setExceptionParameters(idx, i, j, zero_charge_prod, sigma, epsilon)
    force.updateParametersInContext(context)
    vdw_state = context.getState(getEnergy=True, groups=1 << group)
    vdw_energy = vdw_state.getPotentialEnergy()

    for idx, (charge, sigma, _) in enumerate(particle_parameters):
        force.setParticleParameters(idx, charge, sigma, zero_epsilon)
    for idx, (i, j, charge_prod, sigma, _) in enumerate(exception_parameters):
        force.setExceptionParameters(idx, i, j, charge_prod, sigma, zero_epsilon)
    force.updateParametersInContext(context)
    electro_state = context.getState(getEnergy=True, groups=1 << group)
    electro_energy = electro_state.getPotentialEnergy()

    for idx, (charge, sigma, epsilon) in enumerate(particle_parameters):
        force.setParticleParameters(idx, charge, sigma, epsilon)
    for idx, (i, j, charge_prod, sigma, epsilon) in enumerate(exception_parameters):
        force.setExceptionParameters(idx, i, j, charge_prod, sigma, epsilon)
    force.updateParametersInContext(context)

    return vdw_energy, electro_energy


def reduce_reference_terms(
    *,
    context: openmm.Context,
    group_metadata: Mapping[int, Mapping[str, object]],
    group_energies: Mapping[int, unit.Quantity],
    split_nonbonded_energy: Callable[..., tuple[unit.Quantity, unit.Quantity]],
) -> tuple[dict[str, unit.Quantity], dict[str, unit.Quantity]]:
    zero_energy = 0.0 * unit.kilojoules_per_mole
    reference_terms = {term: zero_energy for term in STATE_ENERGY_TERMS}
    other_terms: dict[str, unit.Quantity] = {}

    nonbonded_payload: tuple[openmm.NonbondedForce, int, unit.Quantity] | None = None
    for group_id, metadata in group_metadata.items():
        energy = group_energies[group_id]
        category = metadata["category"]
        if category in {"bond", "angle", "torsion"}:
            reference_terms[category] = reference_terms[category] + energy
        elif category == "nonbonded":
            force = metadata["force"]
            if not isinstance(force, openmm.NonbondedForce):
                raise TypeError(
                    "Nonbonded category expected an openmm.NonbondedForce; "
                    f"received {force.__class__.__name__}."
                )
            if nonbonded_payload is not None:
                raise RuntimeError(
                    "Multiple NonbondedForce instances detected; extend the debug helper "
                    "before attempting to split their energies."
                )
            nonbonded_payload = (force, group_id, energy)
        else:
            other_terms[str(metadata["label"])] = energy

    if nonbonded_payload is not None:
        nonbonded_force, group_id, total_nb = nonbonded_payload
        reference_terms["nonbonded"] = total_nb
    else:
        reference_terms["nonbonded"] = zero_energy

    return reference_terms, other_terms


def print_energy_debug_summary(data: Mapping[str, object], *, state_names: Sequence[str]) -> None:
    zero_energy = 0.0 * unit.kilojoules_per_mole
    breakdown = data["evb_breakdown"]
    shared_total = breakdown["shared_total"]
    evb_state_total = data["evb_state_total"]
    evb_state_total_with_shared = data["evb_state_total_with_shared"]
    shared_term_estimates = data.get("shared_term_estimates", {})
    shared_term_spread = data.get("shared_term_spread", {})
    shared_term_estimate_residual = data.get("shared_term_estimate_residual", zero_energy)
    reference_total = data["reference_grand_total"]
    state_results = data["states"]
    kcal = lambda value: value.value_in_unit(unit.kilocalories_per_mole)

    print("=== EVB Debug Energy Summary ===")
    print("Totals (kcal/mol)")
    print(f"  Shared (once):                {kcal(shared_total): 12.6f}")
    print(f"  Σ EVB states (state-only):    {kcal(evb_state_total): 12.6f}")
    print(f"  Σ EVB states + shared:        {kcal(evb_state_total_with_shared): 12.6f}")
    print(f"  Σ reference states:           {kcal(reference_total): 12.6f}")
    if shared_term_estimates:
        print("\nShared decomposition estimates (kcal/mol)")
        print(f"  {'term':<15s}{'estimate':>14s}{'spread':>14s}")
        for term in STATE_ENERGY_TERMS:
            shared_value = kcal(shared_term_estimates.get(term, zero_energy))
            spread_value = float(shared_term_spread.get(term, 0.0))
            print(f"  {term:<15s}{shared_value: 14.6f}{spread_value: 14.6f}")
        residual_value = kcal(shared_term_estimate_residual)
        print(f"  {'residual':<15s}{residual_value: 14.6f}")

    for state_name in state_names:
        info = state_results[state_name]
        evb_terms = info["evb_terms"]
        evb_terms_with_shared = info.get("evb_terms_with_shared_estimate", {})
        ref_terms = info["reference_terms"]
        deltas = info["term_deltas"]
        corrected = info.get("term_deltas_with_shared_estimate", {})
        print(f"\nState '{state_name}' (indicator lambda)")
        print(
            f"  {'term':<15s}{'EVB(state)':>14s}{'shared(est)':>14s}{'EVB+shared':>14s}"
            f"{'reference':>14s}{'Δ(state)':>12s}{'Δ(with)':>12s}"
        )
        for term in STATE_ENERGY_TERMS:
            evb_val = kcal(evb_terms.get(term, zero_energy))
            evb_with_shared = kcal(evb_terms_with_shared.get(term, zero_energy))
            shared_est = evb_with_shared - evb_val
            ref_val = kcal(ref_terms[term])
            delta_val = kcal(deltas[term])
            corrected_val = kcal(corrected.get(term, zero_energy))
            print(
                f"  {term:<15s}{evb_val: 14.6f}{shared_est: 14.6f}{evb_with_shared: 14.6f}"
                f"{ref_val: 14.6f}{delta_val: 12.6f}{corrected_val: 12.6f}"
            )
        total_delta = kcal(info["total_delta"])
        print(f"  {'total Δ (Ref - [EVB+shared])':<57s}{total_delta: 12.6f}")


__all__ = [
    "build_evb_energy_breakdown",
    "print_energy_debug_summary",
    "reduce_reference_terms",
    "split_nonbonded_energy",
]
