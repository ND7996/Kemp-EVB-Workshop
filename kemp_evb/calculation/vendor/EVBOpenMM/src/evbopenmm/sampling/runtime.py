from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from numbers import Real
from types import MappingProxyType

import openmm
from openmm import app, unit

from ..models import MDEnergySampleRecord


def simulation_timestep_ps(
    simulation: app.Simulation,
    *,
    fallback_timestep: unit.Quantity | float | None = None,
    coerce_timestep: Callable[[unit.Quantity | float], unit.Quantity] | None = None,
) -> float:
    step_size: unit.Quantity | None = None
    integrator = simulation.integrator
    if hasattr(integrator, "getStepSize"):
        try:
            step_size = integrator.getStepSize()
        except openmm.OpenMMException:
            step_size = None
    if step_size is None and fallback_timestep is not None:
        if hasattr(fallback_timestep, "value_in_unit"):
            step_size = fallback_timestep
        elif coerce_timestep is not None:
            step_size = coerce_timestep(fallback_timestep)
    if step_size is None:
        raise ValueError(
            "Unable to determine integrator timestep for coordinate reporting. "
            "Provide a simulation integrator with getStepSize() or pass a valid 'timestep'."
        )
    return float(step_size.value_in_unit(unit.picoseconds))


def steps_from_report_interval_ps(interval_ps: float, timestep_ps: float) -> int:
    if interval_ps <= 0.0:
        raise ValueError("'coordinate_report_interval_ps' must be > 0.")
    if timestep_ps <= 0.0:
        raise ValueError("Integrator timestep must be > 0 ps.")
    raw_steps = interval_ps / timestep_ps
    return max(1, int(round(raw_steps)))


def collect_md_energy_samples(
    simulation: app.Simulation,
    *,
    equilibration_steps: int,
    production_steps: int,
    sample_interval: int,
) -> tuple[MDEnergySampleRecord, ...]:
    if equilibration_steps < 0:
        raise ValueError("'equilibration_steps' must be >= 0.")
    if production_steps <= 0:
        raise ValueError("'production_steps' must be > 0.")
    if sample_interval <= 0:
        raise ValueError("'sample_interval' must be > 0.")

    if equilibration_steps > 0:
        simulation.step(int(equilibration_steps))

    samples: list[MDEnergySampleRecord] = []
    steps_done = 0
    remaining = int(production_steps)
    while remaining > 0:
        advance = min(int(sample_interval), remaining)
        simulation.step(advance)
        remaining -= advance
        steps_done += advance

        if remaining > 0 and steps_done % int(sample_interval) != 0:
            continue

        state = simulation.context.getState(getEnergy=True)
        potential = float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        kinetic = float(state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole))
        samples.append(
            MDEnergySampleRecord(
                md_step=steps_done,
                potential_energy_kj_mol=potential,
                kinetic_energy_kj_mol=kinetic,
                total_energy_kj_mol=potential + kinetic,
            )
        )
    return tuple(samples)


def summarize_md_energy_samples(samples: Sequence[MDEnergySampleRecord]) -> dict[str, float]:
    if not samples:
        raise ValueError("Cannot summarize an empty set of MD samples.")
    weight = 1.0 / float(len(samples))
    potential = sum(sample.potential_energy_kj_mol for sample in samples) * weight
    kinetic = sum(sample.kinetic_energy_kj_mol for sample in samples) * weight
    total = sum(sample.total_energy_kj_mol for sample in samples) * weight
    return {
        "potential_energy_kj_mol": potential,
        "kinetic_energy_kj_mol": kinetic,
        "total_energy_kj_mol": total,
    }


def default_equilibration_schedule() -> tuple[tuple[float, ...], tuple[int, ...], tuple[float, ...]]:
    heating_temperatures = tuple(float(value) for value in range(50, 301, 25))
    timesteps_fs = (0.5,) + tuple(1.0 for _ in heating_temperatures)
    steps_per_stage = (100,) + tuple(100 for _ in heating_temperatures)
    temperatures_k = (25.0,) + heating_temperatures
    return timesteps_fs, steps_per_stage, temperatures_k


def coerce_lambda_stage_sequence(
    lambda_values: Iterable[float] | float | Mapping[str, float],
    *,
    state_names: Sequence[str],
    coerce_lambda_sequence: Callable[[Iterable[float] | float], tuple[float, ...]],
    normalize_full_lambda_sequence: Callable[[Sequence[float]], tuple[float, ...]],
) -> tuple[float, ...]:
    n_states = len(state_names)
    if isinstance(lambda_values, Mapping):
        if all(state in lambda_values for state in state_names):
            full = tuple(float(lambda_values[state]) for state in state_names)
            return normalize_full_lambda_sequence(full)

        missing = [state for state in state_names[:-1] if state not in lambda_values]
        if missing:
            joined = ", ".join(missing)
            raise ValueError(
                "Lambda mapping entries must include either all states or at least "
                f"states {tuple(state_names[:-1])}; missing: {joined}."
            )
        partial = tuple(float(lambda_values[state]) for state in state_names[:-1])
        return coerce_lambda_sequence(partial)

    if isinstance(lambda_values, Real):
        return coerce_lambda_sequence(float(lambda_values))

    sequence = tuple(float(value) for value in lambda_values)
    if len(sequence) == n_states:
        return normalize_full_lambda_sequence(sequence)
    return coerce_lambda_sequence(sequence)


def apply_lambda_stage_sequence(
    simulation: app.Simulation,
    lambda_sequence: Sequence[float],
    *,
    state_names: Sequence[str],
    apply_context_lambdas: Callable[[object, Mapping[str, float]], None],
) -> Mapping[str, float]:
    if len(lambda_sequence) != len(state_names):
        raise ValueError(
            f"Expected {len(state_names)} full-state lambda values, received {len(lambda_sequence)}."
        )
    applied = {state: float(lambda_sequence[idx]) for idx, state in enumerate(state_names)}
    apply_context_lambdas(simulation.context, applied)
    return MappingProxyType(applied)


def estimate_kinetic_dof(system: openmm.System) -> int:
    dof = 0
    for atom_index in range(system.getNumParticles()):
        mass_da = float(system.getParticleMass(atom_index).value_in_unit(unit.dalton))
        if mass_da > 0.0:
            dof += 3
    dof -= int(system.getNumConstraints())
    for force_index in range(system.getNumForces()):
        if isinstance(system.getForce(force_index), openmm.CMMotionRemover):
            dof -= 3
            break
    return max(1, int(dof))


__all__ = [
    "simulation_timestep_ps",
    "steps_from_report_interval_ps",
    "collect_md_energy_samples",
    "summarize_md_energy_samples",
    "default_equilibration_schedule",
    "coerce_lambda_stage_sequence",
    "apply_lambda_stage_sequence",
    "estimate_kinetic_dof",
]
