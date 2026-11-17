from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from types import MappingProxyType

import openmm
from openmm import app, unit

from ..models import MDEnergySampleRecord


def run_endpoint_state_comparison(
    *,
    state_name: str,
    state_seed: int | None,
    evb_source_state: str,
    evb_topology: app.Topology,
    evb_system: openmm.System,
    evb_start_positions: Sequence[openmm.Vec3],
    expected_particles: int,
    reference_topology: app.Topology,
    reference_system: openmm.System,
    temperature: unit.Quantity,
    friction: unit.Quantity,
    timestep: unit.Quantity,
    platform: openmm.Platform | None,
    platform_properties: Mapping[str, str] | None,
    initialize_velocities: bool,
    equilibration_steps: int,
    production_steps: int,
    sample_interval: int,
    create_langevin_middle_simulation: Callable[..., app.Simulation],
    attach_scaas_to_simulation: Callable[[app.Simulation, bool], None],
    apply_indicator_lambdas: Callable[..., None],
    coerce_vec3_tuple: Callable[..., tuple[openmm.Vec3, ...]],
    vec3s_to_quantity: Callable[..., object],
    collect_md_energy_samples: Callable[..., tuple[MDEnergySampleRecord, ...]],
    get_state_lambdas: Callable[[app.Simulation], Mapping[str, float]],
    remap_positions_for_state: Callable[..., tuple[openmm.Vec3, ...]],
    summarize_md_energy_samples: Callable[[Sequence[MDEnergySampleRecord]], Mapping[str, float]],
) -> Mapping[str, object]:
    evb_simulation = create_langevin_middle_simulation(
        topology=evb_topology,
        system=evb_system,
        temperature=temperature,
        friction=friction,
        timestep=timestep,
        random_seed=state_seed,
        platform=platform,
        platform_properties=platform_properties,
    )
    evb_simulation.context.setPositions(vec3s_to_quantity(evb_start_positions))
    attach_scaas_to_simulation(evb_simulation, initialize=True)
    apply_indicator_lambdas(evb_simulation.context, active_state=state_name)

    evb_initial_velocities: tuple[openmm.Vec3, ...] | None = None
    if initialize_velocities:
        if state_seed is None:
            evb_simulation.context.setVelocitiesToTemperature(temperature)
        else:
            evb_simulation.context.setVelocitiesToTemperature(temperature, int(state_seed))
        evb_velocity_quantity = evb_simulation.context.getState(getVelocities=True).getVelocities()
        evb_initial_velocities = coerce_vec3_tuple(
            evb_velocity_quantity,
            quantity_unit=unit.nanometer / unit.picosecond,
            expected_count=expected_particles,
            label="EVB initial velocities",
        )

    evb_samples = collect_md_energy_samples(
        evb_simulation,
        equilibration_steps=equilibration_steps,
        production_steps=production_steps,
        sample_interval=sample_interval,
    )
    evb_lambda_values = dict(get_state_lambdas(evb_simulation))

    reference_simulation = create_langevin_middle_simulation(
        topology=reference_topology,
        system=reference_system,
        temperature=temperature,
        friction=friction,
        timestep=timestep,
        random_seed=state_seed,
        platform=platform,
        platform_properties=platform_properties,
    )

    try:
        reference_start_positions = remap_positions_for_state(
            source_state=evb_source_state,
            target_state=state_name,
            source_positions=evb_start_positions,
        )
    except ValueError as exc:
        raise ValueError(
            f"Failed to remap starting coordinates from '{evb_source_state}' to '{state_name}'. "
            "Ensure atom mapping is complete across states before comparing endpoint trajectories."
        ) from exc
    reference_simulation.context.setPositions(vec3s_to_quantity(reference_start_positions))

    if initialize_velocities:
        if evb_initial_velocities is None:
            raise RuntimeError("Internal error: EVB initial velocities were not generated.")
        try:
            reference_initial_velocities = remap_positions_for_state(
                source_state=evb_source_state,
                target_state=state_name,
                source_positions=evb_initial_velocities,
            )
        except ValueError as exc:
            raise ValueError(
                f"Failed to remap initial velocities from '{evb_source_state}' to '{state_name}'. "
                "Ensure atom mapping is complete across states before comparing endpoint trajectories."
            ) from exc
        reference_simulation.context.setVelocities(
            vec3s_to_quantity(
                reference_initial_velocities,
                quantity_unit=unit.nanometer / unit.picosecond,
            )
        )

    reference_samples = collect_md_energy_samples(
        reference_simulation,
        equilibration_steps=equilibration_steps,
        production_steps=production_steps,
        sample_interval=sample_interval,
    )

    evb_summary = summarize_md_energy_samples(evb_samples)
    reference_summary = summarize_md_energy_samples(reference_samples)
    mean_delta = {
        key: evb_summary[key] - reference_summary[key]
        for key in ("potential_energy_kj_mol", "kinetic_energy_kj_mol", "total_energy_kj_mol")
    }

    return {
        "seed": state_seed,
        "evb_endpoint": {
            "lambda_values": MappingProxyType(evb_lambda_values),
            "samples": evb_samples,
            "mean_energies_kj_mol": MappingProxyType(dict(evb_summary)),
        },
        "reference_ff": {
            "samples": reference_samples,
            "mean_energies_kj_mol": MappingProxyType(dict(reference_summary)),
        },
        "mean_energy_delta_kj_mol": MappingProxyType(mean_delta),
    }


__all__ = ["run_endpoint_state_comparison"]
