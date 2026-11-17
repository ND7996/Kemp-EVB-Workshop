from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from types import MappingProxyType
import csv
import math
import warnings

import openmm
from openmm import app
from openmm import unit

from ..types import PathLike
from .reporting import remove_coordinate_reporter, setup_coordinate_reporting


def validate_velocity_mode(mode: str) -> str:
    velocity_mode = str(mode).strip().lower()
    if velocity_mode not in {"first", "each", "never"}:
        raise ValueError("'initialize_velocities' must be one of: 'first', 'each', 'never'.")
    return velocity_mode


def build_restraint_schedule_kj_mol_nm2(
    *,
    n_stages: int,
    solute_restraint_k_schedule: Sequence[unit.Quantity | float] | None,
    solute_restraint_k_start: unit.Quantity | float,
    solute_restraint_k_end: unit.Quantity | float,
    coerce_force_constant_kj_mol_nm2: Callable[[unit.Quantity | float, str], unit.Quantity],
) -> tuple[float, ...]:
    if solute_restraint_k_schedule is None:
        k_start_q = coerce_force_constant_kj_mol_nm2(solute_restraint_k_start, "solute_restraint_k_start")
        k_end_q = coerce_force_constant_kj_mol_nm2(solute_restraint_k_end, "solute_restraint_k_end")
        k_start = float(k_start_q.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2))
        k_end = float(k_end_q.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2))
        if n_stages == 1:
            return (k_start,)
        return tuple(
            k_start + (k_end - k_start) * (float(stage_idx) / float(n_stages - 1))
            for stage_idx in range(n_stages)
        )

    raw_schedule = tuple(solute_restraint_k_schedule)
    if len(raw_schedule) != n_stages:
        raise ValueError(
            f"'solute_restraint_k_schedule' must contain {n_stages} entries, received {len(raw_schedule)}."
        )
    return tuple(
        float(
            coerce_force_constant_kj_mol_nm2(value, "solute_restraint_k_schedule").value_in_unit(
                unit.kilojoules_per_mole / unit.nanometer**2
            )
        )
        for value in raw_schedule
    )


def state_positions_nm(simulation):
    """Positions as a NumPy ``(N, 3)`` array in nm. Isolated so the monitor's one expensive
    operation has a single call site."""
    state = simulation.context.getState(getPositions=True)
    return state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)


def scan_positions_nm(positions_nm) -> tuple[float, int]:
    """(max |coordinate|, number of non-finite atoms) over a position array.

    Accepts a NumPy ``(N, 3)`` array (the fast path, used by the equilibration loop) or a sequence
    of ``Vec3`` (kept so existing callers and tests are unaffected). The two paths return
    identical values.

    Performance note, because this is called inside the per-step stability monitor: the Vec3 path
    costs three ``float()`` conversions and two ``math.isfinite()`` calls per atom, in pure Python.
    On an 88 715-atom system that is ~440 000 interpreter operations per MD step and made the
    staged heat-up ~140x slower than plain stepping (measured on an H100: 8.8 steps/s inside the
    monitor against 1240 outside it). The NumPy path does the same work vectorised.
    """
    array = getattr(positions_nm, "shape", None)
    if array is not None:                       # NumPy (N, 3) -- vectorised
        import numpy as np

        values = np.asarray(positions_nm, dtype=float)
        finite_rows = np.isfinite(values).all(axis=1)
        nonfinite_coords = int((~finite_rows).sum())
        max_abs = float(np.abs(values[finite_rows]).max()) if finite_rows.any() else 0.0
        return max_abs, nonfinite_coords

    max_abs = 0.0
    nonfinite_coords = 0
    for vec in positions_nm:
        x = float(vec.x)
        y = float(vec.y)
        z = float(vec.z)
        if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
            nonfinite_coords += 1
            continue
        max_abs = max(max_abs, abs(x), abs(y), abs(z))
    return float(max_abs), int(nonfinite_coords)


def detect_stage_instability(
    *,
    potential_kj_mol: float,
    kinetic_kj_mol: float,
    nonfinite_coords: int,
    max_abs_coord_nm: float,
    stage_index: int,
    stage_step: int,
    global_step: int,
    energy_abs_cap_kj_mol: float,
    max_abs_coord_nm_cap: float | None,
) -> str | None:
    if nonfinite_coords > 0:
        return (
            f"non-finite coordinates detected at stage {stage_index + 1} step {stage_step} "
            f"(global step {global_step})."
        )
    if not math.isfinite(potential_kj_mol) or not math.isfinite(kinetic_kj_mol):
        return (
            f"non-finite energy detected at stage {stage_index + 1} step {stage_step} "
            f"(global step {global_step})."
        )
    if abs(potential_kj_mol) > float(energy_abs_cap_kj_mol):
        return (
            f"|potential| exceeded {float(energy_abs_cap_kj_mol):.6g} kJ/mol at stage {stage_index + 1} "
            f"step {stage_step} (global step {global_step})."
        )
    if max_abs_coord_nm_cap is not None and max_abs_coord_nm > float(max_abs_coord_nm_cap):
        return (
            f"Max |coordinate| exceeded {float(max_abs_coord_nm_cap):.6g} nm at stage {stage_index + 1} "
            f"step {stage_step} (global step {global_step})."
        )
    return None


def build_stage_result(
    *,
    stage_index: int,
    stage_name: str,
    target_steps: int,
    completed_steps: int,
    dt_femtoseconds: float,
    target_temperature_kelvin: float,
    friction_per_picosecond: float,
    solute_restraint_k_kj_mol_nm2: float,
    lambda_values: dict[str, float],
    unstable: bool,
    message: str | None,
):
    return MappingProxyType(
        {
            "stage_index": int(stage_index),
            "stage_name": stage_name,
            "target_steps": int(target_steps),
            "completed_steps": int(completed_steps),
            "dt_femtoseconds": float(dt_femtoseconds),
            "target_temperature_kelvin": float(target_temperature_kelvin),
            "friction_per_picosecond": float(friction_per_picosecond),
            "solute_restraint_k_kj_mol_nm2": float(solute_restraint_k_kj_mol_nm2),
            "lambda_values": MappingProxyType(dict(lambda_values)),
            "unstable": bool(unstable),
            "message": message,
        }
    )


def run_equilibration_schedule(
    *,
    simulation: app.Simulation | None,
    evb_topology: app.Topology,
    evb_system_scaffold: openmm.System,
    evb_initial_positions: Sequence[openmm.Vec3],
    state_names: Sequence[str],
    step_counts: Sequence[int],
    timestep_schedule: Sequence[unit.Quantity],
    temperature_schedule: Sequence[unit.Quantity],
    stage_frictions: Sequence[unit.Quantity],
    stage_lambda_values: Sequence[Sequence[float]],
    velocity_mode: str,
    random_seed: int | None,
    platform: openmm.Platform | None,
    platform_properties: Mapping[str, str] | None,
    coordinate_report_interval_ps: float | None,
    coordinate_report_file: PathLike | None,
    coordinate_report_append: bool,
    coordinate_topology_file: PathLike | None,
    summary_file: PathLike | None,
    coordinate_csv_file: PathLike | None,
    coordinate_csv_interval_steps: int | None,
    energy_abs_cap_kj_mol: float,
    max_abs_coord_nm_cap: float | None,
    monitor_interval_steps: int = 1,
    solute_restraints: bool,
    solute_restraint_atom_indices: Sequence[int] | None,
    solute_restraint_k_schedule: Sequence[unit.Quantity | float] | None,
    solute_restraint_k_start: unit.Quantity | float,
    solute_restraint_k_end: unit.Quantity | float,
    solute_restraint_reference_positions: (
        unit.Quantity | Sequence[openmm.Vec3] | Sequence[Sequence[float]] | None
    ),
    solute_restraint_include_hydrogens: bool,
    solute_restraint_solvent_residue_names: Sequence[str] | None,
    solute_restraint_global_parameter: str,
    fail_on_instability: bool,
    create_langevin_middle_simulation: Callable[..., app.Simulation],
    attach_scaas_to_simulation: Callable[..., None],
    context_parameter_exists: Callable[[openmm.Context, str], bool],
    simulation_timestep_ps: Callable[[app.Simulation], float],
    steps_from_report_interval_ps: Callable[[float, float], int],
    estimate_kinetic_dof: Callable[[openmm.System], int],
    apply_lambda_stage_sequence: Callable[[app.Simulation, Sequence[float]], Mapping[str, float]],
    normalize_atom_index_sequence: Callable[..., tuple[int, ...]],
    default_solute_restraint_atom_indices: Callable[..., tuple[int, ...]],
    coerce_vec3_tuple: Callable[..., tuple[openmm.Vec3, ...]],
    coerce_force_constant_kj_mol_nm2: Callable[[unit.Quantity | float, str], unit.Quantity],
    register_force_to_group: Callable[..., None],
) -> Mapping[str, object]:
    local_simulation = simulation
    created_simulation = False
    if local_simulation is None:
        local_simulation = create_langevin_middle_simulation(
            topology=evb_topology,
            system=evb_system_scaffold,
            temperature=temperature_schedule[0],
            friction=stage_frictions[0],
            timestep=timestep_schedule[0],
            random_seed=random_seed,
            platform=platform,
            platform_properties=platform_properties,
        )
        local_simulation.context.setPositions(evb_initial_positions)
        attach_scaas_to_simulation(local_simulation, initialize=True)
        created_simulation = True
    else:
        if local_simulation.system.getNumParticles() != evb_system_scaffold.getNumParticles():
            raise ValueError(
                "Provided simulation uses a system with a different number of particles than this EVB scaffold."
            )
        attach_scaas_to_simulation(local_simulation, initialize=True)

    integrator = local_simulation.integrator
    if not hasattr(integrator, "setStepSize"):
        raise ValueError("Simulation integrator does not support staged timestep updates (missing 'setStepSize').")
    if not hasattr(integrator, "setTemperature"):
        raise ValueError("Simulation integrator does not support thermostat temperature updates (missing 'setTemperature').")
    supports_friction_updates = hasattr(integrator, "setFriction")
    if not supports_friction_updates:
        unique_friction = {round(value.value_in_unit(unit.picosecond**-1), 12) for value in stage_frictions}
        if len(unique_friction) > 1:
            raise ValueError(
                "Simulation integrator does not expose 'setFriction', but multiple stage friction values were provided."
            )

    if random_seed is not None and hasattr(integrator, "setRandomNumberSeed"):
        try:
            integrator.setRandomNumberSeed(int(random_seed))
        except openmm.OpenMMException:
            pass

    n_stages = len(step_counts)
    monitor_stride = int(monitor_interval_steps)
    if monitor_stride < 1:
        raise ValueError(
            f"'monitor_interval_steps' must be >= 1, received {monitor_interval_steps!r}."
        )

    restraint_force_index: int | None = None
    restraint_atom_indices: tuple[int, ...] = tuple()
    restraint_schedule_kj_mol_nm2: tuple[float, ...] = tuple()
    if solute_restraints:
        restraint_schedule_kj_mol_nm2 = build_restraint_schedule_kj_mol_nm2(
            n_stages=n_stages,
            solute_restraint_k_schedule=solute_restraint_k_schedule,
            solute_restraint_k_start=solute_restraint_k_start,
            solute_restraint_k_end=solute_restraint_k_end,
            coerce_force_constant_kj_mol_nm2=coerce_force_constant_kj_mol_nm2,
        )

        restraint_atom_indices = normalize_atom_index_sequence(
            solute_restraint_atom_indices,
            expected_particles=local_simulation.system.getNumParticles(),
        )
        if not restraint_atom_indices:
            restraint_atom_indices = default_solute_restraint_atom_indices(
                topology=local_simulation.topology,
                system=local_simulation.system,
                include_hydrogens=bool(solute_restraint_include_hydrogens),
                solvent_residue_names=solute_restraint_solvent_residue_names,
            )
        if not restraint_atom_indices:
            raise ValueError(
                "No atoms were selected for solute positional restraints. "
                "Provide 'solute_restraint_atom_indices' explicitly or adjust "
                "'solute_restraint_solvent_residue_names'."
            )

        if solute_restraint_reference_positions is None:
            reference_positions = coerce_vec3_tuple(
                local_simulation.context.getState(getPositions=True).getPositions(),
                quantity_unit=unit.nanometer,
                expected_count=local_simulation.system.getNumParticles(),
                label="simulation positions",
            )
        else:
            reference_positions = coerce_vec3_tuple(
                solute_restraint_reference_positions,
                quantity_unit=unit.nanometer,
                expected_count=local_simulation.system.getNumParticles(),
                label="solute_restraint_reference_positions",
            )

        restraint_force = openmm.CustomExternalForce(
            f"0.5*{solute_restraint_global_parameter}*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
        )
        restraint_force.addGlobalParameter(solute_restraint_global_parameter, float(restraint_schedule_kj_mol_nm2[0]))
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")
        for atom_index in restraint_atom_indices:
            reference = reference_positions[atom_index]
            restraint_force.addParticle(atom_index, (float(reference.x), float(reference.y), float(reference.z)))

        evb_system_scaffold.addForce(restraint_force)
        register_force_to_group(restraint_force, term="equilibration_restraint", state=None)
        restraint_force_index = evb_system_scaffold.getNumForces() - 1
        local_simulation.context.reinitialize(preserveState=True)
        attach_scaas_to_simulation(local_simulation, initialize=True)
        if random_seed is not None and hasattr(integrator, "setRandomNumberSeed"):
            try:
                integrator.setRandomNumberSeed(int(random_seed))
            except openmm.OpenMMException:
                pass
        if not context_parameter_exists(local_simulation.context, solute_restraint_global_parameter):
            raise RuntimeError(
                "Failed to initialize solute restraint context parameter "
                f"'{solute_restraint_global_parameter}'."
            )
        local_simulation.context.setParameter(solute_restraint_global_parameter, float(restraint_schedule_kj_mol_nm2[0]))

    summary_path = None if summary_file is None else Path(summary_file)
    if summary_path is not None:
        summary_path.parent.mkdir(parents=True, exist_ok=True)

    coords_path = None if coordinate_csv_file is None else Path(coordinate_csv_file)
    if coords_path is not None:
        coords_path.parent.mkdir(parents=True, exist_ok=True)
        if coordinate_csv_interval_steps is None:
            coordinate_csv_interval_steps = 1
    if coordinate_csv_interval_steps is not None and int(coordinate_csv_interval_steps) <= 0:
        raise ValueError("'coordinate_csv_interval_steps' must be > 0 when provided.")
    coordinate_interval = None if coordinate_csv_interval_steps is None else int(coordinate_csv_interval_steps)

    (
        added_coordinate_reporter,
        report_interval_steps,
        coordinate_path,
        topology_path,
    ) = setup_coordinate_reporting(
        simulation=local_simulation,
        coordinate_report_interval_ps=coordinate_report_interval_ps,
        coordinate_report_file=coordinate_report_file,
        coordinate_report_append=bool(coordinate_report_append),
        coordinate_topology_file=coordinate_topology_file,
        simulation_timestep_ps=simulation_timestep_ps,
        steps_from_report_interval_ps=steps_from_report_interval_ps,
    )

    dof = estimate_kinetic_dof(local_simulation.system)
    gas_constant = float(unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoules_per_mole / unit.kelvin))
    target_total_steps = int(sum(step_counts))
    target_total_time_ps = sum(
        int(stage_steps) * float(stage_dt.value_in_unit(unit.picoseconds))
        for stage_steps, stage_dt in zip(step_counts, timestep_schedule)
    )

    summary_handle = None
    summary_writer = None
    coords_handle = None
    coords_writer = None

    completed_steps = 0
    completed_time_ps = 0.0
    stage_results: list[Mapping[str, object]] = []
    instability_reason: str | None = None

    try:
        if summary_path is not None:
            summary_handle = summary_path.open("w", newline="", encoding="utf-8")
            summary_writer = csv.writer(summary_handle)
            summary_writer.writerow(
                [
                    "global_step",
                    "stage_index",
                    "stage_name",
                    "stage_step",
                    "dt_fs",
                    "target_temp_k",
                    "instant_temp_k",
                    "friction_per_ps",
                    "solute_restraint_k_kj_mol_nm2",
                    *[f"lambda_{state}" for state in state_names],
                    "potential_kj_mol",
                    "kinetic_kj_mol",
                    "total_kj_mol",
                    "max_abs_coord_nm",
                    "n_nonfinite_coords",
                ]
            )

        if coords_path is not None:
            coords_handle = coords_path.open("w", newline="", encoding="utf-8")
            coords_writer = csv.writer(coords_handle)
            coords_writer.writerow(
                ["global_step", "stage_index", "stage_name", "stage_step", "atom_index", "x_nm", "y_nm", "z_nm", "is_finite"]
            )

        for stage_index in range(n_stages):
            stage_name = f"stage_{stage_index + 1}"
            stage_steps = int(step_counts[stage_index])
            stage_dt = timestep_schedule[stage_index]
            stage_temperature = temperature_schedule[stage_index]
            stage_friction = stage_frictions[stage_index]
            stage_lambda = stage_lambda_values[stage_index]
            stage_restraint_k = float(restraint_schedule_kj_mol_nm2[stage_index]) if restraint_schedule_kj_mol_nm2 else 0.0

            integrator.setStepSize(stage_dt)
            integrator.setTemperature(stage_temperature)
            if supports_friction_updates:
                integrator.setFriction(stage_friction)
            if context_parameter_exists(local_simulation.context, solute_restraint_global_parameter):
                local_simulation.context.setParameter(solute_restraint_global_parameter, stage_restraint_k)

            applied_lambdas = apply_lambda_stage_sequence(local_simulation, stage_lambda)

            should_initialize = velocity_mode == "each" or (velocity_mode == "first" and stage_index == 0)
            if should_initialize:
                if random_seed is None:
                    local_simulation.context.setVelocitiesToTemperature(stage_temperature)
                else:
                    local_simulation.context.setVelocitiesToTemperature(stage_temperature, int(random_seed) + stage_index)

            stage_completed = 0
            stage_unstable = False
            stage_message: str | None = None
            stage_dt_fs = float(stage_dt.value_in_unit(unit.femtoseconds))
            stage_dt_ps = float(stage_dt.value_in_unit(unit.picoseconds))
            stage_target_temp = float(stage_temperature.value_in_unit(unit.kelvin))
            stage_friction_ps = float(stage_friction.value_in_unit(unit.picosecond**-1))

            # Step in blocks of `monitor_stride` and inspect the state once per block. At stride 1
            # this is the original per-step monitor. Above 1 it trades detection LATENCY (an
            # instability is caught within `monitor_stride` steps rather than immediately) for
            # throughput: each inspection is a GPU sync plus an O(N) position download, which
            # dominates the step itself on large systems, and stepping one at a time also prevents
            # the integrator from pipelining. Nothing about WHETHER an instability is caught
            # changes -- a NaN never heals -- and the run still stops before the next block.
            stage_step = 0
            while stage_step < stage_steps:
                block = min(monitor_stride, stage_steps - stage_step)
                local_simulation.step(block)
                stage_step += block
                completed_steps += block
                stage_completed += block
                completed_time_ps += stage_dt_ps * block

                # Energy first, positions only if something needs them. Downloading positions is
                # a GPU sync plus an O(N) transfer -- ~78 ms on an 88 715-atom system, which after
                # the scan was vectorised became the dominant remaining cost of the monitor. A
                # blow-up shows up in the energy (NaN or a cap breach) before it needs coordinates
                # to diagnose, so the common healthy path never pays for them.
                state = local_simulation.context.getState(getEnergy=True)
                potential = float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
                kinetic = float(state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole))
                total = potential + kinetic
                instant_temp = (2.0 * kinetic) / (float(dof) * gas_constant)

                energy_suspect = not (math.isfinite(potential) and math.isfinite(kinetic)) or (
                    abs(total) > float(energy_abs_cap_kj_mol)
                )
                needs_coords_now = (
                    coords_writer is not None
                    and coordinate_interval is not None
                    and completed_steps % coordinate_interval == 0
                )
                positions_nm = None
                max_abs, nonfinite_coords = 0.0, 0
                if energy_suspect or needs_coords_now or max_abs_coord_nm_cap is not None:
                    positions_nm = state_positions_nm(local_simulation)
                    max_abs, nonfinite_coords = scan_positions_nm(positions_nm)

                if summary_writer is not None:
                    summary_writer.writerow(
                        [
                            completed_steps,
                            stage_index,
                            stage_name,
                            stage_step,
                            stage_dt_fs,
                            stage_target_temp,
                            instant_temp,
                            stage_friction_ps,
                            stage_restraint_k,
                            *[applied_lambdas[state_name] for state_name in state_names],
                            potential,
                            kinetic,
                            total,
                            max_abs,
                            nonfinite_coords,
                        ]
                    )

                should_write_coords = needs_coords_now

                instability_local = detect_stage_instability(
                    potential_kj_mol=potential,
                    kinetic_kj_mol=kinetic,
                    nonfinite_coords=nonfinite_coords,
                    max_abs_coord_nm=max_abs,
                    stage_index=stage_index,
                    stage_step=stage_step,
                    global_step=completed_steps,
                    energy_abs_cap_kj_mol=float(energy_abs_cap_kj_mol),
                    max_abs_coord_nm_cap=None if max_abs_coord_nm_cap is None else float(max_abs_coord_nm_cap),
                )

                if coords_writer is not None and (should_write_coords or instability_local is not None):
                    if positions_nm is None:            # instability path: fetch for diagnosis
                        positions_nm = state_positions_nm(local_simulation)
                    for atom_index, vec in enumerate(positions_nm):
                        # positions are a NumPy (N, 3) array now; keep the Vec3 attribute access
                        # working for any caller that still passes one.
                        if hasattr(vec, "x"):
                            x, y, z = float(vec.x), float(vec.y), float(vec.z)
                        else:
                            x, y, z = float(vec[0]), float(vec[1]), float(vec[2])
                        is_finite = int(math.isfinite(x) and math.isfinite(y) and math.isfinite(z))
                        coords_writer.writerow(
                            [completed_steps, stage_index, stage_name, stage_step, atom_index, x, y, z, is_finite]
                        )

                if instability_local is not None:
                    stage_unstable = True
                    stage_message = instability_local
                    instability_reason = instability_local
                    break

            stage_results.append(
                build_stage_result(
                    stage_index=stage_index,
                    stage_name=stage_name,
                    target_steps=stage_steps,
                    completed_steps=stage_completed,
                    dt_femtoseconds=stage_dt_fs,
                    target_temperature_kelvin=stage_target_temp,
                    friction_per_picosecond=stage_friction_ps,
                    solute_restraint_k_kj_mol_nm2=float(stage_restraint_k),
                    lambda_values=dict(applied_lambdas),
                    unstable=stage_unstable,
                    message=stage_message,
                )
            )

            if stage_unstable:
                break
    finally:
        if summary_handle is not None:
            summary_handle.close()
        if coords_handle is not None:
            coords_handle.close()
        remove_coordinate_reporter(local_simulation, added_coordinate_reporter)
        if context_parameter_exists(local_simulation.context, solute_restraint_global_parameter):
            try:
                local_simulation.context.setParameter(solute_restraint_global_parameter, 0.0)
            except openmm.OpenMMException as exc:
                warnings.warn(
                    f"Failed to zero '{solute_restraint_global_parameter}' after equilibration: {exc}",
                    RuntimeWarning,
                    stacklevel=2,
                )
        if restraint_force_index is not None:
            try:
                evb_system_scaffold.removeForce(int(restraint_force_index))
                local_simulation.context.reinitialize(preserveState=True)
                attach_scaas_to_simulation(local_simulation, initialize=True)
            except openmm.OpenMMException as exc:
                warnings.warn(
                    "Failed to remove temporary equilibration restraint force cleanly: "
                    f"{exc}",
                    RuntimeWarning,
                    stacklevel=2,
                )

    final_state = local_simulation.context.getState(getEnergy=True)
    final_potential = float(final_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
    final_kinetic = float(final_state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole))
    final_total = final_potential + final_kinetic

    result: dict[str, object] = {
        "created_simulation": bool(created_simulation),
        "target_steps": int(target_total_steps),
        "completed_steps": int(completed_steps),
        "target_time_picoseconds": float(target_total_time_ps),
        "completed_time_picoseconds": float(completed_time_ps),
        "stages": tuple(stage_results),
        "unstable": bool(instability_reason is not None),
        "instability_reason": instability_reason,
        "final_potential_energy_kj_mol": float(final_potential),
        "final_kinetic_energy_kj_mol": float(final_kinetic),
        "final_total_energy_kj_mol": float(final_total),
        "summary_file": None if summary_path is None else str(summary_path),
        "coordinate_csv_file": None if coords_path is None else str(coords_path),
        "coordinate_report_file": None if coordinate_path is None else str(coordinate_path),
        "coordinate_topology_file": None if topology_path is None else str(topology_path),
        "coordinate_report_interval_steps": None if report_interval_steps is None else int(report_interval_steps),
        "solute_restraints": bool(solute_restraints),
        "solute_restraint_atom_indices": tuple(int(idx) for idx in restraint_atom_indices),
        "solute_restraint_schedule_kj_mol_nm2": tuple(float(value) for value in restraint_schedule_kj_mol_nm2),
    }

    if instability_reason is not None and bool(fail_on_instability):
        raise RuntimeError(instability_reason)
    return MappingProxyType(result)


def equilibrate_reactive_evb(
    evb_system: object,
    simulation: app.Simulation,
    *,
    lambda_state: object,
    temperatures_k: Sequence[float] = (50.0, 150.0, 300.0),
    ps_per_stage: Sequence[float] = (2.0, 2.0, 10.0),
    timestep_fs: float = 1.0,
    minimize_iterations: int = 5000,
    random_seed: int = 1,
    nan_retry_timestep_fs: float | None = 0.5,
    max_nan_retries: int = 1,
    chunk_ps: float = 1.0,
) -> dict:
    """Robustly equilibrate a reactive EVB ``Simulation`` at a fixed mapping ``lambda_state``.

    A one-call convenience wrapping the pattern hand-rolled in every reactive-EVB runner:
    minimize, then a graded heat-up (e.g. 50 -> 150 -> 300 K) while holding the mapping lambda
    fixed -- with an automatic NaN guard. The stiff reactive Morse coordinate can blow up at
    intermediate lambda / aggressive timesteps; if non-finite coordinates appear during the
    heat-up, the helper restores the minimized restart point, drops to ``nan_retry_timestep_fs``
    and restarts the whole heat-up (up to ``max_nan_retries`` times) instead of crashing.

    ``lambda_state`` is passed verbatim to ``evb_system.set_state_lambdas(simulation, lambda_state)``
    (a float, or a per-state sequence). Returns a dict::

        {"stable": True, "final_timestep_fs": <fs used>, "nan_retries_used": int,
         "stages_completed": int, "lambda_state": <as given>}

    Raises ``RuntimeError`` if the heat-up cannot be stabilized within the retry budget.
    """
    temps = [float(t) for t in temperatures_k]
    stage_ps = [float(p) for p in ps_per_stage]
    if len(temps) != len(stage_ps):
        raise ValueError(
            f"temperatures_k ({len(temps)}) and ps_per_stage ({len(stage_ps)}) must have equal length."
        )
    integrator = simulation.integrator
    if not hasattr(integrator, "setStepSize"):
        raise ValueError("Simulation integrator does not support timestep updates (missing 'setStepSize').")

    # minimize once at the target lambda; this is the restart point for every retry.
    evb_system.set_state_lambdas(simulation, lambda_state)
    if minimize_iterations and int(minimize_iterations) > 0:
        simulation.minimizeEnergy(maxIterations=int(minimize_iterations))
    else:
        simulation.minimizeEnergy()
    saved_positions = simulation.context.getState(getPositions=True).getPositions()

    current_dt_fs = float(timestep_fs)
    retries_used = 0
    while True:
        integrator.setStepSize(current_dt_fs * unit.femtosecond)
        simulation.context.setPositions(saved_positions)
        evb_system.set_state_lambdas(simulation, lambda_state)
        stages_completed = 0
        instability: tuple[int, int, float] | None = None
        for stage_index, (temp_k, this_ps) in enumerate(zip(temps, stage_ps)):
            simulation.context.setVelocitiesToTemperature(temp_k * unit.kelvin, int(random_seed) + stage_index)
            if hasattr(integrator, "setTemperature"):
                integrator.setTemperature(temp_k * unit.kelvin)
            evb_system.set_state_lambdas(simulation, lambda_state)
            total_steps = max(1, int(round(this_ps * 1000.0 / current_dt_fs)))
            chunk_steps = max(1, int(round(float(chunk_ps) * 1000.0 / current_dt_fs)))
            done = 0
            while done < total_steps:
                n = min(chunk_steps, total_steps - done)
                simulation.step(n)
                done += n
                positions = simulation.context.getState(getPositions=True).getPositions()
                _, nonfinite = scan_positions_nm(positions)
                if nonfinite > 0:
                    instability = (stage_index, done, current_dt_fs)
                    break
            if instability is not None:
                break
            stages_completed += 1

        if instability is None:
            return {
                "stable": True,
                "final_timestep_fs": current_dt_fs,
                "nan_retries_used": retries_used,
                "stages_completed": stages_completed,
                "lambda_state": lambda_state,
            }

        stage_index, step_in_stage, dt_fs = instability
        can_retry = (
            retries_used < int(max_nan_retries)
            and nan_retry_timestep_fs is not None
            and float(nan_retry_timestep_fs) < current_dt_fs
        )
        if not can_retry:
            reason = "retries exhausted" if retries_used >= int(max_nan_retries) else "no smaller timestep available"
            raise RuntimeError(
                f"equilibrate_reactive_evb: non-finite coordinates at stage {stage_index + 1} "
                f"step {step_in_stage} with timestep {dt_fs:g} fs ({reason})."
            )
        retries_used += 1
        current_dt_fs = float(nan_retry_timestep_fs)
        # loop restarts the entire graded heat-up from the minimized point at the smaller timestep.


__all__ = [
    "validate_velocity_mode",
    "build_restraint_schedule_kj_mol_nm2",
    "scan_positions_nm",
    "detect_stage_instability",
    "build_stage_result",
    "run_equilibration_schedule",
    "equilibrate_reactive_evb",
]
