from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence

import openmm
from openmm import app

from ..models import FEPSampleRecord, FEPWindowRecord


def collect_fep_window_records(
    *,
    simulation: app.Simulation,
    windows: Sequence[Sequence[float]],
    state_names: Sequence[str],
    equilibration_steps: int,
    production_steps: int,
    sample_interval: int,
    coefficients: Sequence[float],
    set_state_lambdas: Callable[[app.Simulation, Iterable[float] | float], Mapping[str, float]],
    sample_diabatic_state_energies_kj_mol: Callable[[app.Simulation], Sequence[float]],
    effective_energy_from_state_energies: Callable[[Sequence[float], Sequence[float]], float],
    reaction_coordinate_from_state_energies: Callable[[Sequence[float], Sequence[float]], float],
    sample_coupling_distance_nm: Callable[[app.Simulation], float] | None = None,
) -> tuple[FEPWindowRecord, ...]:
    records: list[FEPWindowRecord] = []
    for window_index, lambda_values in enumerate(windows):
        partial_values: Iterable[float] | float
        if len(lambda_values) == 2:
            partial_values = lambda_values[0]
        else:
            partial_values = lambda_values[:-1]
        applied_lambdas = set_state_lambdas(simulation, partial_values)
        lambda_sequence = tuple(float(applied_lambdas[state]) for state in state_names)

        if equilibration_steps > 0:
            try:
                simulation.step(equilibration_steps)
            except openmm.OpenMMException as exc:
                raise RuntimeError(
                    "Mapping-FEP instability during per-window equilibration "
                    f"(window_index={window_index}, lambda_values={lambda_sequence}, "
                    f"equilibration_steps={int(equilibration_steps)}): {exc}"
                ) from exc

        samples: list[FEPSampleRecord] = []
        window_step = 0
        remaining = int(production_steps)
        while remaining > 0:
            advance = min(int(sample_interval), remaining)
            try:
                simulation.step(advance)
            except openmm.OpenMMException as exc:
                failed_step = window_step + advance
                raise RuntimeError(
                    "Mapping-FEP instability during production "
                    f"(window_index={window_index}, lambda_values={lambda_sequence}, "
                    f"window_step={failed_step}, sample_interval={int(sample_interval)}): {exc}"
                ) from exc
            remaining -= advance
            window_step += advance

            if remaining > 0 and window_step % int(sample_interval) != 0:
                continue

            state_energies = tuple(float(value) for value in sample_diabatic_state_energies_kj_mol(simulation))
            effective_energy = float(effective_energy_from_state_energies(lambda_sequence, state_energies))
            reaction_coordinate = float(reaction_coordinate_from_state_energies(coefficients, state_energies))
            coupling_distance = (
                float(sample_coupling_distance_nm(simulation))
                if sample_coupling_distance_nm is not None
                else None
            )
            samples.append(
                FEPSampleRecord(
                    md_step=int(window_step),
                    effective_energy_kj_mol=effective_energy,
                    reaction_coordinate_kj_mol=reaction_coordinate,
                    state_energies_kj_mol=state_energies,
                    coupling_distance_nm=coupling_distance,
                )
            )

        records.append(
            FEPWindowRecord(
                index=int(window_index),
                lambda_values=lambda_sequence,
                samples=tuple(samples),
            )
        )
    return tuple(records)


__all__ = ["collect_fep_window_records"]
