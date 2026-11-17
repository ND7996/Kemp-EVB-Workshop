from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from openmm import app

from ..types import PathLike


def setup_coordinate_reporting(
    *,
    simulation: app.Simulation,
    coordinate_report_interval_ps: float | None,
    coordinate_report_file: PathLike | None,
    coordinate_report_append: bool,
    coordinate_topology_file: PathLike | None,
    simulation_timestep_ps: Callable[[app.Simulation], float],
    steps_from_report_interval_ps: Callable[[float, float], int],
) -> tuple[object | None, int | None, Path | None, Path | None]:
    if coordinate_report_interval_ps is None:
        return None, None, None, None

    if coordinate_report_file is None:
        raise ValueError("Set 'coordinate_report_file' when using 'coordinate_report_interval_ps'.")

    interval_ps = float(coordinate_report_interval_ps)
    timestep_ps = simulation_timestep_ps(simulation)
    report_interval_steps = steps_from_report_interval_ps(interval_ps, timestep_ps)
    coordinate_path = Path(coordinate_report_file)
    coordinate_path.parent.mkdir(parents=True, exist_ok=True)
    reporter = app.DCDReporter(
        str(coordinate_path),
        reportInterval=report_interval_steps,
        append=bool(coordinate_report_append),
    )
    simulation.reporters.append(reporter)

    topology_path = (
        Path(coordinate_topology_file) if coordinate_topology_file is not None else coordinate_path.with_suffix(".pdb")
    )
    topology_path.parent.mkdir(parents=True, exist_ok=True)
    start_state = simulation.context.getState(getPositions=True)
    start_positions = start_state.getPositions()
    with topology_path.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(simulation.topology, start_positions, handle)

    return reporter, int(report_interval_steps), coordinate_path, topology_path


def remove_coordinate_reporter(simulation: app.Simulation, reporter: object | None) -> None:
    if reporter is None:
        return
    try:
        simulation.reporters.remove(reporter)
    except ValueError:
        pass


__all__ = ["setup_coordinate_reporting", "remove_coordinate_reporter"]
