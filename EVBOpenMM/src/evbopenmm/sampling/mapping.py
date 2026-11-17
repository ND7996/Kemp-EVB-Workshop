from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

from ..types import PathLike


def resolve_mapping_coordinate_outputs(
    *,
    coordinate_report_interval_ps: float | None,
    coordinate_report_file: PathLike | None,
    coordinate_topology_file: PathLike | None,
    output_directory: PathLike | None,
) -> tuple[PathLike | None, PathLike | None]:
    resolved_coordinate_file = coordinate_report_file
    if (
        coordinate_report_interval_ps is not None
        and resolved_coordinate_file is None
        and output_directory is not None
    ):
        resolved_coordinate_file = Path(output_directory) / "coordinates.dcd"

    resolved_coordinate_topology_file = coordinate_topology_file
    if (
        coordinate_report_interval_ps is not None
        and resolved_coordinate_topology_file is None
        and resolved_coordinate_file is not None
    ):
        resolved_coordinate_topology_file = Path(resolved_coordinate_file).with_suffix(".pdb")

    return resolved_coordinate_file, resolved_coordinate_topology_file


def run_mapping_window_workflow(
    *,
    coordinate_report_interval_ps: float | None,
    coordinate_report_file: PathLike | None,
    coordinate_topology_file: PathLike | None,
    output_directory: PathLike | None,
    collect_records: Callable[[PathLike | None, PathLike | None], Sequence[object]],
    analyze_records: Callable[[Sequence[object]], Mapping[str, object]],
    write_outputs: Sequence[Callable[[Mapping[str, object], PathLike], None]] = tuple(),
) -> dict[str, object]:
    resolved_coordinate_file, resolved_coordinate_topology_file = resolve_mapping_coordinate_outputs(
        coordinate_report_interval_ps=coordinate_report_interval_ps,
        coordinate_report_file=coordinate_report_file,
        coordinate_topology_file=coordinate_topology_file,
        output_directory=output_directory,
    )
    records = collect_records(resolved_coordinate_file, resolved_coordinate_topology_file)
    analysis = dict(analyze_records(records))
    if output_directory is not None:
        for write_output in write_outputs:
            write_output(analysis, output_directory)
    return analysis


__all__ = [
    "resolve_mapping_coordinate_outputs",
    "run_mapping_window_workflow",
]
