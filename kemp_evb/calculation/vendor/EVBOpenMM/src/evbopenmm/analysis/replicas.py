from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

from ..models import FEPWindowRecord
from ..types import PathLike


def analyze_evb_profile_replicas_result(
    replica_window_records: Sequence[Sequence[FEPWindowRecord]],
    *,
    analyze_profile: Callable[[Sequence[FEPWindowRecord]], Mapping[str, object]],
    aggregate_profiles: Callable[..., Mapping[str, object]],
    minimum_replicas_per_bin: int = 1,
    output_directory: PathLike | None = None,
    write_mapping_fep_outputs: Callable[[Mapping[str, object], PathLike], None] | None = None,
    write_evb_profile_outputs: Callable[[Mapping[str, object], PathLike], None] | None = None,
    write_replica_outputs: Callable[[Mapping[str, object], PathLike], None] | None = None,
) -> dict[str, object]:
    if not replica_window_records:
        raise ValueError("At least one replica window-record sequence is required.")

    replica_analyses = tuple(analyze_profile(records) for records in replica_window_records)
    aggregate = aggregate_profiles(
        replica_analyses,
        minimum_replicas_per_bin=minimum_replicas_per_bin,
    )
    result = {
        "replicas": replica_analyses,
        "aggregate": aggregate,
    }

    if output_directory is not None:
        base = Path(output_directory)
        base.mkdir(parents=True, exist_ok=True)
        for idx, analysis in enumerate(replica_analyses):
            replica_out = base / f"rep_{idx:03d}"
            if write_mapping_fep_outputs is not None:
                write_mapping_fep_outputs(analysis, replica_out)
            if write_evb_profile_outputs is not None:
                write_evb_profile_outputs(analysis, replica_out)
        if write_replica_outputs is not None:
            write_replica_outputs(aggregate, base)

    return result


__all__ = ["analyze_evb_profile_replicas_result"]
