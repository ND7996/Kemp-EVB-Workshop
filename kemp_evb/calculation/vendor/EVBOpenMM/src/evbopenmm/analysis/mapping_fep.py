from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import MappingProxyType

from openmm import unit

from ..constants import REPORTING_ENERGY_UNIT
from ..models import FEPSampleRecord, FEPWindowRecord
from ..types import PathLike
from .binning import resolve_gap_bounds
from .energetics import effective_energy_from_state_energies, energy_scale_to_kj, kj_to_kcal, log_mean_exp


def analyze_mapping_fep_umbrella_records(
    window_records: Sequence[FEPWindowRecord],
    *,
    state_names: Sequence[str],
    kT_kj_mol: float,
    skip_samples: int,
    n_gap_bins: int,
    min_points_per_bin: int,
    input_energy_unit: str,
    gap_range: tuple[float, float] | None = None,
    gap_percentile_range: tuple[float, float] | None = None,
) -> dict[str, object]:
    if skip_samples < 0:
        raise ValueError("'skip_samples' must be >= 0.")
    if n_gap_bins <= 0:
        raise ValueError("'n_gap_bins' must be > 0.")
    if min_points_per_bin <= 0:
        raise ValueError("'min_points_per_bin' must be > 0.")
    if len(window_records) < 2:
        raise ValueError("At least two windows are required for FEP analysis.")

    n_states = len(state_names)
    if n_states < 2:
        raise ValueError("At least two states are required for FEP analysis.")
    sample_scale_to_kj = energy_scale_to_kj(input_energy_unit)

    kT = float(kT_kj_mol)
    if kT <= 0.0:
        raise ValueError("kT must be > 0.")
    beta = 1.0 / kT

    selected_windows: list[tuple[FEPWindowRecord, tuple[FEPSampleRecord, ...]]] = []
    for window in window_records:
        selected = tuple(window.samples[skip_samples:])
        if not selected:
            raise ValueError(f"Window {window.index} has no samples after skipping {skip_samples} entries.")
        selected_windows.append((window, selected))

    n_windows = len(selected_windows)
    forward_edges: list[float] = []
    reverse_edges: list[float] = []
    for edge in range(n_windows - 1):
        window_i, samples_i = selected_windows[edge]
        window_j, samples_j = selected_windows[edge + 1]

        lambda_i = window_i.lambda_values
        lambda_j = window_j.lambda_values
        if len(lambda_i) != n_states or len(lambda_j) != n_states:
            raise ValueError("Window lambda vectors must match the number of EVB states.")

        forward_dv = [
            effective_energy_from_state_energies(
                lambda_j,
                tuple(sample_scale_to_kj * value for value in sample.state_energies_kj_mol),
            )
            - (sample_scale_to_kj * sample.effective_energy_kj_mol)
            for sample in samples_i
        ]
        reverse_dv = [
            effective_energy_from_state_energies(
                lambda_i,
                tuple(sample_scale_to_kj * value for value in sample.state_energies_kj_mol),
            )
            - (sample_scale_to_kj * sample.effective_energy_kj_mol)
            for sample in samples_j
        ]

        dG_forward = -kT * log_mean_exp(tuple(-beta * delta for delta in forward_dv))
        dG_reverse = -kT * log_mean_exp(tuple(-beta * delta for delta in reverse_dv))
        forward_edges.append(dG_forward)
        reverse_edges.append(dG_reverse)

    cumulative_forward = [0.0] * n_windows
    cumulative_reverse = [0.0] * n_windows
    cumulative_average = [0.0] * n_windows
    for edge, value in enumerate(forward_edges):
        cumulative_forward[edge + 1] = cumulative_forward[edge] + value
    for edge in range(n_windows - 2, -1, -1):
        cumulative_reverse[edge] = cumulative_reverse[edge + 1] + reverse_edges[edge]
    for edge in range(n_windows - 1):
        cumulative_average[edge + 1] = cumulative_average[edge] + 0.5 * (
            forward_edges[edge] - reverse_edges[edge]
        )

    all_samples = [sample for _, samples in selected_windows for sample in samples]
    gap_values = [sample_scale_to_kj * sample.reaction_coordinate_kj_mol for sample in all_samples]
    # Bounds restrict the HISTOGRAM only -- the forward/reverse exponential averages above are
    # already accumulated over every sample, so the free-energy estimator is untouched.
    gap_min, gap_max, n_gap_excluded, gap_rule = resolve_gap_bounds(
        gap_values,
        gap_range=gap_range,
        gap_percentile_range=gap_percentile_range,
        scale_to_internal=sample_scale_to_kj,
    )
    gap_width = (gap_max - gap_min) / float(n_gap_bins)

    internal_bin_count = n_gap_bins + 1

    def _gap_bin_index(value: float) -> int:
        raw = int((value - gap_min) / gap_width)
        if raw < 0:
            return 0
        if raw > n_gap_bins:
            return n_gap_bins
        return raw

    weighted_state_free_energies = [[0.0 for _ in range(internal_bin_count)] for _ in range(n_states)]
    total_bin_counts = [0 for _ in range(internal_bin_count)]

    for window_idx, (_, samples) in enumerate(selected_windows):
        per_bin_counts = [0 for _ in range(internal_bin_count)]
        mean_dvv = [[0.0 for _ in range(internal_bin_count)] for _ in range(n_states)]
        for sample in samples:
            reaction_coordinate = sample_scale_to_kj * sample.reaction_coordinate_kj_mol
            if reaction_coordinate < gap_min or reaction_coordinate > gap_max:
                continue        # outside the requested histogram bounds; excluded, not clamped
            effective_energy = sample_scale_to_kj * sample.effective_energy_kj_mol
            state_energies = tuple(sample_scale_to_kj * value for value in sample.state_energies_kj_mol)
            bin_idx = _gap_bin_index(reaction_coordinate)
            per_bin_counts[bin_idx] += 1
            for state_idx in range(n_states):
                mean_dvv[state_idx][bin_idx] += state_energies[state_idx] - effective_energy

        for bin_idx, count in enumerate(per_bin_counts):
            if count == 0:
                continue
            inv_count = 1.0 / float(count)
            for state_idx in range(n_states):
                mean_dvv[state_idx][bin_idx] *= inv_count

        exp_arguments: list[list[list[float]]] = [[[] for _ in range(internal_bin_count)] for _ in range(n_states)]
        for sample in samples:
            reaction_coordinate = sample_scale_to_kj * sample.reaction_coordinate_kj_mol
            if reaction_coordinate < gap_min or reaction_coordinate > gap_max:
                continue        # outside the requested histogram bounds; excluded, not clamped
            effective_energy = sample_scale_to_kj * sample.effective_energy_kj_mol
            state_energies = tuple(sample_scale_to_kj * value for value in sample.state_energies_kj_mol)
            bin_idx = _gap_bin_index(reaction_coordinate)
            for state_idx in range(n_states):
                dvv_centered = state_energies[state_idx] - effective_energy - mean_dvv[state_idx][bin_idx]
                exp_arguments[state_idx][bin_idx].append(-dvv_centered / kT)

        for bin_idx in range(n_gap_bins):
            count = per_bin_counts[bin_idx]
            if count < min_points_per_bin:
                continue
            total_bin_counts[bin_idx] += count
            for state_idx in range(n_states):
                values = exp_arguments[state_idx][bin_idx]
                if not values:
                    continue
                log_avg_exp = log_mean_exp(tuple(values))
                free_energy_offset = -kT * log_avg_exp + mean_dvv[state_idx][bin_idx]
                state_free_energy = cumulative_average[window_idx] + free_energy_offset
                weighted_state_free_energies[state_idx][bin_idx] += state_free_energy * float(count)

    state_profiles: dict[str, list[float | None]] = {state: [None] * n_gap_bins for state in state_names}
    state_profiles_normalized: dict[str, list[float | None]] = {state: [None] * n_gap_bins for state in state_names}
    for state_idx, state_name in enumerate(state_names):
        for bin_idx in range(n_gap_bins):
            count = total_bin_counts[bin_idx]
            if count < min_points_per_bin:
                continue
            state_profiles[state_name][bin_idx] = weighted_state_free_energies[state_idx][bin_idx] / float(count)

        finite_values = [value for value in state_profiles[state_name] if value is not None]
        if not finite_values:
            continue
        profile_min = min(finite_values)
        for bin_idx, value in enumerate(state_profiles[state_name]):
            if value is None:
                continue
            state_profiles_normalized[state_name][bin_idx] = value - profile_min

    gap_edges = [gap_min + gap_width * float(idx) for idx in range(n_gap_bins + 1)]
    gap_centers = [gap_min + gap_width * (float(idx) + 0.5) for idx in range(n_gap_bins)]

    window_summaries: list[dict[str, object]] = []
    for idx, (window, samples) in enumerate(selected_windows):
        lambda_map = {state: window.lambda_values[i] for i, state in enumerate(state_names)}
        mean_effective = sum(sample_scale_to_kj * sample.effective_energy_kj_mol for sample in samples) / float(
            len(samples)
        )
        mean_coordinate = sum(sample_scale_to_kj * sample.reaction_coordinate_kj_mol for sample in samples) / float(
            len(samples)
        )
        window_summaries.append(
            {
                "index": idx,
                "lambda_values": lambda_map,
                "n_samples": len(samples),
                "mean_effective_energy_kcal_mol": kj_to_kcal(mean_effective),
                "mean_reaction_coordinate_kcal_mol": kj_to_kcal(mean_coordinate),
                "cumulative_forward_kcal_mol": kj_to_kcal(cumulative_forward[idx]),
                "cumulative_reverse_kcal_mol": kj_to_kcal(cumulative_reverse[idx]),
                "cumulative_average_kcal_mol": kj_to_kcal(cumulative_average[idx]),
            }
        )

    forward_edges_kcal = tuple(kj_to_kcal(value) for value in forward_edges)
    reverse_edges_kcal = tuple(kj_to_kcal(value) for value in reverse_edges)
    cumulative_forward_kcal = tuple(kj_to_kcal(value) for value in cumulative_forward)
    cumulative_reverse_kcal = tuple(kj_to_kcal(value) for value in cumulative_reverse)
    cumulative_average_kcal = tuple(kj_to_kcal(value) for value in cumulative_average)
    gap_edges_kcal = tuple(kj_to_kcal(value) for value in gap_edges)
    gap_centers_kcal = tuple(kj_to_kcal(value) for value in gap_centers)
    state_profiles_kcal = {
        state: tuple(None if value is None else kj_to_kcal(value) for value in values)
        for state, values in state_profiles.items()
    }
    state_profiles_normalized_kcal = {
        state: tuple(None if value is None else kj_to_kcal(value) for value in values)
        for state, values in state_profiles_normalized.items()
    }

    return {
        "state_names": tuple(state_names),
        "energy_unit": REPORTING_ENERGY_UNIT,
        "kT_kcal_mol": kj_to_kcal(kT),
        "skip_samples": int(skip_samples),
        "windows": tuple(window_records),
        "window_summaries": tuple(window_summaries),
        "fep": {
            "forward_edge_dg_kcal_mol": forward_edges_kcal,
            "reverse_edge_dg_kcal_mol": reverse_edges_kcal,
            "cumulative_forward_kcal_mol": cumulative_forward_kcal,
            "cumulative_reverse_kcal_mol": cumulative_reverse_kcal,
            "cumulative_average_kcal_mol": cumulative_average_kcal,
        },
        "umbrella": {
            "gap_edges_kcal_mol": gap_edges_kcal,
            "gap_centers_kcal_mol": gap_centers_kcal,
            "bin_counts": tuple(total_bin_counts[:n_gap_bins]),
            "state_free_energies_kcal_mol": MappingProxyType(state_profiles_kcal),
            "state_free_energies_normalized_kcal_mol": MappingProxyType(state_profiles_normalized_kcal),
        },
    }


def write_mapping_fep_outputs(analysis: Mapping[str, object], *, output_directory: PathLike) -> None:
    outdir = Path(output_directory)
    outdir.mkdir(parents=True, exist_ok=True)

    state_names = tuple(analysis["state_names"])
    window_summaries = tuple(analysis["window_summaries"])
    fep_data = analysis["fep"]
    umbrella_data = analysis["umbrella"]

    def _fmt(value: float, *, digits: int = 8) -> str:
        return f"{float(value):.{digits}f}"

    def _fmt_optional(value: float | None, *, digits: int = 8) -> str:
        if value is None:
            return ""
        return _fmt(value, digits=digits)

    window_path = outdir / "fep_windows.csv"
    with window_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header = ["window_index", *[f"lambda_{name}" for name in state_names], "n_samples"]
        header.extend(
            [
                "mean_effective_energy",
                "mean_reaction_coordinate",
                "cumulative_forward",
                "cumulative_reverse",
                "cumulative_average",
            ]
        )
        writer.writerow(header)
        for summary in window_summaries:
            row = [summary["index"]]
            lambda_values = summary["lambda_values"]
            for state in state_names:
                row.append(_fmt(lambda_values[state], digits=6))
            row.extend(
                [
                    summary["n_samples"],
                    _fmt(summary["mean_effective_energy_kcal_mol"]),
                    _fmt(summary["mean_reaction_coordinate_kcal_mol"]),
                    _fmt(summary["cumulative_forward_kcal_mol"]),
                    _fmt(summary["cumulative_reverse_kcal_mol"]),
                    _fmt(summary["cumulative_average_kcal_mol"]),
                ]
            )
            writer.writerow(row)

    edges_path = outdir / "fep_edges.csv"
    with edges_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["edge_index", "from_window", "to_window", "forward_dg", "reverse_dg", "average_increment"])
        forward = tuple(fep_data["forward_edge_dg_kcal_mol"])
        reverse = tuple(fep_data["reverse_edge_dg_kcal_mol"])
        for idx, (dgf, dgr) in enumerate(zip(forward, reverse)):
            writer.writerow([idx, idx, idx + 1, _fmt(dgf), _fmt(dgr), _fmt(0.5 * (dgf - dgr))])

    umbrella_path = outdir / "umbrella_profile.csv"
    with umbrella_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header = ["bin_index", "gap_center", "n_points"]
        header.extend([f"{state}_free_energy" for state in state_names])
        header.extend([f"{state}_free_energy_norm" for state in state_names])
        writer.writerow(header)

        centers = tuple(umbrella_data["gap_centers_kcal_mol"])
        counts = tuple(umbrella_data["bin_counts"])
        profiles = umbrella_data["state_free_energies_kcal_mol"]
        profiles_norm = umbrella_data["state_free_energies_normalized_kcal_mol"]

        for idx, center in enumerate(centers):
            row: list[object] = [idx, _fmt(center), counts[idx]]
            for state in state_names:
                row.append(_fmt_optional(profiles[state][idx]))
            for state in state_names:
                row.append(_fmt_optional(profiles_norm[state][idx]))
            writer.writerow(row)


__all__ = ["analyze_mapping_fep_umbrella_records", "write_mapping_fep_outputs"]
