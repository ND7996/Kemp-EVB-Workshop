from __future__ import annotations

import warnings

import csv
import json
from bisect import bisect_left
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import MappingProxyType

from ..constants import KJ_PER_KCAL, REPORTING_ENERGY_UNIT
from ..models import FEPSampleRecord, FEPWindowRecord
from ..types import PathLike
from .binning import resolve_gap_bounds
from .energetics import (
    energy_scale_to_kj,
    evb_ground_state_energy_distance_coupling_kj_mol,
    evb_ground_state_energy_gaussian_coupling_nstate_kj_mol,
    evb_ground_state_energy_nstate_kj_mol,
    kj_to_kcal,
    log_mean_exp,
    mean_std_sem,
)


def interpolate_profile_to_grid(
    source_x: Sequence[float],
    source_y: Sequence[float | None],
    target_x: Sequence[float],
    *,
    tolerance: float,
) -> tuple[float | None, ...]:
    if len(source_x) != len(source_y):
        raise ValueError("Source profile grid and values must have identical lengths.")
    if not source_x:
        return tuple(None for _ in target_x)

    finite_pairs = [
        (float(x_value), float(y_value))
        for x_value, y_value in zip(source_x, source_y)
        if y_value is not None
    ]
    if not finite_pairs:
        return tuple(None for _ in target_x)

    finite_pairs.sort(key=lambda item: item[0])
    xs = [item[0] for item in finite_pairs]
    ys = [item[1] for item in finite_pairs]

    if len(xs) == 1:
        single_x = xs[0]
        single_y = ys[0]
        out_values: list[float | None] = []
        for x_target in target_x:
            x_value = float(x_target)
            if abs(x_value - single_x) <= tolerance:
                out_values.append(single_y)
            else:
                out_values.append(None)
        return tuple(out_values)

    out_values: list[float | None] = []
    x_min = xs[0]
    x_max = xs[-1]

    for x_target in target_x:
        x_value = float(x_target)
        if x_value < (x_min - tolerance) or x_value > (x_max + tolerance):
            out_values.append(None)
            continue
        if abs(x_value - x_min) <= tolerance:
            out_values.append(ys[0])
            continue
        if abs(x_value - x_max) <= tolerance:
            out_values.append(ys[-1])
            continue

        right = bisect_left(xs, x_value)
        if right <= 0:
            out_values.append(ys[0])
            continue
        if right >= len(xs):
            out_values.append(ys[-1])
            continue

        left = right - 1
        x_left = xs[left]
        x_right = xs[right]
        y_left = ys[left]
        y_right = ys[right]

        if abs(x_value - x_left) <= tolerance:
            out_values.append(y_left)
            continue
        if abs(x_value - x_right) <= tolerance:
            out_values.append(y_right)
            continue

        dx = x_right - x_left
        if abs(dx) <= 1e-14:
            out_values.append(0.5 * (y_left + y_right))
            continue

        fraction = (x_value - x_left) / dx
        out_values.append(y_left + fraction * (y_right - y_left))

    return tuple(out_values)


def detect_profile_thermochemistry(
    gap_values: Sequence[float],
    profile_values: Sequence[float],
) -> dict[str, object]:
    if len(gap_values) != len(profile_values):
        raise ValueError("Gap and profile arrays must have identical lengths.")
    if len(gap_values) < 2:
        return {
            "reactant": None,
            "product": None,
            "transition_states": tuple(),
            "warnings": ("Not enough profile points for extrema detection.",),
        }

    derivatives = [profile_values[idx + 1] - profile_values[idx] for idx in range(len(profile_values) - 1)]
    warnings_list: list[str] = []

    rs_index = 0
    rs_warning = True
    for idx in range(len(derivatives) - 1):
        if derivatives[idx] < 0.0 and derivatives[idx + 1] >= 0.0:
            rs_index = idx + 1
            rs_warning = False
            break

    ps_index = len(profile_values) - 1
    ps_warning = True
    for idx in range(len(derivatives) - 1, 0, -1):
        if derivatives[idx] > 0.0 and derivatives[idx - 1] <= 0.0:
            ps_index = idx
            ps_warning = False
            break

    if rs_warning:
        warnings_list.append("Reactant minimum could not be located robustly; using first profile point.")
    if ps_warning:
        warnings_list.append("Product minimum could not be located robustly; using last profile point.")

    ts_indices: list[int] = []
    if len(derivatives) >= 2 and ps_index > rs_index:
        stop = min(ps_index, len(derivatives) - 1)
        for idx in range(rs_index, stop):
            if (derivatives[idx] > 0.0 and derivatives[idx + 1] <= 0.0) or (
                derivatives[idx] >= 0.0 and derivatives[idx + 1] < 0.0
            ):
                ts_indices.append(idx + 1)

    while len(ts_indices) > 1:
        removed = False
        for idx in range(len(ts_indices) - 1):
            delta = profile_values[ts_indices[idx + 1]] - profile_values[ts_indices[idx]]
            if 0.0 <= delta < 0.5:
                ts_indices.pop(idx)
                removed = True
                break
            if -0.5 < delta < 0.0:
                ts_indices.pop(idx + 1)
                removed = True
                break
        if not removed:
            break

    while len(ts_indices) > 1:
        first_barrier = profile_values[ts_indices[0]] - profile_values[rs_index]
        last_barrier = profile_values[ts_indices[-1]] - profile_values[ps_index]
        if first_barrier < last_barrier and first_barrier < 0.5:
            ts_indices.pop(0)
        elif first_barrier > last_barrier and last_barrier < 0.5:
            ts_indices.pop(-1)
        else:
            break

    if gap_values[rs_index] > 0.0:
        warnings_list.append("Reactant basin appears on the product side of the gap coordinate.")
    if gap_values[ps_index] < 0.0:
        warnings_list.append("Product basin appears on the reactant side of the gap coordinate.")

    reactant = {
        "index": int(rs_index),
        "gap_kcal_mol": float(gap_values[rs_index]),
        "free_energy_kcal_mol": float(profile_values[rs_index]),
        "warning": bool(rs_warning),
    }
    product = {
        "index": int(ps_index),
        "gap_kcal_mol": float(gap_values[ps_index]),
        "free_energy_kcal_mol": float(profile_values[ps_index]),
        "warning": bool(ps_warning),
    }
    transition_states = tuple(
        {
            "index": int(idx),
            "gap_kcal_mol": float(gap_values[idx]),
            "free_energy_kcal_mol": float(profile_values[idx]),
        }
        for idx in ts_indices
    )
    return {
        "reactant": reactant,
        "product": product,
        "transition_states": transition_states,
        "warnings": tuple(warnings_list),
    }


def build_diabatic_profile(
    window_records: Sequence[FEPWindowRecord],
    *,
    state_names: Sequence[str],
    base_analysis: Mapping[str, object],
    reaction_coordinate_states: tuple[int, int] = (0, 1),
    skip_samples: int = 0,
    n_gap_bins: int = 40,
    min_points_per_bin: int = 20,
    input_energy_unit: str = "kcal/mol",
    gap_range: tuple[float, float] | None = None,
    gap_percentile_range: tuple[float, float] | None = None,
) -> dict[str, object]:
    """Build the alpha/H12-INDEPENDENT sampling product of an EVB mapping-FEP run.

    Returns the RAW per-state diabatic free-energy curves ``F_k(gap)`` binned by
    the RAW energy gap ``eps_i - eps_j`` (``(i, j) = reaction_coordinate_states``),
    plus bin edges/centers, bin counts and the per-bin mean reactive-pair distance
    (for the distance-coupling analyser). It applies **no** diabatic offset (alpha),
    **no** coupling (H12) and computes **no** adiabatic ground curve -- those are the
    job of the EVB analyser :func:`analyze_evb_ground`. Because nothing alpha/H12
    dependent is folded in, the stored curves are re-analysable at any (alpha, H12)
    without the double-count trap.
    """
    n_states = len(state_names)
    if n_states < 2:
        raise ValueError("EVB profile analysis requires at least two diabatic states.")
    ref_i, ref_j = reaction_coordinate_states
    if ref_i < 0 or ref_i >= n_states or ref_j < 0 or ref_j >= n_states or ref_i == ref_j:
        raise ValueError(
            f"'reaction_coordinate_states' must be two distinct indices in [0, {n_states}); "
            f"received {reaction_coordinate_states}."
        )

    sample_scale_to_kj = energy_scale_to_kj(input_energy_unit)

    selected_windows: list[tuple[FEPWindowRecord, tuple[FEPSampleRecord, ...]]] = []
    for window in window_records:
        selected = tuple(window.samples[skip_samples:])
        if not selected:
            continue
        selected_windows.append((window, selected))
    if not selected_windows:
        raise ValueError("No samples are available for EVB profile analysis.")

    cumulative_average_kj = tuple(
        float(value) * KJ_PER_KCAL for value in base_analysis["fep"]["cumulative_average_kcal_mol"]
    )
    if len(cumulative_average_kj) != len(selected_windows):
        raise RuntimeError("Internal mismatch between FEP summary and selected EVB windows.")

    def _raw_diabatic_kj(sample: FEPSampleRecord) -> tuple[float, ...]:
        raw = sample.state_energies_kj_mol
        if len(raw) != n_states:
            raise ValueError(
                f"EVB profile analysis requires exactly {n_states} diabatic-state energies "
                f"per sample; received {len(raw)}."
            )
        return tuple(sample_scale_to_kj * float(raw[k]) for k in range(n_states))

    def _reaction_coordinate_kj(raw: Sequence[float]) -> float:
        return raw[ref_i] - raw[ref_j]

    # RAW gap range (no offsets). Bounds restrict the histogram only.
    _raw_gaps = [
        _reaction_coordinate_kj(_raw_diabatic_kj(sample))
        for _, samples in selected_windows
        for sample in samples
    ]
    gap_min, gap_max, n_gap_excluded, gap_rule = resolve_gap_bounds(
        _raw_gaps,
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

    # Per-state RAW free-energy curves (n_states; NO ground curve here).
    weighted_free_energies = [[0.0 for _ in range(internal_bin_count)] for _ in range(n_states)]
    total_bin_counts = [0 for _ in range(internal_bin_count)]
    distance_sum = [0.0 for _ in range(internal_bin_count)]
    distance_count = [0 for _ in range(internal_bin_count)]
    kT_kj_mol = float(base_analysis["kT_kcal_mol"]) * KJ_PER_KCAL

    for window_idx, (_, samples) in enumerate(selected_windows):
        per_bin_counts = [0 for _ in range(internal_bin_count)]
        mean_dvv = [[0.0 for _ in range(internal_bin_count)] for _ in range(n_states)]

        for sample in samples:
            raw = _raw_diabatic_kj(sample)
            coordinate = _reaction_coordinate_kj(raw)
            if coordinate < gap_min or coordinate > gap_max:
                continue        # outside the requested histogram bounds; excluded, not clamped
            effective_energy = sample_scale_to_kj * sample.effective_energy_kj_mol
            bin_idx = _gap_bin_index(coordinate)
            per_bin_counts[bin_idx] += 1
            for k in range(n_states):
                mean_dvv[k][bin_idx] += raw[k] - effective_energy
            distance = sample.coupling_distance_nm
            if distance is not None:
                distance_sum[bin_idx] += float(distance)
                distance_count[bin_idx] += 1

        for bin_idx, count in enumerate(per_bin_counts):
            if count == 0:
                continue
            inv_count = 1.0 / float(count)
            for k in range(n_states):
                mean_dvv[k][bin_idx] *= inv_count

        exp_arguments: list[list[list[float]]] = [[[] for _ in range(internal_bin_count)] for _ in range(n_states)]
        for sample in samples:
            raw = _raw_diabatic_kj(sample)
            effective_energy = sample_scale_to_kj * sample.effective_energy_kj_mol
            bin_idx = _gap_bin_index(_reaction_coordinate_kj(raw))
            for k in range(n_states):
                centered = (raw[k] - effective_energy) - mean_dvv[k][bin_idx]
                exp_arguments[k][bin_idx].append(-centered / kT_kj_mol)

        for bin_idx in range(n_gap_bins):
            count = per_bin_counts[bin_idx]
            if count < min_points_per_bin:
                continue
            total_bin_counts[bin_idx] += count
            for k in range(n_states):
                values = exp_arguments[k][bin_idx]
                if not values:
                    continue
                log_avg_exp = log_mean_exp(tuple(values))
                free_energy_offset = -kT_kj_mol * log_avg_exp + mean_dvv[k][bin_idx]
                free_energy = cumulative_average_kj[window_idx] + free_energy_offset
                weighted_free_energies[k][bin_idx] += free_energy * float(count)

    curves_kj: list[list[float | None]] = [[None] * n_gap_bins for _ in range(n_states)]
    for bin_idx in range(n_gap_bins):
        count = total_bin_counts[bin_idx]
        if count < min_points_per_bin:
            continue
        for k in range(n_states):
            curves_kj[k][bin_idx] = weighted_free_energies[k][bin_idx] / float(count)

    curves_kcal: list[tuple[float | None, ...]] = [
        tuple(None if value is None else kj_to_kcal(value) for value in curve) for curve in curves_kj
    ]
    mean_distance_nm = tuple(
        (distance_sum[idx] / distance_count[idx]) if distance_count[idx] > 0 else None
        for idx in range(n_gap_bins)
    )
    gap_bounds_rule = gap_rule
    n_samples_excluded_by_bounds = int(n_gap_excluded)
    gap_edges_kcal = tuple(kj_to_kcal(gap_min + gap_width * float(idx)) for idx in range(n_gap_bins + 1))
    gap_centers_kcal = tuple(kj_to_kcal(gap_min + gap_width * (float(idx) + 0.5)) for idx in range(n_gap_bins))

    state_curves_kcal = {name: curves_kcal[k] for k, name in enumerate(state_names)}

    return {
        "state_names": tuple(state_names),
        "reaction_coordinate_states": (int(ref_i), int(ref_j)),
        "n_gap_bins": int(n_gap_bins),
        "min_points_per_bin": int(min_points_per_bin),
        "skip_samples": int(skip_samples),
        "kT_kcal_mol": float(base_analysis["kT_kcal_mol"]),
        "gap_min_kj_mol": float(gap_min),
        "gap_max_kj_mol": float(gap_max),
        "gap_width_kj_mol": float(gap_width),
        "gap_bounds_rule": gap_bounds_rule,
        "n_samples_excluded_by_bounds": n_samples_excluded_by_bounds,
        "gap_edges_kcal_mol": gap_edges_kcal,
        "gap_centers_kcal_mol": gap_centers_kcal,
        "bin_counts": tuple(total_bin_counts[:n_gap_bins]),
        "state_free_energies_kcal_mol": state_curves_kcal,
        "state_free_energies_kj_mol": [tuple(curve) for curve in curves_kj],
        "mean_coupling_distance_nm": mean_distance_nm,
        "windows": tuple(window_records),
        "window_summaries": tuple(base_analysis["window_summaries"]),
        "fep": base_analysis["fep"],
        "umbrella": base_analysis["umbrella"],
        "storage": "raw",
    }


def resolve_curve_reference(reference: str) -> str:
    """Normalise a curve-reference choice to one of {'raw','reactant','product'}.

    ``'substrate'`` is accepted as an alias for ``'reactant'``.
    """
    normalized = str(reference).lower()
    if normalized == "substrate":
        normalized = "reactant"
    if normalized not in {"raw", "reactant", "product"}:
        raise ValueError(
            "'reference' must be 'raw', 'reactant'/'substrate', or 'product'; "
            f"received {reference!r}."
        )
    return normalized


def _curve_reference_shift_kcal(
    reference: str,
    reactant: Mapping[str, object] | None,
    product: Mapping[str, object] | None,
    ground_curve: Sequence[float | None],
) -> float:
    """Zero-reference shift (kcal/mol) to subtract from every reported curve.

    ``'raw'`` -> 0; ``'reactant'``/``'product'`` -> that basin's detected minimum
    free energy (falling back to the lowest ground point if the basin was not
    identified). dG-dagger / dG0 are differences and so are invariant to this.
    """
    if reference == "raw":
        return 0.0
    detail = reactant if reference == "reactant" else product
    if detail is not None:
        return float(detail["free_energy_kcal_mol"])
    finite = [value for value in ground_curve if value is not None]
    return min(finite) if finite else 0.0


def analyze_evb_ground(
    raw_profile: Mapping[str, object],
    *,
    diabatic_offsets_kcal_mol: Sequence[float] | None = None,
    coupling_matrix_kcal_mol: Sequence[Sequence[float]] | None = None,
    coupling_form: str = "constant",
    coupling_width_matrix_per_kcal2_mol2: Sequence[Sequence[float]] | None = None,
    coupling_decay_per_nm: float = 0.0,
    coupling_reference_distance_nm: float = 0.0,
    reference: str = "raw",
) -> dict[str, object]:
    """EVB analyser (N2): apply diabatic offsets (alpha) + coupling (H12) to the
    RAW per-state free-energy curves from :func:`build_diabatic_profile`.

    The adiabatic ground curve is formed POINTWISE per gap bin as the lowest
    eigenvalue of ``diag(F_k(gap) + offset_k) + H_offdiag`` (the standard
    Warshel free-energy-curve construction), then thermochemistry (dG-dagger,
    dG0) is detected on it. This is the ONLY place alpha/H12 enter, so the same
    raw profile can be re-analysed at any (alpha, H12) without double-counting.
    """
    state_names = tuple(raw_profile["state_names"])
    n_states = len(state_names)
    ref_i, ref_j = tuple(raw_profile["reaction_coordinate_states"])  # type: ignore[misc]
    n_gap_bins = int(raw_profile["n_gap_bins"])
    min_points_per_bin = int(raw_profile["min_points_per_bin"])
    curves_kj = [list(curve) for curve in raw_profile["state_free_energies_kj_mol"]]  # type: ignore[arg-type]
    gap_centers_kcal = tuple(raw_profile["gap_centers_kcal_mol"])  # type: ignore[arg-type]
    gap_min = float(raw_profile["gap_min_kj_mol"])
    gap_max = float(raw_profile["gap_max_kj_mol"])
    mean_distance_nm = tuple(raw_profile.get("mean_coupling_distance_nm") or [None] * n_gap_bins)

    if diabatic_offsets_kcal_mol is None:
        offsets_kcal = [0.0] * n_states
    else:
        offsets_kcal = [float(value) for value in diabatic_offsets_kcal_mol]
        if len(offsets_kcal) != n_states:
            raise ValueError(
                f"'diabatic_offsets_kcal_mol' must have length {n_states} (one per state); "
                f"received {len(offsets_kcal)}."
            )

    if coupling_matrix_kcal_mol is None:
        coupling_kcal = [[0.0] * n_states for _ in range(n_states)]
    else:
        coupling_kcal = [list(map(float, row)) for row in coupling_matrix_kcal_mol]
        if len(coupling_kcal) != n_states or any(len(row) != n_states for row in coupling_kcal):
            raise ValueError(f"'coupling_matrix_kcal_mol' must be {n_states}x{n_states}.")

    coupling_form_normalized = str(coupling_form).lower()
    if coupling_form_normalized not in {"constant", "gaussian", "distance"}:
        raise ValueError(
            f"'coupling_form' must be 'constant', 'gaussian' or 'distance'; received {coupling_form!r}."
        )
    if coupling_form_normalized == "distance" and n_states != 2:
        raise ValueError(
            "coupling_form='distance' (the Warshel/Q6 H12(r)=A*exp(-mu*(r-r0)) form) is "
            f"defined for two diabatic states only; received {n_states}."
        )

    if coupling_width_matrix_per_kcal2_mol2 is None:
        width_kcal: list[list[float]] | None = None
    else:
        width_kcal = [list(map(float, row)) for row in coupling_width_matrix_per_kcal2_mol2]
        if len(width_kcal) != n_states or any(len(row) != n_states for row in width_kcal):
            raise ValueError(
                f"'coupling_width_matrix_per_kcal2_mol2' must be {n_states}x{n_states}."
            )
    if coupling_form_normalized == "gaussian" and width_kcal is None:
        raise ValueError(
            "'coupling_width_matrix_per_kcal2_mol2' is required when coupling_form='gaussian'."
        )

    offsets_kj = tuple(value * KJ_PER_KCAL for value in offsets_kcal)
    coupling_kj = tuple(tuple(value * KJ_PER_KCAL for value in row) for row in coupling_kcal)
    if width_kcal is None:
        width_kj: tuple[tuple[float, ...], ...] | None = None
    else:
        inv_kj_sq = 1.0 / (KJ_PER_KCAL * KJ_PER_KCAL)
        width_kj = tuple(tuple(value * inv_kj_sq for value in row) for row in width_kcal)

    if coupling_form_normalized == "gaussian":
        def _ground_kj(diabatic_kj: Sequence[float], bin_idx: int) -> float | None:
            return evb_ground_state_energy_gaussian_coupling_nstate_kj_mol(diabatic_kj, coupling_kj, width_kj)
    elif coupling_form_normalized == "distance":
        amplitude_kj = coupling_kj[0][1]
        decay_per_nm = float(coupling_decay_per_nm)
        reference_distance_nm = float(coupling_reference_distance_nm)
        if all(distance is None for distance in mean_distance_nm):
            raise ValueError(
                "coupling_form='distance' requires a per-sample reactive-pair distance, but every "
                "bin has coupling_distance_nm=None. Re-run the mapping FEP with "
                "collect_mapping_fep_samples(coupling_distance_atom_pair=(i, j)) so the bond being "
                "made/broken is recorded each frame."
            )

        def _ground_kj(diabatic_kj: Sequence[float], bin_idx: int) -> float | None:
            distance = mean_distance_nm[bin_idx]
            if distance is None:
                return None
            return evb_ground_state_energy_distance_coupling_kj_mol(
                diabatic_kj[0], diabatic_kj[1],
                alpha_kj_mol=0.0,
                coupling_distance_nm=float(distance),
                amplitude_kj_mol=amplitude_kj,
                decay_per_nm=decay_per_nm,
                reference_distance_nm=reference_distance_nm,
            )
    else:
        def _ground_kj(diabatic_kj: Sequence[float], bin_idx: int) -> float | None:
            return evb_ground_state_energy_nstate_kj_mol(diabatic_kj, coupling_kj)

    ground_curve: list[float | None] = [None] * n_gap_bins
    for bin_idx in range(n_gap_bins):
        if any(curves_kj[k][bin_idx] is None for k in range(n_states)):
            continue
        diabatic_kj = [float(curves_kj[k][bin_idx]) + offsets_kj[k] for k in range(n_states)]
        ground_kj = _ground_kj(diabatic_kj, bin_idx)
        ground_curve[bin_idx] = None if ground_kj is None else kj_to_kcal(ground_kj)

    valid_bins = [idx for idx, value in enumerate(ground_curve) if value is not None]
    thermo_warnings: list[str] = []
    reactant: dict[str, object] | None = None
    product: dict[str, object] | None = None
    transition_states: tuple[dict[str, object], ...] = tuple()
    dg_dagger: float | None = None
    dg_zero: float | None = None
    profile_shape: str = "undetermined"
    rate_limiting_ts: dict[str, object] | None = None

    if len(valid_bins) >= 3:
        local_gaps = [gap_centers_kcal[idx] for idx in valid_bins]
        local_profile = [float(ground_curve[idx]) for idx in valid_bins]
        extrema = detect_profile_thermochemistry(local_gaps, local_profile)
        thermo_warnings.extend(extrema["warnings"])

        local_reactant = extrema["reactant"]
        if local_reactant is not None:
            global_idx = valid_bins[int(local_reactant["index"])]
            reactant = {
                "index": int(global_idx),
                "gap_kcal_mol": float(local_reactant["gap_kcal_mol"]),
                "free_energy_kcal_mol": float(local_reactant["free_energy_kcal_mol"]),
                "warning": bool(local_reactant["warning"]),
            }
        local_product = extrema["product"]
        if local_product is not None:
            global_idx = valid_bins[int(local_product["index"])]
            product = {
                "index": int(global_idx),
                "gap_kcal_mol": float(local_product["gap_kcal_mol"]),
                "free_energy_kcal_mol": float(local_product["free_energy_kcal_mol"]),
                "warning": bool(local_product["warning"]),
            }
        ts_entries: list[dict[str, object]] = []
        for entry in extrema["transition_states"]:
            global_idx = valid_bins[int(entry["index"])]
            ts_entries.append(
                {
                    "index": int(global_idx),
                    "gap_kcal_mol": float(entry["gap_kcal_mol"]),
                    "free_energy_kcal_mol": float(entry["free_energy_kcal_mol"]),
                }
            )
        transition_states = tuple(ts_entries)
    else:
        thermo_warnings.append("Not enough populated bins for robust RS/TS/PS detection.")

    if reactant is not None and product is not None and len(transition_states) >= 1:
        # Rate-limiting TS = highest free-energy maximum between the reactant and
        # product basins. With a single TS this is that TS (unchanged behaviour);
        # with several surviving maxima we take the HIGHEST (rate-limiting) one
        # rather than returning None, so an ill-conditioned (noisy multi-max)
        # ground profile still yields a well-defined forward barrier.
        rate_limiting_ts = max(
            transition_states, key=lambda ts: float(ts["free_energy_kcal_mol"])
        )
        dg_dagger = float(rate_limiting_ts["free_energy_kcal_mol"] - reactant["free_energy_kcal_mol"])
        dg_zero = float(product["free_energy_kcal_mol"] - reactant["free_energy_kcal_mol"])
        profile_shape = "double_well"
        if len(transition_states) > 1:
            thermo_warnings.append(
                f"Multiple ({len(transition_states)}) transition-state maxima detected; "
                "using the highest (rate-limiting) one for dG-dagger."
            )
    elif reactant is None or product is None:
        profile_shape = "unresolved"
        thermo_warnings.append("Reactant and/or product basin not resolved; dG-dagger undefined.")
    else:
        # A reactant and product basin exist but no interior maximum separates them:
        # the ground-state profile is flat / monotonic. dG-dagger is genuinely
        # undefined (no barrier) -- report it as such rather than as a number.
        profile_shape = "flat_or_monotonic"
        thermo_warnings.append(
            "No transition-state maximum between the reactant and product basins "
            "(flat/monotonic ground-state profile); dG-dagger is undefined."
        )

    # --- gap-coverage / truncation diagnostic (see analyze_evb_profile_records history) ---
    gap_coverage: dict[str, object] | None = None
    if valid_bins:
        resolved_lo, resolved_hi = valid_bins[0], valid_bins[-1]
        n_dropped_low = int(resolved_lo)
        n_dropped_high = int((n_gap_bins - 1) - resolved_hi)
        truncated_basins: list[str] = []
        if product is not None and int(product["index"]) >= resolved_hi and n_dropped_high > 0:
            msg = (
                f"Product basin at the high-gap edge of the resolved range "
                f"(bin {resolved_hi}/{n_gap_bins - 1}); {n_dropped_high} higher-gap bin(s) dropped "
                f"for <{min_points_per_bin} pts. dG0 may be over-estimated (truncated product well) -- "
                f"sample the product side more densely/longer."
            )
            thermo_warnings.append(msg)
            warnings.warn(msg, RuntimeWarning, stacklevel=2)
            truncated_basins.append("product")
        if reactant is not None and int(reactant["index"]) <= resolved_lo and n_dropped_low > 0:
            msg = (
                f"Reactant basin at the low-gap edge of the resolved range (bin {resolved_lo}); "
                f"{n_dropped_low} lower-gap bin(s) dropped. Reactant well may be truncated."
            )
            thermo_warnings.append(msg)
            warnings.warn(msg, RuntimeWarning, stacklevel=2)
            truncated_basins.append("reactant")
        gap_coverage = {
            "sampled_gap_kcal_mol": (float(kj_to_kcal(gap_min)), float(kj_to_kcal(gap_max))),
            "resolved_gap_kcal_mol": (
                float(gap_centers_kcal[resolved_lo]), float(gap_centers_kcal[resolved_hi])
            ),
            "n_bins": int(n_gap_bins),
            "n_bins_resolved": int(len(valid_bins)),
            "n_bins_dropped_low": n_dropped_low,
            "n_bins_dropped_high": n_dropped_high,
            "truncated_basins": tuple(truncated_basins),
            "gap_bounds_rule": raw_profile.get("gap_bounds_rule"),
            "n_samples_excluded_by_bounds": int(
                raw_profile.get("n_samples_excluded_by_bounds", 0)
            ),
        }

    reference_normalized = resolve_curve_reference(reference)
    reference_shift = _curve_reference_shift_kcal(
        reference_normalized, reactant, product, ground_curve
    )
    evb_ground_referenced = tuple(
        None if value is None else float(value) - reference_shift for value in ground_curve
    )
    return {
        "evb_ground_kcal_mol": tuple(ground_curve),
        "evb_ground_referenced_kcal_mol": evb_ground_referenced,
        "reference": reference_normalized,
        "reference_shift_kcal_mol": float(reference_shift),
        "dg_dagger_kcal_mol": dg_dagger,
        "dg0_kcal_mol": dg_zero,
        "profile_shape": profile_shape,
        "has_barrier": dg_dagger is not None,
        "rate_limiting_transition_state": rate_limiting_ts,
        "reactant": reactant,
        "product": product,
        "transition_states": transition_states,
        "gap_coverage": gap_coverage,
        "warnings": tuple(thermo_warnings),
        "evb_parameters": {
            "diabatic_offsets_kcal_mol": tuple(offsets_kcal),
            "coupling_matrix_kcal_mol": tuple(tuple(row) for row in coupling_kcal),
            "coupling_form": coupling_form_normalized,
            "coupling_width_matrix_per_kcal2_mol2": (
                None if width_kcal is None else tuple(tuple(row) for row in width_kcal)
            ),
            "reaction_coordinate_states": (int(ref_i), int(ref_j)),
            "reference": reference_normalized,
        },
    }


def analyze_evb_profile_records(
    window_records: Sequence[FEPWindowRecord],
    *,
    state_names: Sequence[str],
    base_analysis: Mapping[str, object],
    diabatic_offsets_kcal_mol: Sequence[float] | None = None,
    coupling_matrix_kcal_mol: Sequence[Sequence[float]] | None = None,
    reaction_coordinate_states: tuple[int, int] = (0, 1),
    skip_samples: int = 0,
    n_gap_bins: int = 40,
    min_points_per_bin: int = 20,
    input_energy_unit: str = "kcal/mol",
    gap_range: tuple[float, float] | None = None,
    gap_percentile_range: tuple[float, float] | None = None,
    coupling_form: str = "constant",
    coupling_width_matrix_per_kcal2_mol2: Sequence[Sequence[float]] | None = None,
    coupling_decay_per_nm: float = 0.0,
    coupling_reference_distance_nm: float = 0.0,
    reference: str = "raw",
) -> dict[str, object]:
    """Analyze an EVB profile for an arbitrary number (N >= 2) of diabatic states.

    Thin wrapper over the split builder/analyser:
    :func:`build_diabatic_profile` produces the alpha/H12-INDEPENDENT raw per-state
    free-energy curves binned by the RAW gap ``eps_i - eps_j``, and
    :func:`analyze_evb_ground` applies the diabatic offsets (alpha) + coupling (H12)
    to form the adiabatic ground curve (pointwise lowest eigenvalue of
    ``diag(E + offset) + H_offdiag``) and its thermochemistry.

    The stored per-state curves and ``gap_center`` are therefore RAW (no alpha);
    ``evb_parameters["storage"] == "raw"`` records this. Re-analysing the same
    profile at a different (alpha, H12) -- e.g. via ``analyze_evb_ground`` -- is
    exact and never double-counts the offset.

    For a two-state run, pass ``diabatic_offsets_kcal_mol = [0, alpha]``,
    ``coupling_matrix_kcal_mol = [[0, Hij], [Hij, 0]]``, and
    ``reaction_coordinate_states = (0, 1)``. ``coupling_form`` selects the
    off-diagonal model ("constant", "gaussian", "distance"); see
    :func:`analyze_evb_ground`.
    """
    if str(coupling_form).lower() == "distance" and len(state_names) != 2:
        raise ValueError(
            "coupling_form='distance' (the Warshel/Q6 H12(r)=A*exp(-mu*(r-r0)) form) is "
            f"defined for two diabatic states only; received {len(state_names)}."
        )
    raw = build_diabatic_profile(
        window_records,
        state_names=state_names,
        base_analysis=base_analysis,
        reaction_coordinate_states=reaction_coordinate_states,
        skip_samples=skip_samples,
        n_gap_bins=n_gap_bins,
        min_points_per_bin=min_points_per_bin,
        input_energy_unit=input_energy_unit,
        gap_range=gap_range,
        gap_percentile_range=gap_percentile_range,
    )
    return analyze_evb_profile_from_raw(
        raw,
        window_records=window_records,
        base_analysis=base_analysis,
        state_names=state_names,
        skip_samples=skip_samples,
        diabatic_offsets_kcal_mol=diabatic_offsets_kcal_mol,
        coupling_matrix_kcal_mol=coupling_matrix_kcal_mol,
        coupling_form=coupling_form,
        coupling_width_matrix_per_kcal2_mol2=coupling_width_matrix_per_kcal2_mol2,
        coupling_decay_per_nm=coupling_decay_per_nm,
        coupling_reference_distance_nm=coupling_reference_distance_nm,
        reference=reference,
    )


def analyze_evb_profile_from_raw(
    raw_profile: Mapping[str, object],
    *,
    window_records: Sequence[FEPWindowRecord],
    base_analysis: Mapping[str, object],
    state_names: Sequence[str],
    skip_samples: int = 0,
    diabatic_offsets_kcal_mol: Sequence[float] | None = None,
    coupling_matrix_kcal_mol: Sequence[Sequence[float]] | None = None,
    coupling_form: str = "constant",
    coupling_width_matrix_per_kcal2_mol2: Sequence[Sequence[float]] | None = None,
    coupling_decay_per_nm: float = 0.0,
    coupling_reference_distance_nm: float = 0.0,
    reference: str = "raw",
) -> dict[str, object]:
    """Apply the EVB analyser (alpha/H12 + zero-reference) to an ALREADY-BUILT raw
    profile from :func:`build_diabatic_profile`, returning the full analysis dict --
    identical to :func:`analyze_evb_profile_records` for the same inputs.

    Split out so an optimiser (EVB parameter fit) can build the alpha/H12-INDEPENDENT
    raw profile ONCE and re-run only this cheap step per ``(offsets, coupling)`` trial,
    instead of re-binning every iteration. ``window_records`` / ``base_analysis`` feed
    only the metadata blocks (windows, FEP/umbrella summaries, kT).
    """
    raw = raw_profile
    ground = analyze_evb_ground(
        raw,
        diabatic_offsets_kcal_mol=diabatic_offsets_kcal_mol,
        coupling_matrix_kcal_mol=coupling_matrix_kcal_mol,
        coupling_form=coupling_form,
        coupling_width_matrix_per_kcal2_mol2=coupling_width_matrix_per_kcal2_mol2,
        coupling_decay_per_nm=coupling_decay_per_nm,
        coupling_reference_distance_nm=coupling_reference_distance_nm,
        reference=reference,
    )

    state_names_t = tuple(state_names)
    ground_curve = ground["evb_ground_kcal_mol"]
    reactant = ground["reactant"]

    curves_kcal: list[tuple[float | None, ...]] = [
        raw["state_free_energies_kcal_mol"][name] for name in state_names_t
    ]
    curves_kcal.append(tuple(ground_curve))

    # Single, user-selectable zero reference (default 'raw' = no shift). dG-dagger
    # / dG0 are differences and are invariant to it; this only sets the reported
    # curve baseline. See analyze_evb_ground / resolve_curve_reference.
    reference_shift = float(ground["reference_shift_kcal_mol"])

    def _apply_reference(values: Sequence[float | None]) -> tuple[float | None, ...]:
        if reference_shift == 0.0:
            return tuple(values)
        return tuple(None if value is None else float(value) - reference_shift for value in values)

    curves_referenced = [_apply_reference(curve) for curve in curves_kcal]

    curve_names = (*state_names_t, "evb_ground")
    state_profiles_dict: dict[str, tuple[float | None, ...]] = {}
    state_profiles_referenced_dict: dict[str, tuple[float | None, ...]] = {}
    for idx, name in enumerate(curve_names):
        state_profiles_dict[name] = curves_kcal[idx]
        state_profiles_referenced_dict[name] = curves_referenced[idx]

    evb_parameters = dict(ground["evb_parameters"])
    evb_parameters["storage"] = "raw"

    return {
        "state_names": state_names_t,
        "energy_unit": REPORTING_ENERGY_UNIT,
        "kT_kcal_mol": float(base_analysis["kT_kcal_mol"]),
        "skip_samples": int(skip_samples),
        "windows": tuple(window_records),
        "window_summaries": tuple(base_analysis["window_summaries"]),
        "fep": base_analysis["fep"],
        "umbrella": base_analysis["umbrella"],
        "evb_parameters": evb_parameters,
        "profile": {
            "gap_edges_kcal_mol": raw["gap_edges_kcal_mol"],
            "gap_centers_kcal_mol": raw["gap_centers_kcal_mol"],
            "bin_counts": raw["bin_counts"],
            "state_free_energies_kcal_mol": MappingProxyType(state_profiles_dict),
            "state_free_energies_referenced_kcal_mol": MappingProxyType(state_profiles_referenced_dict),
        },
        "dg_dagger_kcal_mol": ground["dg_dagger_kcal_mol"],
        "dg0_kcal_mol": ground["dg0_kcal_mol"],
        "profile_shape": ground["profile_shape"],
        "has_barrier": ground["has_barrier"],
        "transition_states": ground["transition_states"],
        "gap_coverage": ground["gap_coverage"],
        "thermochemistry": {
            "reactant": ground["reactant"],
            "product": ground["product"],
            "transition_states": ground["transition_states"],
            "rate_limiting_transition_state": ground["rate_limiting_transition_state"],
            "dg_dagger_kcal_mol": ground["dg_dagger_kcal_mol"],
            "dg0_kcal_mol": ground["dg0_kcal_mol"],
            "profile_shape": ground["profile_shape"],
            "has_barrier": ground["has_barrier"],
            "gap_coverage": ground["gap_coverage"],
            "warnings": ground["warnings"],
        },
    }


def write_evb_profile_outputs(analysis: Mapping[str, object], *, output_directory: PathLike) -> None:
    """Write CSV + thermochemistry outputs for an EVB profile analysis.

    Accepts analyses produced by ``analyze_evb_profile_records``: the CSV has
    one column per diabatic state (plus ``evb_ground``) per block, and the
    text summary reports diabatic offsets + coupling matrix.
    """
    outdir = Path(output_directory)
    outdir.mkdir(parents=True, exist_ok=True)

    profile_data = analysis["profile"]
    state_names = tuple(analysis["state_names"])
    curve_names = (*state_names, "evb_ground")
    raw_profiles = profile_data["state_free_energies_kcal_mol"]
    centers = tuple(profile_data["gap_centers_kcal_mol"])
    counts = tuple(profile_data["bin_counts"])

    # The per-replica CSV is a RAW store: alpha/H12-independent diabatic curves +
    # raw gap (the EVB ground column is at the analyser's alpha/H12). Any zero-
    # reference shift is a reporting choice applied downstream from these.
    profile_path = outdir / "evb_profile.csv"
    with profile_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header: list[str] = ["bin_index", "gap_center", "n_points"]
        for name in curve_names:
            header.append(f"{name}_free_energy")
        writer.writerow(header)

        for idx, center in enumerate(centers):
            row: list[object] = [idx, f"{float(center):.8f}", counts[idx]]
            for name in curve_names:
                value = raw_profiles[name][idx]
                row.append("" if value is None else f"{float(value):.8f}")
            writer.writerow(row)

    thermo_path = outdir / "evb_thermochemistry.txt"
    thermo = analysis["thermochemistry"]
    params = analysis["evb_parameters"]
    with thermo_path.open("w", encoding="utf-8") as handle:
        handle.write("N-state EVB profile summary\n")
        handle.write("===========================\n")
        offsets = tuple(float(value) for value in params["diabatic_offsets_kcal_mol"])
        coupling = tuple(tuple(float(value) for value in row) for row in params["coupling_matrix_kcal_mol"])
        rc_states = tuple(int(value) for value in params["reaction_coordinate_states"])
        handle.write(f"states: {', '.join(state_names)}\n")
        handle.write(
            "diabatic offsets (kcal/mol): "
            + ", ".join(f"{name}={offset:.6f}" for name, offset in zip(state_names, offsets))
            + "\n"
        )
        handle.write("coupling matrix (kcal/mol):\n")
        for row in coupling:
            handle.write("  [" + ", ".join(f"{value:.6f}" for value in row) + "]\n")
        handle.write(
            f"reaction coordinate: E[{state_names[rc_states[0]]}] - E[{state_names[rc_states[1]]}]\n\n"
        )

        reactant = thermo["reactant"]
        if reactant is None:
            handle.write("Reactant: not identified\n")
        else:
            handle.write(
                "Reactant: "
                f"bin={int(reactant['index'])}, "
                f"gap={float(reactant['gap_kcal_mol']):.6f}, "
                f"G={float(reactant['free_energy_kcal_mol']):.6f}\n"
            )

        product = thermo["product"]
        if product is None:
            handle.write("Product: not identified\n")
        else:
            handle.write(
                "Product: "
                f"bin={int(product['index'])}, "
                f"gap={float(product['gap_kcal_mol']):.6f}, "
                f"G={float(product['free_energy_kcal_mol']):.6f}\n"
            )

        ts_entries = tuple(thermo["transition_states"])
        if not ts_entries:
            handle.write("Transition states: none\n")
        else:
            handle.write("Transition states:\n")
            for entry in ts_entries:
                handle.write(
                    "  "
                    f"bin={int(entry['index'])}, "
                    f"gap={float(entry['gap_kcal_mol']):.6f}, "
                    f"G={float(entry['free_energy_kcal_mol']):.6f}\n"
                )

        dg_dagger = thermo["dg_dagger_kcal_mol"]
        dg_zero = thermo["dg0_kcal_mol"]
        handle.write("\n")
        handle.write(f"dG‡ (kcal/mol): {'' if dg_dagger is None else f'{float(dg_dagger):.6f}'}\n")
        handle.write(f"dG0 (kcal/mol): {'' if dg_zero is None else f'{float(dg_zero):.6f}'}\n")

        gc = thermo.get("gap_coverage") if isinstance(thermo, Mapping) else None
        if gc is not None:
            sg, rg = gc["sampled_gap_kcal_mol"], gc["resolved_gap_kcal_mol"]
            handle.write(
                f"\ngap coverage: sampled [{sg[0]:.1f}, {sg[1]:.1f}] -> resolved "
                f"[{rg[0]:.1f}, {rg[1]:.1f}] kcal/mol "
                f"({gc['n_bins_resolved']}/{gc['n_bins']} bins; dropped "
                f"{gc['n_bins_dropped_low']} low / {gc['n_bins_dropped_high']} high)\n"
            )
            if gc.get("gap_bounds_rule"):
                handle.write(
                    f"gap bounds: {gc['gap_bounds_rule']}; "
                    f"{gc.get('n_samples_excluded_by_bounds', 0)} sample(s) outside the "
                    f"histogram (free energies still use every sample)\n"
                )
            if gc["truncated_basins"]:
                handle.write(f"** TRUNCATED BASIN(S): {', '.join(gc['truncated_basins'])} "
                             "-- dG may be biased; sample that side more densely/longer **\n")

        warnings_list = tuple(thermo["warnings"])
        if warnings_list:
            handle.write("\nWarnings:\n")
            for warning in warnings_list:
                handle.write(f"- {warning}\n")

    # Sidecar so the RAW evb_profile.csv is self-identifying: a reader knows the
    # diabatic curves are stored RAW (no alpha baked in) and can re-analyse at any
    # (alpha, H12) via load_profile_raw -> analyze_evb_ground without double-counting.
    gap_edges = tuple(float(value) for value in profile_data["gap_edges_kcal_mol"])
    meta = {
        "storage": "raw",
        "reference": params.get("reference", "raw"),
        "state_names": list(state_names),
        "reaction_coordinate_states": [int(value) for value in params["reaction_coordinate_states"]],
        "diabatic_offsets_kcal_mol": [float(value) for value in params["diabatic_offsets_kcal_mol"]],
        "coupling_matrix_kcal_mol": [
            [float(value) for value in row] for row in params["coupling_matrix_kcal_mol"]
        ],
        "coupling_form": params.get("coupling_form", "constant"),
        "n_gap_bins": len(centers),
        "gap_min_kj_mol": float(gap_edges[0]) * KJ_PER_KCAL if gap_edges else 0.0,
        "gap_max_kj_mol": float(gap_edges[-1]) * KJ_PER_KCAL if gap_edges else 0.0,
    }
    (outdir / "evb_profile_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")


def combine_evb_profiles(
    raw_profiles: Sequence[Mapping[str, object]],
    *,
    reference_index: int = 0,
    gap_grid_kcal_mol: Sequence[float] | None = None,
    weight_by_counts: bool = True,
) -> dict[str, object]:
    """Combine several RAW EVB profiles (bidirectional and/or multi-replica) into one.

    Each FEP run fixes its free energy only up to ONE additive constant (its
    reference window), so two runs of the same system have curves that differ by a
    per-run offset even where they sample the same gap. Combination therefore =
    align each run to ``raw_profiles[reference_index]`` by a single constant fit over
    the overlapping gap region (shared across all diabatic states, since the offset
    is one per run), resample onto a common gap grid, then pool.

    Inputs are raw-profile dicts as produced by :func:`build_diabatic_profile` /
    :func:`load_profile_raw`; the result is the same shape (``storage='raw'``) so it
    feeds straight into :func:`analyze_evb_ground`. All inputs must share
    ``state_names`` and ``reaction_coordinate_states``.

    ``weight_by_counts`` weights each profile's contribution to a combined bin by its
    per-bin sample count when available (else equal weights). The default common grid
    spans the union of the inputs' populated gap ranges at the reference's bin spacing.
    """
    profiles = list(raw_profiles)
    if not profiles:
        raise ValueError("combine_evb_profiles requires at least one raw profile.")
    if reference_index < 0 or reference_index >= len(profiles):
        raise ValueError(f"reference_index {reference_index} out of range for {len(profiles)} profiles.")
    ref = profiles[reference_index]
    state_names = tuple(ref["state_names"])
    n_states = len(state_names)
    rcs = tuple(int(v) for v in ref["reaction_coordinate_states"])  # type: ignore[arg-type]

    for idx, prof in enumerate(profiles):
        if tuple(prof["state_names"]) != state_names:
            raise ValueError(
                f"profile {idx} state_names {tuple(prof['state_names'])} != reference {state_names}."
            )
        if tuple(int(v) for v in prof["reaction_coordinate_states"]) != rcs:  # type: ignore[arg-type]
            raise ValueError(
                f"profile {idx} reaction_coordinate_states differ from the reference {rcs}."
            )
        if len(prof["state_free_energies_kj_mol"]) != n_states:  # type: ignore[arg-type]
            raise ValueError(
                f"profile {idx} has {len(prof['state_free_energies_kj_mol'])} state curves; expected {n_states}."
            )

    def _centers(prof: Mapping[str, object]) -> list[float]:
        return [float(x) for x in prof["gap_centers_kcal_mol"]]  # type: ignore[arg-type]

    def _curves_kj(prof: Mapping[str, object]) -> list[list[float | None]]:
        return [[None if v is None else float(v) for v in curve] for curve in prof["state_free_energies_kj_mol"]]  # type: ignore[arg-type]

    # populated gap range (state 0 finite == all states finite, by construction in build_diabatic_profile)
    def _finite_centers(prof: Mapping[str, object]) -> list[float]:
        cen = _centers(prof)
        cur0 = _curves_kj(prof)[0]
        return [cen[i] for i in range(len(cen)) if i < len(cur0) and cur0[i] is not None]

    all_finite = [x for prof in profiles for x in _finite_centers(prof)]
    if not all_finite:
        raise ValueError("none of the profiles has any populated gap bins to combine.")
    gmin, gmax = min(all_finite), max(all_finite)

    # ---- common gap grid (kcal) ----
    if gap_grid_kcal_mol is not None:
        grid = [float(x) for x in gap_grid_kcal_mol]
        spacing = (grid[-1] - grid[0]) / (len(grid) - 1) if len(grid) > 1 else 1.0
    else:
        ref_c = _centers(ref)
        if len(ref_c) >= 2:
            spacing = ref_c[1] - ref_c[0]
        else:
            spacing = next((_centers(p)[1] - _centers(p)[0] for p in profiles if len(_centers(p)) >= 2), 1.0)
        grid = list(ref_c)  # keep the reference centers EXACTLY (so the reference resamples to itself)
        x = ref_c[0] - spacing
        while x >= gmin - abs(spacing) * 0.5:
            grid.insert(0, x); x -= spacing
        x = ref_c[-1] + spacing
        while x <= gmax + abs(spacing) * 0.5:
            grid.append(x); x += spacing
    n_grid = len(grid)
    tol = abs(spacing) * 0.5 + 1e-9

    # ---- resample every profile onto the grid (None outside its own coverage) ----
    gridded_kj: list[list[tuple[float | None, ...]]] = []
    gridded_cnt: list[tuple[float | None, ...]] = []
    gridded_dist: list[tuple[float | None, ...]] = []
    for prof in profiles:
        cen = _centers(prof)
        cur = _curves_kj(prof)
        gridded_kj.append([interpolate_profile_to_grid(cen, cur[k], grid, tolerance=tol) for k in range(n_states)])
        bc = prof.get("bin_counts")
        if bc is not None:
            cnt_curve = [
                (float(bc[i]) if (i < len(cur[0]) and cur[0][i] is not None and i < len(bc)) else None)
                for i in range(len(cen))
            ]
            gridded_cnt.append(interpolate_profile_to_grid(cen, cnt_curve, grid, tolerance=tol))
        else:
            gridded_cnt.append(tuple(None for _ in grid))
        dd = list(prof.get("mean_coupling_distance_nm") or [None] * len(cen))
        dcurve = [(float(dd[i]) if (i < len(dd) and dd[i] is not None) else None) for i in range(len(cen))]
        gridded_dist.append(interpolate_profile_to_grid(cen, dcurve, grid, tolerance=tol))

    # ---- align each profile to the reference by one constant (shared across states) ----
    ref_grid = gridded_kj[reference_index]
    shifts = [0.0] * len(profiles)
    for i in range(len(profiles)):
        if i == reference_index:
            continue
        diffs = [
            ref_grid[k][g] - gridded_kj[i][k][g]
            for k in range(n_states)
            for g in range(n_grid)
            if ref_grid[k][g] is not None and gridded_kj[i][k][g] is not None
        ]
        if not diffs:
            raise ValueError(
                f"profile {i} has no gap overlap with the reference (index {reference_index}); cannot align it."
            )
        shifts[i] = sum(diffs) / len(diffs)
    aligned = [
        [[None if v is None else v + shifts[i] for v in gridded_kj[i][k]] for k in range(n_states)]
        for i in range(len(profiles))
    ]

    # ---- pool ----
    combined_kj: list[list[float | None]] = [[None] * n_grid for _ in range(n_states)]
    combined_counts: list[float] = [0.0] * n_grid
    combined_dist: list[float | None] = [None] * n_grid
    for g in range(n_grid):
        for k in range(n_states):
            num = 0.0
            den = 0.0
            for i in range(len(profiles)):
                v = aligned[i][k][g]
                if v is None:
                    continue
                w = 1.0
                if weight_by_counts:
                    cnt = gridded_cnt[i][g]
                    if cnt is not None and cnt > 0.0:
                        w = float(cnt)
                num += v * w
                den += w
            if den > 0.0:
                combined_kj[k][g] = num / den
        csum = 0.0
        has_cnt = False
        dnum = 0.0
        dden = 0.0
        for i in range(len(profiles)):
            if aligned[i][0][g] is None:
                continue
            cnt = gridded_cnt[i][g]
            if cnt is not None:
                csum += cnt
                has_cnt = True
            dist = gridded_dist[i][g]
            if dist is not None:
                w = float(cnt) if (cnt is not None and cnt > 0.0) else 1.0
                dnum += dist * w
                dden += w
        combined_counts[g] = csum if has_cnt else 0.0
        if dden > 0.0:
            combined_dist[g] = dnum / dden

    combined_kcal = {
        state_names[k]: tuple(None if v is None else kj_to_kcal(v) for v in combined_kj[k])
        for k in range(n_states)
    }
    gap_edges = tuple(grid[0] + (idx - 0.5) * spacing for idx in range(n_grid + 1))
    min_points = max(int(prof.get("min_points_per_bin", 1)) for prof in profiles)  # type: ignore[arg-type]
    return {
        "state_names": state_names,
        "reaction_coordinate_states": (rcs[0], rcs[1]),
        "n_gap_bins": int(n_grid),
        "min_points_per_bin": int(min_points),
        "kT_kcal_mol": float(ref.get("kT_kcal_mol", 0.0)),  # type: ignore[arg-type]
        "gap_min_kj_mol": float((grid[0] - 0.5 * spacing) * KJ_PER_KCAL),
        "gap_max_kj_mol": float((grid[-1] + 0.5 * spacing) * KJ_PER_KCAL),
        "gap_width_kj_mol": float(abs(spacing) * KJ_PER_KCAL),
        "gap_edges_kcal_mol": gap_edges,
        "gap_centers_kcal_mol": tuple(grid),
        "bin_counts": tuple(int(round(c)) for c in combined_counts),
        "state_free_energies_kcal_mol": combined_kcal,
        "state_free_energies_kj_mol": [tuple(curve) for curve in combined_kj],
        "mean_coupling_distance_nm": tuple(combined_dist),
        "alignment_shifts_kj_mol": tuple(shifts),
        "reference_index": int(reference_index),
        "n_combined": len(profiles),
        "storage": "raw",
    }


def load_profile_raw(
    csv_path: PathLike,
    *,
    alpha_run: float | None = None,
    meta_path: PathLike | None = None,
    state_names: Sequence[str] | None = None,
    reaction_coordinate_states: tuple[int, int] = (0, 1),
    min_points_per_bin: int = 1,
) -> dict[str, object]:
    """Load a stored ``evb_profile.csv`` back into the RAW-profile dict shape that
    :func:`analyze_evb_ground` consumes, so it can be re-analysed at any (alpha, H12).

    Two storage conventions are handled:

    * **RAW** (current code) -- detected when a sibling ``evb_profile_meta.json`` (or
      ``meta_path``) reports ``storage == "raw"``. Diabatic columns are used as-is.
    * **OLD / alpha-baked** (pre-refactor CSVs, e.g. existing chorismate runs) -- no
      meta marker. The old analysis baked the run's alpha into storage:
      ``product_free_energy = eps2 + alpha_run`` and ``gap_center = eps1 - eps2 - alpha_run``.
      ``alpha_run`` is then REQUIRED; the product diabat (index ``reaction_coordinate_states[1]``)
      is de-baked by subtracting it and the gap by adding it. Pass ``alpha_run=0.0`` for a
      raw CSV that simply lacks the meta sidecar.

    Old CSVs carry no per-bin coupling distance, so distance-coupling re-analysis is
    unavailable for them (``mean_coupling_distance_nm`` is all ``None``).
    """
    csv_path = Path(csv_path)
    meta: Mapping[str, object] | None = None
    meta_file = Path(meta_path) if meta_path is not None else csv_path.with_name("evb_profile_meta.json")
    if meta_file.exists():
        meta = json.loads(meta_file.read_text(encoding="utf-8"))
    is_raw = bool(meta) and str(meta.get("storage")) == "raw"

    if meta is not None:
        rcs = tuple(int(v) for v in meta.get("reaction_coordinate_states", reaction_coordinate_states))
        meta_state_names = meta.get("state_names")
    else:
        rcs = tuple(int(v) for v in reaction_coordinate_states)
        meta_state_names = None

    with csv_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    # Diabatic state columns = "<state>_free_energy" excluding the derived evb_ground.
    csv_states = [
        name[: -len("_free_energy")]
        for name in fieldnames
        if name.endswith("_free_energy") and name != "evb_ground_free_energy"
    ]
    names = list(state_names) if state_names is not None else (list(meta_state_names) if meta_state_names else csv_states)
    n_states = len(names)
    if n_states < 2:
        raise ValueError(f"{csv_path}: could not identify >=2 diabatic state columns (found {csv_states}).")

    def _f(value: str | None) -> float | None:
        if value is None or value == "" or value == "nan":
            return None
        return float(value)

    gap_centers_kcal = [float(row["gap_center"]) for row in rows]
    curves_kcal: list[list[float | None]] = [[_f(row[f"{name}_free_energy"]) for row in rows] for name in names]

    if not is_raw:
        if alpha_run is None:
            raise ValueError(
                f"{csv_path} has no evb_profile_meta.json with storage='raw' and is treated as an "
                f"OLD alpha-baked profile; pass alpha_run=<the run's alpha_kcal_mol> to de-bake it "
                f"(or alpha_run=0.0 if it is actually raw but lacks the meta sidecar)."
            )
        prod_idx = int(rcs[1])
        curves_kcal[prod_idx] = [None if v is None else v - float(alpha_run) for v in curves_kcal[prod_idx]]
        gap_centers_kcal = [g + float(alpha_run) for g in gap_centers_kcal]

    n_gap_bins = len(rows)
    if meta is not None and "gap_min_kj_mol" in meta:
        gap_min_kj = float(meta["gap_min_kj_mol"])
        gap_max_kj = float(meta["gap_max_kj_mol"])
    elif n_gap_bins >= 2:
        width = gap_centers_kcal[1] - gap_centers_kcal[0]
        gap_min_kj = (gap_centers_kcal[0] - 0.5 * width) * KJ_PER_KCAL
        gap_max_kj = (gap_centers_kcal[-1] + 0.5 * width) * KJ_PER_KCAL
    else:
        gap_min_kj = gap_max_kj = (gap_centers_kcal[0] if gap_centers_kcal else 0.0) * KJ_PER_KCAL

    curves_kj = [
        tuple(None if v is None else float(v) * KJ_PER_KCAL for v in curve) for curve in curves_kcal
    ]
    return {
        "state_names": tuple(names),
        "reaction_coordinate_states": (int(rcs[0]), int(rcs[1])),
        "n_gap_bins": int(n_gap_bins),
        "min_points_per_bin": int(meta.get("min_points_per_bin", min_points_per_bin)) if meta else int(min_points_per_bin),
        "gap_min_kj_mol": float(gap_min_kj),
        "gap_max_kj_mol": float(gap_max_kj),
        "gap_centers_kcal_mol": tuple(gap_centers_kcal),
        "state_free_energies_kj_mol": curves_kj,
        "mean_coupling_distance_nm": tuple([None] * n_gap_bins),
        "storage": "raw",
    }


def aggregate_evb_profiles(
    analyses: Sequence[Mapping[str, object]],
    *,
    minimum_replicas_per_bin: int = 1,
    gap_tolerance: float = 1e-8,
) -> dict[str, object]:
    """Aggregate per-replica EVB profile analyses across an arbitrary number of diabatic states.

    Accepts analyses produced by ``analyze_evb_profile_records``. Per-state
    curves are keyed by state name; the mixed-state curve is keyed ``"evb_ground"``.
    """
    if minimum_replicas_per_bin <= 0:
        raise ValueError("'minimum_replicas_per_bin' must be > 0.")
    if not analyses:
        raise ValueError("At least one replica analysis is required for aggregation.")

    first = analyses[0]
    state_names = tuple(first["state_names"])
    if len(state_names) < 2:
        raise ValueError("Replica aggregation requires at least two EVB states.")

    if str(first.get("energy_unit", REPORTING_ENERGY_UNIT)) != REPORTING_ENERGY_UNIT:
        raise ValueError(f"Replica aggregation expects '{REPORTING_ENERGY_UNIT}' analyses.")

    reference_centers = tuple(float(value) for value in first["profile"]["gap_centers_kcal_mol"])
    n_bins = len(reference_centers)
    if n_bins == 0:
        raise ValueError("Replica profile contains no populated gap bins.")

    profile_blocks = (
        "state_free_energies_kcal_mol",
        "state_free_energies_referenced_kcal_mol",
    )
    curve_names = (*state_names, "evb_ground")

    resampled_profiles: list[dict[str, dict[str, tuple[float | None, ...]]]] = []
    for replica_idx, analysis in enumerate(analyses):
        if tuple(analysis["state_names"]) != state_names:
            raise ValueError(f"Replica {replica_idx} state order does not match the first replica.")
        if str(analysis.get("energy_unit", REPORTING_ENERGY_UNIT)) != REPORTING_ENERGY_UNIT:
            raise ValueError(f"Replica {replica_idx} uses unsupported energy unit '{analysis.get('energy_unit')}'.")
        centers = tuple(float(value) for value in analysis["profile"]["gap_centers_kcal_mol"])
        if len(centers) != n_bins:
            raise ValueError(f"Replica {replica_idx} has {len(centers)} profile bins; expected {n_bins}.")
        profile = analysis["profile"]
        replica_profile: dict[str, dict[str, tuple[float | None, ...]]] = {}
        for block in profile_blocks:
            block_data = profile[block]
            replica_profile[block] = {}
            for curve in curve_names:
                values = tuple(None if value is None else float(value) for value in block_data[curve])
                if len(values) != n_bins:
                    raise ValueError(f"Replica {replica_idx} curve '{curve}' in '{block}' has wrong length.")
                if replica_idx == 0:
                    replica_profile[block][curve] = values
                else:
                    replica_profile[block][curve] = interpolate_profile_to_grid(
                        centers,
                        values,
                        reference_centers,
                        tolerance=float(gap_tolerance),
                    )
        resampled_profiles.append(replica_profile)

    aggregated_profiles: dict[str, dict[str, Mapping[str, tuple[float | None, ...] | tuple[int, ...]]]] = {}
    for block in profile_blocks:
        mean_map: dict[str, tuple[float | None, ...]] = {}
        std_map: dict[str, tuple[float | None, ...]] = {}
        sem_map: dict[str, tuple[float | None, ...]] = {}
        count_map: dict[str, tuple[int, ...]] = {}

        for curve in curve_names:
            mean_values: list[float | None] = []
            std_values: list[float | None] = []
            sem_values: list[float | None] = []
            count_values: list[int] = []

            for bin_idx in range(n_bins):
                replica_values: list[float] = []
                for replica_profile in resampled_profiles:
                    value = replica_profile[block][curve][bin_idx]
                    if value is None:
                        continue
                    replica_values.append(float(value))

                n_rep = len(replica_values)
                count_values.append(n_rep)
                if n_rep < minimum_replicas_per_bin:
                    mean_values.append(None)
                    std_values.append(None)
                    sem_values.append(None)
                    continue

                mean, std, sem = mean_std_sem(replica_values)
                mean_values.append(mean)
                std_values.append(std)
                sem_values.append(sem)

            mean_map[curve] = tuple(mean_values)
            std_map[curve] = tuple(std_values)
            sem_map[curve] = tuple(sem_values)
            count_map[curve] = tuple(count_values)

        aggregated_profiles[block] = {
            "mean_kcal_mol": MappingProxyType(mean_map),
            "std_kcal_mol": MappingProxyType(std_map),
            "sem_kcal_mol": MappingProxyType(sem_map),
            "n_replicas": MappingProxyType(count_map),
        }

    replica_summaries: list[dict[str, object]] = []
    dg_dagger_values: list[float] = []
    dg0_values: list[float] = []
    warning_messages: list[str] = []

    for idx, analysis in enumerate(analyses):
        thermo = analysis["thermochemistry"]
        dg_dagger = thermo.get("dg_dagger_kcal_mol")
        dg0 = thermo.get("dg0_kcal_mol")
        warnings_tuple = tuple(thermo.get("warnings", tuple()))
        if dg_dagger is not None:
            dg_dagger_values.append(float(dg_dagger))
        if dg0 is not None:
            dg0_values.append(float(dg0))
        warning_messages.extend(str(value) for value in warnings_tuple)
        replica_summaries.append(
            {
                "replica_index": int(idx),
                "dg_dagger_kcal_mol": None if dg_dagger is None else float(dg_dagger),
                "dg0_kcal_mol": None if dg0 is None else float(dg0),
                "warnings": warnings_tuple,
            }
        )

    dg_dagger_mean = dg_dagger_std = dg_dagger_sem = None
    if dg_dagger_values:
        dg_dagger_mean, dg_dagger_std, dg_dagger_sem = mean_std_sem(dg_dagger_values)

    dg0_mean = dg0_std = dg0_sem = None
    if dg0_values:
        dg0_mean, dg0_std, dg0_sem = mean_std_sem(dg0_values)

    return {
        "state_names": state_names,
        "energy_unit": REPORTING_ENERGY_UNIT,
        "replica_count": len(analyses),
        "minimum_replicas_per_bin": int(minimum_replicas_per_bin),
        "gap_centers_kcal_mol": reference_centers,
        "profiles": MappingProxyType(aggregated_profiles),
        "thermochemistry": {
            "dg_dagger_values_kcal_mol": tuple(dg_dagger_values),
            "dg_dagger_mean_kcal_mol": dg_dagger_mean,
            "dg_dagger_std_kcal_mol": dg_dagger_std,
            "dg_dagger_sem_kcal_mol": dg_dagger_sem,
            "dg0_values_kcal_mol": tuple(dg0_values),
            "dg0_mean_kcal_mol": dg0_mean,
            "dg0_std_kcal_mol": dg0_std,
            "dg0_sem_kcal_mol": dg0_sem,
            "warnings": tuple(warning_messages),
        },
        "replica_summaries": tuple(replica_summaries),
    }


def write_evb_replica_outputs(aggregation: Mapping[str, object], *, output_directory: PathLike) -> None:
    """Write replica ensemble outputs (CSV + thermochemistry summary).

    Emits one column per diabatic state (plus ``evb_ground``) per block per
    statistic, deriving the curve list from ``aggregation["state_names"]``.
    """
    outdir = Path(output_directory)
    outdir.mkdir(parents=True, exist_ok=True)

    state_names = tuple(aggregation["state_names"])
    curve_names = (*state_names, "evb_ground")
    centers = tuple(aggregation["gap_centers_kcal_mol"])
    profiles = aggregation["profiles"]
    block_alias = (
        ("state_free_energies_kcal_mol", "raw"),
        ("state_free_energies_referenced_kcal_mol", "referenced"),
    )

    csv_path = outdir / "evb_replica_profile_ensemble.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header: list[str] = ["bin_index", "gap_center"]
        for curve in curve_names:
            for _, alias in block_alias:
                header.extend(
                    [
                        f"{curve}_{alias}_mean",
                        f"{curve}_{alias}_std",
                        f"{curve}_{alias}_sem",
                        f"{curve}_{alias}_n",
                    ]
                )
        writer.writerow(header)

        for idx, center in enumerate(centers):
            row: list[object] = [idx, f"{float(center):.8f}"]
            for curve in curve_names:
                for block, _ in block_alias:
                    block_data = profiles[block]
                    mean = block_data["mean_kcal_mol"][curve][idx]
                    std = block_data["std_kcal_mol"][curve][idx]
                    sem = block_data["sem_kcal_mol"][curve][idx]
                    count = int(block_data["n_replicas"][curve][idx])
                    row.append("" if mean is None else f"{float(mean):.8f}")
                    row.append("" if std is None else f"{float(std):.8f}")
                    row.append("" if sem is None else f"{float(sem):.8f}")
                    row.append(count)
            writer.writerow(row)

    thermo = aggregation["thermochemistry"]
    txt_path = outdir / "evb_replica_thermochemistry.txt"
    with txt_path.open("w", encoding="utf-8") as handle:
        def _fmt_optional(value: float | None) -> str:
            if value is None:
                return ""
            return f"{float(value):.6f}"

        handle.write("EVB replica thermochemistry summary\n")
        handle.write("==================================\n")
        handle.write(f"states: {', '.join(state_names)}\n")
        handle.write(f"replicas: {int(aggregation['replica_count'])}\n")
        handle.write(f"minimum replicas per profile bin: {int(aggregation['minimum_replicas_per_bin'])}\n\n")
        handle.write(
            "dG‡ (kcal/mol): "
            f"mean={_fmt_optional(thermo['dg_dagger_mean_kcal_mol'])}, "
            f"std={_fmt_optional(thermo['dg_dagger_std_kcal_mol'])}, "
            f"sem={_fmt_optional(thermo['dg_dagger_sem_kcal_mol'])}\n"
        )
        handle.write(
            "dG0 (kcal/mol): "
            f"mean={_fmt_optional(thermo['dg0_mean_kcal_mol'])}, "
            f"std={_fmt_optional(thermo['dg0_std_kcal_mol'])}, "
            f"sem={_fmt_optional(thermo['dg0_sem_kcal_mol'])}\n"
        )

        warning_list = tuple(thermo.get("warnings", tuple()))
        if warning_list:
            handle.write("\nWarnings:\n")
            for warning in warning_list:
                handle.write(f"- {warning}\n")


__all__ = [
    "aggregate_evb_profiles",
    "analyze_evb_profile_records",
    "analyze_evb_profile_from_raw",
    "build_diabatic_profile",
    "analyze_evb_ground",
    "load_profile_raw",
    "combine_evb_profiles",
    "resolve_curve_reference",
    "detect_profile_thermochemistry",
    "interpolate_profile_to_grid",
    "write_evb_profile_outputs",
    "write_evb_replica_outputs",
]
