"""Reaction-coordinate (gap) histogram bounds for umbrella/EVB profile analysis.

Both the mapping-FEP umbrella analysis and the EVB ground-profile analysis bin samples along
the diabatic energy gap. By default the histogram spans the full sampled range, which is not
robust: the energy of a diabatic state evaluated on configurations sampled by a *different*
state can have a heavy tail whenever the two states differ in their nonbonded parameters (a
transferring atom that is a nonpolar hydrogen in one state and a polar hydrogen with no
Lennard-Jones sphere in the other is the usual case). A handful of such samples stretches the
range until nearly every bin is empty and no basin can be resolved.

This module resolves the bounds used for that histogram. It deliberately does **not** touch the
free-energy estimator: the forward/reverse exponential averages between windows must be taken
over every sample, so restricting the range affects binning only.
"""

from __future__ import annotations

from collections.abc import Sequence

__all__ = ["percentile", "resolve_gap_bounds"]


def percentile(sorted_values: Sequence[float], fraction: float) -> float:
    """Linear-interpolated percentile of an already-sorted sequence.

    ``fraction`` is in percent (0-100). Implemented locally to keep this module dependency-free.
    """
    if not sorted_values:
        raise ValueError("Cannot take a percentile of an empty sequence.")
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    position = (float(fraction) / 100.0) * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = position - lower
    return float(sorted_values[lower]) * (1.0 - weight) + float(sorted_values[upper]) * weight


def resolve_gap_bounds(
    gap_values: Sequence[float],
    *,
    gap_range: tuple[float, float] | Sequence[float] | None = None,
    gap_percentile_range: tuple[float, float] | Sequence[float] | None = None,
    scale_to_internal: float = 1.0,
) -> tuple[float, float, int, str]:
    """Resolve the histogram bounds for a set of reaction-coordinate values.

    Parameters
    ----------
    gap_values :
        All sampled gap values, already converted to the internal energy unit.
    gap_range :
        Explicit ``(low, high)`` bounds in the *input* energy unit (they are multiplied by
        ``scale_to_internal``). Mutually exclusive with ``gap_percentile_range``.
    gap_percentile_range :
        Robust ``(low, high)`` percentiles in 0-100, e.g. ``(1.0, 99.0)``. Mutually exclusive
        with ``gap_range``. **The bound must undercut the outlier fraction**: percentiles are
        linearly interpolated, so with 1% of samples in the tail a 99.5th percentile still lands
        inside it and clips nothing useful. Choose the bound from how many samples you expect to
        be pathological, not from habit.
    scale_to_internal :
        Factor converting the input energy unit to the internal one, applied to ``gap_range``.

    Returns
    -------
    tuple
        ``(low, high, n_excluded, rule)`` -- the bounds in the internal unit, how many samples
        fall outside them, and a short description of the rule applied.

    Notes
    -----
    Passing neither argument reproduces the historical behaviour exactly: the full sampled
    range, with nothing excluded.
    """
    if not gap_values:
        raise ValueError("No gap values were supplied; cannot resolve histogram bounds.")
    if gap_range is not None and gap_percentile_range is not None:
        raise ValueError("Specify only one of 'gap_range' or 'gap_percentile_range'.")

    if gap_range is not None:
        bounds = tuple(float(value) for value in gap_range)
        if len(bounds) != 2:
            raise ValueError("'gap_range' must be a (low, high) pair.")
        low, high = bounds[0] * scale_to_internal, bounds[1] * scale_to_internal
        if not (low < high):
            raise ValueError(f"'gap_range' requires low < high; received {bounds}.")
        if any(value != value or value in (float("inf"), float("-inf")) for value in (low, high)):
            raise ValueError("'gap_range' bounds must be finite.")
        rule = f"explicit range [{bounds[0]:.6g}, {bounds[1]:.6g}]"
    elif gap_percentile_range is not None:
        bounds = tuple(float(value) for value in gap_percentile_range)
        if len(bounds) != 2:
            raise ValueError("'gap_percentile_range' must be a (low, high) pair.")
        if not (0.0 <= bounds[0] < bounds[1] <= 100.0):
            raise ValueError(
                f"'gap_percentile_range' requires 0 <= low < high <= 100; received {bounds}."
            )
        ordered = sorted(float(value) for value in gap_values)
        low = percentile(ordered, bounds[0])
        high = percentile(ordered, bounds[1])
        rule = f"percentile range [{bounds[0]:.6g}%, {bounds[1]:.6g}%]"
    else:
        low = min(gap_values)
        high = max(gap_values)
        rule = "full sampled range"

    if abs(high - low) <= 1e-12:
        high = low + 1e-6

    n_excluded = sum(1 for value in gap_values if value < low or value > high)
    if n_excluded == len(gap_values):
        raise ValueError(
            "The requested gap bounds exclude every sample; widen 'gap_range' / "
            "'gap_percentile_range'."
        )
    return low, high, n_excluded, rule
