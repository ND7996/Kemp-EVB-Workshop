"""Plotting helpers that sit on top of :class:`~evbopenmm.system.EVBSystem`.

These are diagnostic/visualization utilities kept out of ``system.py`` because
they are pure helpers around an :class:`EVBSystem` instance — they do not
mutate state and they depend on optional plotting/analysis libraries.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


__all__ = [
    "plot_evb_profile_replicas",
]


def plot_evb_profile_replicas(
    replica_result: Mapping[str, object],
    *,
    profile_block: str = "state_free_energies_kcal_mol",
    curves: Sequence[str] | None = None,
    include_aggregate_mean: bool = True,
    replica_alpha: float = 0.45,
    replica_linewidth: float = 1.1,
    mean_linewidth: float = 2.2,
    figure_size: Tuple[float, float] | None = None,
    title: str | None = None,
    show_legend: bool = True,
) -> Dict[str, object]:
    """Plot per-replica EVB profile curves for an arbitrary number of diabatic states.

    Accepts the output of :func:`analyze_evb_profile_replicas_result`. Per-state
    curves are keyed by state name; the mixed curve is keyed ``"evb_ground"``.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "matplotlib is required for plotting replica curves. Install matplotlib to use this helper."
        ) from exc

    replicas = tuple(replica_result.get("replicas", tuple()))
    if not replicas:
        raise ValueError("'replica_result' does not contain replica analyses.")

    first = replicas[0]
    state_names = tuple(first["state_names"])
    if len(state_names) < 2:
        raise ValueError("Replica plotting requires at least two EVB states.")

    valid_profile_blocks = (
        "state_free_energies_kcal_mol",
        "state_free_energies_referenced_kcal_mol",
    )
    if profile_block not in valid_profile_blocks:
        allowed = ", ".join(valid_profile_blocks)
        raise ValueError(f"Unsupported profile_block '{profile_block}'. Allowed values: {allowed}.")

    default_curves = (*state_names, "evb_ground")
    curve_names = tuple(default_curves if curves is None else curves)
    if not curve_names:
        raise ValueError("'curves' must contain at least one curve name.")
    for curve in curve_names:
        if curve not in default_curves:
            allowed_curves = ", ".join(default_curves)
            raise ValueError(f"Unsupported curve '{curve}'. Allowed values: {allowed_curves}.")

    n_curves = len(curve_names)
    if figure_size is None:
        figure_size = (7.0, 3.2 * float(n_curves))
    fig, axes = plt.subplots(n_curves, 1, figsize=figure_size, sharex=True, squeeze=False)
    flat_axes = tuple(axes[:, 0])

    n_replicas = len(replicas)
    color_map = plt.get_cmap("viridis")
    curve_labels: dict[str, str] = {name: f"{name} diabatic" for name in state_names}
    curve_labels["evb_ground"] = "EVB ground"

    for axis, curve_name in zip(flat_axes, curve_names):
        for replica_idx, analysis in enumerate(replicas):
            profile_data = analysis["profile"]
            gap_centers = tuple(float(value) for value in profile_data["gap_centers_kcal_mol"])
            values_raw = profile_data[profile_block][curve_name]
            values = tuple(math.nan if value is None else float(value) for value in values_raw)

            if n_replicas == 1:
                color_value = color_map(0.55)
            else:
                color_value = color_map(float(replica_idx) / float(max(n_replicas - 1, 1)))

            line_label = None
            if show_legend and n_replicas <= 12:
                line_label = f"replica {replica_idx + 1}"

            axis.plot(
                gap_centers,
                values,
                color=color_value,
                alpha=float(replica_alpha),
                linewidth=float(replica_linewidth),
                label=line_label,
            )

        if include_aggregate_mean:
            aggregate = replica_result.get("aggregate")
            if aggregate is None:
                raise ValueError("Replica aggregate data is required when include_aggregate_mean=True.")
            aggregate_profiles = aggregate["profiles"][profile_block]["mean_kcal_mol"]
            aggregate_centers = tuple(float(value) for value in aggregate["gap_centers_kcal_mol"])
            aggregate_values = tuple(
                math.nan if value is None else float(value)
                for value in aggregate_profiles[curve_name]
            )
            axis.plot(
                aggregate_centers,
                aggregate_values,
                color="black",
                linewidth=float(mean_linewidth),
                label="replica mean" if show_legend else None,
            )

        axis.set_ylabel("dG (kcal/mol)")
        axis.set_title(curve_labels[curve_name])
        axis.grid(alpha=0.3)
        if show_legend:
            axis.legend(loc="best", fontsize=8, ncol=2)

    flat_axes[-1].set_xlabel("Reaction coordinate (kcal/mol)")
    if title:
        fig.suptitle(str(title))
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    else:
        fig.tight_layout()

    return {
        "figure": fig,
        "axes": flat_axes,
        "curve_names": curve_names,
    }
