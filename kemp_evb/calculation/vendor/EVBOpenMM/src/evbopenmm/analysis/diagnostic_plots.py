"""Matplotlib render helpers for per-replica EVB diagnostic panels.

Each helper consumes one of the frozen records in
:mod:`evbopenmm.analysis.diagnostic_records` (no recomputation, no tie to
any specific engine or simulation output) and draws onto a Matplotlib
``Axes``. Follow the same ``ax=None → auto-create`` convention as other
helpers in this package:

* Pass ``ax`` to draw onto an existing axis (returns that axis).
* Pass ``ax=None`` to let the helper create a fresh ``Figure``/``Axes``.

:func:`render_diagnostic_grid` assembles a five-panel per-replica health
check figure from a single ``analyze_evb_profile_records`` result — the
EVB-literature-standard diagnostic view.

Matplotlib is imported lazily inside each entry point so importing this
module stays cheap and keeps Matplotlib a soft runtime dependency.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Mapping

import numpy as np

from .diagnostic_records import (
    ForwardReverseResidual,
    ProfileCurve,
    StateEnergyVsLambda,
    WindowOverlap,
)


def _import_mpl():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - runtime dependency
        raise ImportError(
            "matplotlib is required for diagnostic plots. "
            "Install it with `pip install matplotlib`."
        ) from exc
    return plt


def _nan_fill(values):
    return [float("nan") if v is None else float(v) for v in values]


def plot_profile_curve(
    record: ProfileCurve,
    *,
    ax=None,
    label: str | None = None,
    **plot_kwargs: Any,
):
    """Plot one state's free-energy profile vs energy-gap reaction coordinate."""
    plt = _import_mpl()
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(
        record.gap_centers_kcal_mol,
        _nan_fill(record.free_energy_kcal_mol),
        label=label or record.state_name,
        **plot_kwargs,
    )
    ax.set_xlabel(r"Energy gap $\Delta E$ (kcal/mol)")
    ax.set_ylabel(f"Free energy (kcal/mol, {record.normalization}-norm)")
    ax.set_title(f"Free-energy profile: {record.state_name}")
    return ax


def plot_gap_histogram(
    record: ProfileCurve,
    *,
    ax=None,
    **bar_kwargs: Any,
):
    """Plot per-bin sample counts along the gap reaction coordinate."""
    plt = _import_mpl()
    if ax is None:
        _, ax = plt.subplots()
    centers = record.gap_centers_kcal_mol
    width = centers[1] - centers[0] if len(centers) >= 2 else 1.0
    ax.bar(centers, record.bin_counts, width=width, **bar_kwargs)
    ax.set_xlabel(r"Energy gap $\Delta E$ (kcal/mol)")
    ax.set_ylabel("Samples per bin")
    ax.set_title("Gap histogram")
    return ax


def plot_window_overlap(
    record: WindowOverlap,
    *,
    ax=None,
    **plot_kwargs: Any,
):
    """Plot per-λ-window mean energy gap — the binning-overlap diagnostic."""
    plt = _import_mpl()
    if ax is None:
        _, ax = plt.subplots()
    plot_kwargs.setdefault("marker", "o")
    ax.plot(record.lambda_values, record.mean_gap_kcal_mol, **plot_kwargs)
    ax.set_xlabel(f"λ ({record.reference_state_name})")
    ax.set_ylabel(r"$\langle \Delta E \rangle$ per window (kcal/mol)")
    ax.set_title("Window overlap")
    return ax


def plot_forward_reverse_residual(
    record: ForwardReverseResidual,
    *,
    ax=None,
    **plot_kwargs: Any,
):
    """Plot the per-window ``ΔGf − ΔGr`` residual (windows with poor overlap spike)."""
    plt = _import_mpl()
    if ax is None:
        _, ax = plt.subplots()
    plot_kwargs.setdefault("marker", "o")
    ax.axhline(0.0, color="0.5", linewidth=0.5)
    ax.plot(record.window_index, record.residual_kcal_mol, **plot_kwargs)
    ax.set_xlabel("Window edge index")
    ax.set_ylabel(r"$\Delta G_f - \Delta G_r$ (kcal/mol)")
    ax.set_title("Forward/reverse residual")
    return ax


def plot_state_energy_vs_lambda(
    record: StateEnergyVsLambda,
    *,
    ax=None,
    state_names: Sequence[str] | None = None,
):
    """Plot per-state mean potential energy across the λ ladder."""
    plt = _import_mpl()
    if ax is None:
        _, ax = plt.subplots()
    requested = tuple(state_names) if state_names is not None else record.state_names
    for name in requested:
        if name not in record.mean_energies_kcal_mol:
            raise ValueError(
                f"state_name '{name}' not present; available: {list(record.state_names)}."
            )
        ax.plot(
            record.lambda_values,
            record.mean_energies_kcal_mol[name],
            marker="o",
            label=name,
        )
    ax.set_xlabel(f"λ ({record.reference_state_name})")
    ax.set_ylabel(r"$\langle E_i \rangle$ (kcal/mol)")
    ax.set_title("Per-state mean energy vs λ")
    ax.legend(loc="best")
    return ax


def plot_replica_profile_overlay(
    records: Sequence[ProfileCurve],
    *,
    ax=None,
    show_mean: bool = True,
    alpha_replicas: float = 0.35,
    replica_color: str = "steelblue",
    mean_color: str = "black",
):
    """Overlay multiple replica ``ProfileCurve`` records with optional mean ± std.

    The mean and ±1σ band are computed per gap bin across replicas, ignoring
    bins where any replica is missing data (``nanmean``/``nanstd``).
    """
    plt = _import_mpl()
    if ax is None:
        _, ax = plt.subplots()
    if not records:
        raise ValueError("records must be non-empty.")
    ref_x = records[0].gap_centers_kcal_mol
    for rec in records[1:]:
        if rec.gap_centers_kcal_mol != ref_x:
            raise ValueError(
                "All records must share identical gap_centers_kcal_mol to overlay."
            )

    matrix = np.array(
        [_nan_fill(rec.free_energy_kcal_mol) for rec in records],
        dtype=float,
    )
    for row in matrix:
        ax.plot(ref_x, row, color=replica_color, alpha=alpha_replicas, linewidth=1.0)

    if show_mean:
        with np.errstate(invalid="ignore"):
            mean = np.nanmean(matrix, axis=0)
            std = np.nanstd(matrix, axis=0)
        ax.plot(
            ref_x,
            mean,
            color=mean_color,
            linewidth=1.5,
            label=f"mean (N={len(records)})",
        )
        ax.fill_between(ref_x, mean - std, mean + std, color=mean_color, alpha=0.15)
        ax.legend(loc="best")

    ax.set_xlabel(r"Energy gap $\Delta E$ (kcal/mol)")
    ax.set_ylabel(f"Free energy (kcal/mol, {records[0].normalization}-norm)")
    ax.set_title(f"Replica overlay: {records[0].state_name}")
    return ax


def render_diagnostic_grid(
    analysis: Mapping[str, Any],
    *,
    figure=None,
):
    """Render a five-panel diagnostic grid from one profile analysis dict.

    Layout (3 rows × 2 cols, bottom row spans):

    ======================================
    |  Profile curve  |  Gap histogram   |
    |  Window overlap |  Fwd/rev resid   |
    |       Per-state energy vs λ        |
    ======================================
    """
    plt = _import_mpl()
    profile = ProfileCurve.from_analysis(analysis)
    overlap = WindowOverlap.from_analysis(analysis)
    fwd_rev = ForwardReverseResidual.from_analysis(analysis)
    state_e = StateEnergyVsLambda.from_analysis(analysis)

    if figure is None:
        figure = plt.figure(figsize=(12, 9))
    else:
        figure.clear()

    gs = figure.add_gridspec(3, 2)
    ax_profile = figure.add_subplot(gs[0, 0])
    ax_hist = figure.add_subplot(gs[0, 1])
    ax_overlap = figure.add_subplot(gs[1, 0])
    ax_residual = figure.add_subplot(gs[1, 1])
    ax_state = figure.add_subplot(gs[2, :])

    plot_profile_curve(profile, ax=ax_profile)
    plot_gap_histogram(profile, ax=ax_hist)
    plot_window_overlap(overlap, ax=ax_overlap)
    plot_forward_reverse_residual(fwd_rev, ax=ax_residual)
    plot_state_energy_vs_lambda(state_e, ax=ax_state)

    figure.tight_layout()
    return figure


__all__ = [
    "plot_forward_reverse_residual",
    "plot_gap_histogram",
    "plot_profile_curve",
    "plot_replica_profile_overlay",
    "plot_state_energy_vs_lambda",
    "plot_window_overlap",
    "render_diagnostic_grid",
]
