"""Campaign-level analyzer for multi-environment, multi-parameter EVB studies.

A *campaign* is a 2-D matrix of cells indexed by ``(environment, parameter_lane)``
— for example ``("gas", "q6_legacy")``, ``("aqueous", "openff")``, etc. Each
cell carries ``N`` replica profile-analysis dicts (as produced by
:func:`analyze_evb_profile_records`), one per independent trajectory.

:class:`CampaignAnalyzer` is a thin, property-based convenience layer on top of
the pure-function analysis core. Every method either:

* delegates to one of the diagnostic-record factories or plot helpers, or
* aggregates per-replica statistics (ΔG‡ / ΔG₀ mean and SEM).

Caching is explicit and per-cell — recomputation costs are modest
(dataclasses from already-computed dicts) but the typical notebook workflow
re-reads the same cells repeatedly, so one-time construction per cell is
enough. No analysis logic lives inside the class — just plumbing.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import numpy as np

from .diagnostic_plots import (
    plot_replica_profile_overlay,
    render_diagnostic_grid,
)
from .diagnostic_records import (
    ForwardReverseResidual,
    ProfileCurve,
    StateEnergyVsLambda,
    WindowOverlap,
)

CellKey = tuple[str, str]


@dataclass(frozen=True)
class ReplicaDiagnostics:
    """The four diagnostic records for one replica's analysis, bundled."""

    profile: ProfileCurve
    window_overlap: WindowOverlap
    forward_reverse: ForwardReverseResidual
    state_energy: StateEnergyVsLambda


@dataclass(frozen=True)
class CellStatistics:
    """Mean/SEM statistics across replicas for a single cell."""

    n_replicas: int
    dg_dagger_mean_kcal_mol: float
    dg_dagger_sem_kcal_mol: float
    dg0_mean_kcal_mol: float
    dg0_sem_kcal_mol: float


def _mean_sem(values: Sequence[float | None]) -> tuple[float, float]:
    valid = np.array([v for v in values if v is not None], dtype=float)
    if valid.size == 0:
        return float("nan"), float("nan")
    mean = float(valid.mean())
    if valid.size < 2:
        return mean, float("nan")
    sem = float(valid.std(ddof=1) / np.sqrt(valid.size))
    return mean, sem


class CampaignAnalyzer:
    """Thin property-based analyzer for a multi-cell EVB campaign.

    Construct with a mapping ``{(env, lane): [analysis_dict, ...]}``. Each
    analysis dict must be the output of
    :func:`analyze_evb_profile_records` for one replica. Environments and
    parameter lanes are inferred from the keys; order is preserved.
    """

    def __init__(
        self,
        cells: Mapping[CellKey, Sequence[Mapping[str, Any]]],
    ) -> None:
        if not cells:
            raise ValueError("cells must contain at least one (env, lane) entry.")
        normalised: dict[CellKey, tuple[Mapping[str, Any], ...]] = {}
        envs_in_order: list[str] = []
        lanes_in_order: list[str] = []
        for key, analyses in cells.items():
            if not isinstance(key, tuple) or len(key) != 2:
                raise ValueError(
                    f"Cell key {key!r} must be a 2-tuple (environment, parameter_lane)."
                )
            env, lane = str(key[0]), str(key[1])
            if not analyses:
                raise ValueError(f"Cell ({env!r}, {lane!r}) has no replica analyses.")
            normalised[(env, lane)] = tuple(analyses)
            if env not in envs_in_order:
                envs_in_order.append(env)
            if lane not in lanes_in_order:
                lanes_in_order.append(lane)
        self._cells: Mapping[CellKey, tuple[Mapping[str, Any], ...]] = normalised
        self._environments: tuple[str, ...] = tuple(envs_in_order)
        self._parameter_lanes: tuple[str, ...] = tuple(lanes_in_order)

    # ------------------------------------------------------------------ descriptors

    @property
    def environments(self) -> tuple[str, ...]:
        return self._environments

    @property
    def parameter_lanes(self) -> tuple[str, ...]:
        return self._parameter_lanes

    @property
    def cells(self) -> tuple[CellKey, ...]:
        return tuple(self._cells)

    def n_replicas(self, environment: str, parameter_lane: str) -> int:
        return len(self._require_cell(environment, parameter_lane))

    # ------------------------------------------------------------------ extraction

    def replica_profiles(
        self,
        environment: str,
        parameter_lane: str,
        *,
        state_name: str = "evb_ground",
        normalization: str = "raw",
    ) -> tuple[ProfileCurve, ...]:
        """One ``ProfileCurve`` per replica for the given cell."""
        return self._cached_profiles(
            environment, parameter_lane, state_name, normalization
        )

    @lru_cache(maxsize=None)
    def _cached_profiles(
        self,
        environment: str,
        parameter_lane: str,
        state_name: str,
        normalization: str,
    ) -> tuple[ProfileCurve, ...]:
        analyses = self._require_cell(environment, parameter_lane)
        return tuple(
            ProfileCurve.from_analysis(
                a, state_name=state_name, normalization=normalization  # type: ignore[arg-type]
            )
            for a in analyses
        )

    def replica_diagnostics(
        self,
        environment: str,
        parameter_lane: str,
    ) -> tuple[ReplicaDiagnostics, ...]:
        """One :class:`ReplicaDiagnostics` bundle per replica."""
        return self._cached_diagnostics(environment, parameter_lane)

    @lru_cache(maxsize=None)
    def _cached_diagnostics(
        self, environment: str, parameter_lane: str
    ) -> tuple[ReplicaDiagnostics, ...]:
        analyses = self._require_cell(environment, parameter_lane)
        return tuple(
            ReplicaDiagnostics(
                profile=ProfileCurve.from_analysis(a),
                window_overlap=WindowOverlap.from_analysis(a),
                forward_reverse=ForwardReverseResidual.from_analysis(a),
                state_energy=StateEnergyVsLambda.from_analysis(a),
            )
            for a in analyses
        )

    def cell_statistics(
        self, environment: str, parameter_lane: str
    ) -> CellStatistics:
        """Mean and SEM of ΔG‡ and ΔG₀ across the cell's replicas."""
        analyses = self._require_cell(environment, parameter_lane)
        dg_daggers = tuple(
            a["thermochemistry"].get("dg_dagger_kcal_mol") for a in analyses
        )
        dg0s = tuple(a["thermochemistry"].get("dg0_kcal_mol") for a in analyses)
        dg_dagger_mean, dg_dagger_sem = _mean_sem(dg_daggers)
        dg0_mean, dg0_sem = _mean_sem(dg0s)
        return CellStatistics(
            n_replicas=len(analyses),
            dg_dagger_mean_kcal_mol=dg_dagger_mean,
            dg_dagger_sem_kcal_mol=dg_dagger_sem,
            dg0_mean_kcal_mol=dg0_mean,
            dg0_sem_kcal_mol=dg0_sem,
        )

    # ------------------------------------------------------------------ plot helpers

    def plot_replica_overlay(
        self,
        environment: str,
        parameter_lane: str,
        *,
        ax=None,
        show_mean: bool = True,
        state_name: str = "evb_ground",
        normalization: str = "raw",
    ):
        """Overlay all replicas' free-energy profiles for the given cell."""
        records = self.replica_profiles(
            environment, parameter_lane,
            state_name=state_name, normalization=normalization,
        )
        ax = plot_replica_profile_overlay(records, ax=ax, show_mean=show_mean)
        ax.set_title(
            f"Replica overlay — env={environment!r}, lane={parameter_lane!r} "
            f"(N={len(records)})"
        )
        return ax

    def plot_diagnostic_grid(
        self,
        environment: str,
        parameter_lane: str,
        *,
        replica_index: int = 0,
        figure=None,
    ):
        """Render a Fig-10-style 5-panel grid for one replica of the cell."""
        analyses = self._require_cell(environment, parameter_lane)
        if not 0 <= int(replica_index) < len(analyses):
            raise ValueError(
                f"replica_index={replica_index} out of range "
                f"[0, {len(analyses) - 1}] for cell ({environment!r}, {parameter_lane!r})."
            )
        return render_diagnostic_grid(analyses[int(replica_index)], figure=figure)

    # ------------------------------------------------------------------ helpers

    def _require_cell(
        self, environment: str, parameter_lane: str
    ) -> tuple[Mapping[str, Any], ...]:
        try:
            return self._cells[(environment, parameter_lane)]
        except KeyError as exc:
            raise KeyError(
                f"No cell for (environment={environment!r}, "
                f"parameter_lane={parameter_lane!r}); "
                f"available cells: {list(self._cells)}."
            ) from exc


__all__ = [
    "CampaignAnalyzer",
    "CellKey",
    "CellStatistics",
    "ReplicaDiagnostics",
]
