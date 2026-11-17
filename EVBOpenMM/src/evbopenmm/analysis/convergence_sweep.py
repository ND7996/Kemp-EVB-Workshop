"""Convergence-sweep tooling for EVB studies.

The EVB literature (e.g. Purg & Kamerlin, *Methods in Enzymology* 2018 Fig 9)
recommends establishing a set of convergence curves before production runs:

* ``ΔG‡`` vs sampling time per λ window (pick a plateau);
* ``ΔG‡`` vs equilibration time;
* ``ΔG‡`` vs position-restraint strength on reactive atoms;
* standard-error-of-the-mean vs number of replicas (to pick how many to
  run, given a target statistical precision).

This module offers two complementary primitives:

* :class:`ConvergenceSweepRecord` — a 1-D sweep of any observable across any
  control knob. One ``x_value`` per point, each with a vector of per-replica
  observable values. Exposes ``means``/``stds``/``sems`` as derived
  properties. Agnostic about what's being varied.
* :class:`SEMvsReplicateCountRecord` — the empirical SEM → N_replicas curve
  (Fig 9 panel D). Bootstrap subsampling of a fixed pool of replicas.

Plus paired plot helpers, following the same ``ax=None → auto-create``
convention used elsewhere in ``evbopenmm.analysis``.
"""
from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ConvergenceSweepRecord:
    """1-D sweep of an observable (e.g. ΔG‡) vs a control knob (e.g. ps/window).

    Per-point replica values are allowed to contain ``None`` to mark failed
    replicas; they are excluded from the derived statistics.
    """

    x_label: str
    x_values: tuple[float, ...]
    observable_label: str
    replica_values_per_point: tuple[tuple[float | None, ...], ...]

    def __post_init__(self) -> None:
        if len(self.x_values) != len(self.replica_values_per_point):
            raise ValueError(
                f"x_values has length {len(self.x_values)} but "
                f"replica_values_per_point has length "
                f"{len(self.replica_values_per_point)}."
            )

    def _valid_per_point(self) -> list[np.ndarray]:
        return [
            np.array(
                [float(v) for v in replicas if v is not None], dtype=float
            )
            for replicas in self.replica_values_per_point
        ]

    @property
    def means(self) -> tuple[float, ...]:
        return tuple(
            float("nan") if arr.size == 0 else float(arr.mean())
            for arr in self._valid_per_point()
        )

    @property
    def stds(self) -> tuple[float, ...]:
        return tuple(
            float("nan") if arr.size < 2 else float(arr.std(ddof=1))
            for arr in self._valid_per_point()
        )

    @property
    def sems(self) -> tuple[float, ...]:
        return tuple(
            float("nan") if arr.size < 2 else float(arr.std(ddof=1) / np.sqrt(arr.size))
            for arr in self._valid_per_point()
        )

    @property
    def n_replicas_per_point(self) -> tuple[int, ...]:
        return tuple(arr.size for arr in self._valid_per_point())


def run_convergence_sweep(
    x_values: Sequence[float],
    analyze_fn: Callable[[float], Sequence[float | None]],
    *,
    x_label: str,
    observable_label: str,
) -> ConvergenceSweepRecord:
    """Run ``analyze_fn(x)`` for each x and aggregate into a sweep record.

    ``analyze_fn(x)`` must return a sequence of per-replica observable values
    at that x (one per replica). ``None`` entries mark missing/failed replicas.
    """
    x_tuple = tuple(float(x) for x in x_values)
    replica_tuples: list[tuple[float | None, ...]] = []
    for x in x_tuple:
        values = analyze_fn(x)
        replica_tuples.append(
            tuple(None if v is None else float(v) for v in values)
        )
    return ConvergenceSweepRecord(
        x_label=str(x_label),
        x_values=x_tuple,
        observable_label=str(observable_label),
        replica_values_per_point=tuple(replica_tuples),
    )


@dataclass(frozen=True)
class SEMvsReplicateCountRecord:
    """Empirical SEM of an observable as a function of number of replicas."""

    n_values: tuple[int, ...]
    sem_values: tuple[float, ...]
    observable_label: str

    def __post_init__(self) -> None:
        if len(self.n_values) != len(self.sem_values):
            raise ValueError(
                f"n_values (len={len(self.n_values)}) and sem_values "
                f"(len={len(self.sem_values)}) must match."
            )


def sem_vs_replicate_count(
    replica_values: Sequence[float],
    *,
    n_max: int | None = None,
    n_bootstrap: int = 500,
    random_seed: int | None = None,
    observable_label: str = "observable",
) -> SEMvsReplicateCountRecord:
    """Empirical SEM vs number of replicas via bootstrap subsampling.

    For each N in ``2..n_max``, draws ``n_bootstrap`` size-N subsamples (with
    replacement) from ``replica_values``, computes each subsample's mean, and
    reports the std of those means as the empirical SEM at N.

    In the infinite-bootstrap limit this recovers ``σ_pop / √N`` when the pool
    is representative of the population.
    """
    values = np.asarray(replica_values, dtype=float)
    if values.size < 2:
        raise ValueError("replica_values must contain at least 2 values.")
    max_n = int(n_max) if n_max is not None else int(values.size)
    if max_n < 2 or max_n > values.size:
        raise ValueError(
            f"n_max must be in [2, {values.size}]; got {max_n}."
        )
    if n_bootstrap < 2:
        raise ValueError(f"n_bootstrap must be >= 2; got {n_bootstrap}.")

    rng = np.random.default_rng(random_seed)
    ns = list(range(2, max_n + 1))
    sems: list[float] = []
    for n in ns:
        samples = rng.choice(values, size=(n_bootstrap, n), replace=True)
        boot_means = samples.mean(axis=1)
        sems.append(float(boot_means.std(ddof=0)))
    return SEMvsReplicateCountRecord(
        n_values=tuple(ns),
        sem_values=tuple(sems),
        observable_label=str(observable_label),
    )


def _import_mpl():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "matplotlib is required for convergence-sweep plots."
        ) from exc
    return plt


def plot_convergence_sweep(
    record: ConvergenceSweepRecord,
    *,
    ax=None,
    **plot_kwargs: Any,
):
    """Plot sweep means with SEM error bars vs ``record.x_values``."""
    plt = _import_mpl()
    if ax is None:
        _, ax = plt.subplots()
    means = np.array(record.means)
    sems = np.array(record.sems)
    plot_kwargs.setdefault("marker", "o")
    plot_kwargs.setdefault("capsize", 3)
    ax.errorbar(record.x_values, means, yerr=sems, **plot_kwargs)
    ax.set_xlabel(record.x_label)
    ax.set_ylabel(record.observable_label)
    ax.set_title(f"Convergence: {record.observable_label} vs {record.x_label}")
    return ax


def plot_sem_vs_replicate_count(
    record: SEMvsReplicateCountRecord,
    *,
    ax=None,
    show_expected_scaling: bool = True,
    **plot_kwargs: Any,
):
    """Plot empirical SEM vs N, with optional ``c/√N`` reference overlay."""
    plt = _import_mpl()
    if ax is None:
        _, ax = plt.subplots()
    plot_kwargs.setdefault("marker", "o")
    ax.plot(record.n_values, record.sem_values, label="empirical", **plot_kwargs)
    if show_expected_scaling and record.n_values:
        n_arr = np.asarray(record.n_values, dtype=float)
        # Anchor the c/√N curve at the first empirical point.
        anchor_c = float(record.sem_values[0]) * np.sqrt(float(record.n_values[0]))
        ax.plot(
            n_arr,
            anchor_c / np.sqrt(n_arr),
            linestyle="--",
            color="0.5",
            label=r"$c / \sqrt{N}$",
        )
        ax.legend(loc="best")
    ax.set_xlabel("Number of replicas N")
    ax.set_ylabel(f"SEM of {record.observable_label}")
    ax.set_title("SEM vs replicate count")
    return ax


__all__ = [
    "ConvergenceSweepRecord",
    "SEMvsReplicateCountRecord",
    "plot_convergence_sweep",
    "plot_sem_vs_replicate_count",
    "run_convergence_sweep",
    "sem_vs_replicate_count",
]
