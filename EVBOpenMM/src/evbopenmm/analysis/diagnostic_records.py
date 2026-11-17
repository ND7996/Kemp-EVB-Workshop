"""Diagnostic records for an EVB mapping-FEP analysis.

Thin, frozen dataclasses that pick the information each diagnostic panel
needs out of an existing ``analyze_evb_profile_records`` result. Each record
is pure data — plotting lives in a separate module — so a campaign can
freeze these to JSON and a notebook can render them later without rerunning
analysis.

Four panels are covered here; together they reproduce the diagnostic grid
commonly shown in the EVB literature for per-replica health checks:

* :class:`ProfileCurve` — the free-energy profile vs energy-gap reaction
  coordinate (the headline curve).
* :class:`WindowOverlap` — per-λ-window mean energy gap, for verifying that
  the λ ladder spans the reaction coordinate smoothly.
* :class:`ForwardReverseResidual` — per-window ``ΔGf − ΔGr`` residual, for
  flagging windows with poor forward/reverse overlap.
* :class:`StateEnergyVsLambda` — per-state mean potential energy across the λ
  ladder, for flagging rare events that jerk a state energy.

Three of the four panels read directly from the analysis dict (no
recomputation). Only :class:`StateEnergyVsLambda` aggregates the raw sample
data — a single pass over ``analysis['windows']``.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Literal

from .energetics import energy_scale_to_kj, kj_to_kcal

NormalizationMode = Literal["raw", "referenced"]


def _profile_key_for_normalization(normalization: NormalizationMode) -> str:
    if normalization == "referenced":
        return "state_free_energies_referenced_kcal_mol"
    if normalization == "raw":
        return "state_free_energies_kcal_mol"
    raise ValueError(
        f"Unknown normalization '{normalization}'; use one of 'raw', 'referenced'."
    )


@dataclass(frozen=True)
class ProfileCurve:
    """A single state's free-energy profile along the energy-gap RC.

    ``state_name='evb_ground'`` is the adiabatic ground-state curve; the
    other valid keys are the diabatic-state names from the analysis.
    """

    gap_centers_kcal_mol: tuple[float, ...]
    free_energy_kcal_mol: tuple[float | None, ...]
    bin_counts: tuple[int, ...]
    state_name: str
    normalization: NormalizationMode

    @classmethod
    def from_analysis(
        cls,
        analysis: Mapping[str, Any],
        *,
        state_name: str = "evb_ground",
        normalization: NormalizationMode = "raw",
    ) -> "ProfileCurve":
        profile = analysis["profile"]
        curves = profile[_profile_key_for_normalization(normalization)]
        if state_name not in curves:
            raise ValueError(
                f"state_name '{state_name}' not found in profile; "
                f"available: {sorted(curves)}."
            )
        return cls(
            gap_centers_kcal_mol=tuple(float(x) for x in profile["gap_centers_kcal_mol"]),
            free_energy_kcal_mol=tuple(
                None if v is None else float(v) for v in curves[state_name]
            ),
            bin_counts=tuple(int(c) for c in profile["bin_counts"]),
            state_name=str(state_name),
            normalization=normalization,
        )


@dataclass(frozen=True)
class WindowOverlap:
    """Per-λ-window mean energy gap — the "binning overlap" diagnostic.

    Plotted against ``lambda_values`` it should sweep smoothly from the
    reactant-basin regime (strongly negative gap) to the product-basin
    regime (strongly positive gap) without gaps or back-tracking. Jumps or
    non-monotonic regions flag a poorly sampled window.
    """

    lambda_values: tuple[float, ...]
    mean_gap_kcal_mol: tuple[float, ...]
    reference_state_name: str

    @classmethod
    def from_analysis(
        cls,
        analysis: Mapping[str, Any],
        *,
        reference_state_index: int = 0,
    ) -> "WindowOverlap":
        state_names = analysis["state_names"]
        if not 0 <= int(reference_state_index) < len(state_names):
            raise ValueError(
                f"reference_state_index={reference_state_index} out of range "
                f"[0, {len(state_names) - 1}]."
            )
        ref_name = state_names[int(reference_state_index)]
        lambdas: list[float] = []
        gaps: list[float] = []
        for summary in analysis["window_summaries"]:
            lambdas.append(float(summary["lambda_values"][ref_name]))
            gaps.append(float(summary["mean_reaction_coordinate_kcal_mol"]))
        return cls(
            lambda_values=tuple(lambdas),
            mean_gap_kcal_mol=tuple(gaps),
            reference_state_name=str(ref_name),
        )


@dataclass(frozen=True)
class ForwardReverseResidual:
    """Per-window ``ΔGf − ΔGr`` — the forward/reverse FEP residual.

    For a converged window the forward and reverse edge free-energy
    increments are equal (up to kT noise), so their difference should sit
    flat near zero. A spike in the residual at any window marks poor phase
    overlap — often a rare event or a gap in the λ ladder.
    """

    window_index: tuple[int, ...]
    forward_edge_dg_kcal_mol: tuple[float, ...]
    reverse_edge_dg_kcal_mol: tuple[float, ...]
    residual_kcal_mol: tuple[float, ...]

    @classmethod
    def from_analysis(cls, analysis: Mapping[str, Any]) -> "ForwardReverseResidual":
        fep = analysis["fep"]
        fwd = tuple(float(v) for v in fep["forward_edge_dg_kcal_mol"])
        rev = tuple(float(v) for v in fep["reverse_edge_dg_kcal_mol"])
        if len(fwd) != len(rev):
            raise ValueError(
                f"forward_edge_dg_kcal_mol (len={len(fwd)}) and "
                f"reverse_edge_dg_kcal_mol (len={len(rev)}) must have equal length."
            )
        return cls(
            window_index=tuple(range(len(fwd))),
            forward_edge_dg_kcal_mol=fwd,
            reverse_edge_dg_kcal_mol=rev,
            residual_kcal_mol=tuple(f - r for f, r in zip(fwd, rev)),
        )


@dataclass(frozen=True)
class StateEnergyVsLambda:
    """Per-state mean potential energy across the λ ladder.

    The fundamental EVB-health plot: each diabatic state's mean energy
    should march monotonically along λ without discontinuities. Kinks
    indicate rare events (e.g. a flipping water or a residue that
    reorganised mid-trajectory in only a subset of replicas).

    The ``reference_state_index`` selects which state's λ component defines
    the x-axis; in a 2-state run passing 0 recovers the usual ``λ``.
    """

    lambda_values: tuple[float, ...]
    state_names: tuple[str, ...]
    mean_energies_kcal_mol: Mapping[str, tuple[float, ...]]
    reference_state_name: str

    @classmethod
    def from_analysis(
        cls,
        analysis: Mapping[str, Any],
        *,
        input_energy_unit: str = "kcal/mol",
        reference_state_index: int = 0,
    ) -> "StateEnergyVsLambda":
        state_names = tuple(analysis["state_names"])
        n_states = len(state_names)
        if not 0 <= int(reference_state_index) < n_states:
            raise ValueError(
                f"reference_state_index={reference_state_index} out of range "
                f"[0, {n_states - 1}]."
            )
        ref_name = state_names[int(reference_state_index)]
        skip = int(analysis.get("skip_samples", 0))
        sample_scale_to_kj = energy_scale_to_kj(input_energy_unit)

        lambdas: list[float] = []
        per_state: dict[str, list[float]] = {name: [] for name in state_names}

        windows = analysis["windows"]
        summaries = analysis["window_summaries"]
        if len(windows) != len(summaries):
            raise ValueError(
                f"'windows' (len={len(windows)}) and 'window_summaries' "
                f"(len={len(summaries)}) length mismatch."
            )

        for window, summary in zip(windows, summaries):
            samples = tuple(window.samples[skip:])
            if not samples:
                continue
            lambdas.append(float(summary["lambda_values"][ref_name]))
            sums = [0.0] * n_states
            for sample in samples:
                raw = sample.state_energies_kj_mol
                if len(raw) != n_states:
                    raise ValueError(
                        f"Sample has {len(raw)} state energies; expected {n_states}."
                    )
                for k in range(n_states):
                    sums[k] += float(raw[k])
            n = float(len(samples))
            for k, name in enumerate(state_names):
                mean_kj = sample_scale_to_kj * sums[k] / n
                per_state[name].append(float(kj_to_kcal(mean_kj)))
        return cls(
            lambda_values=tuple(lambdas),
            state_names=state_names,
            mean_energies_kcal_mol={
                name: tuple(values) for name, values in per_state.items()
            },
            reference_state_name=str(ref_name),
        )


__all__ = [
    "ForwardReverseResidual",
    "NormalizationMode",
    "ProfileCurve",
    "StateEnergyVsLambda",
    "WindowOverlap",
]
