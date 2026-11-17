"""N-state EVB parameter fitting via scipy.optimize.minimize.

Two-state fits are the common special case: pass
``coupling_pattern=[(reactant_name, product_name)]`` and one or both of
``target_dg_dagger_kcal_mol`` / ``target_dg0_kcal_mol``. For N >= 3 diabatic
states, list every coupling pair you want to vary; all others stay at zero.

The optimizer is scipy.optimize.minimize. The default method is ``Powell``
because the objective is non-differentiable (binned free-energy profile) but
smooth enough. Grid search was retired — it does not scale past 2 parameters.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np

try:  # pragma: no cover - scipy is a required dependency; guard for clearer errors
    from scipy import optimize as _scipy_optimize
except ImportError as exc:  # pragma: no cover
    _scipy_optimize = None
    _SCIPY_IMPORT_ERROR: Exception | None = exc
else:
    _SCIPY_IMPORT_ERROR = None


CouplingPair = tuple[str, str]
AnalysisResult = Mapping[str, object]


def _canonical_pair(pair: tuple[str, str]) -> CouplingPair:
    i, j = pair
    return (i, j) if i <= j else (j, i)


def _require_scipy() -> None:
    if _scipy_optimize is None:
        raise ImportError(
            "scipy is required for fit_evb_parameters. "
            f"Underlying import error: {_SCIPY_IMPORT_ERROR!r}"
        )


def _validate_coupling_pattern(
    coupling_pattern: Sequence[tuple[str, str]],
    state_names: Sequence[str],
) -> tuple[CouplingPair, ...]:
    state_set = set(state_names)
    seen: set[CouplingPair] = set()
    normalized: list[CouplingPair] = []
    for pair in coupling_pattern:
        if len(pair) != 2:
            raise ValueError(f"Coupling pair {pair!r} must have exactly two entries.")
        i, j = pair
        if i not in state_set or j not in state_set:
            raise ValueError(f"Coupling pair ({i!r}, {j!r}) references unknown state.")
        if i == j:
            raise ValueError(f"Coupling pair ({i!r}, {j!r}) cannot be a self-coupling.")
        canonical = _canonical_pair((i, j))
        if canonical in seen:
            raise ValueError(f"Duplicate coupling pair in pattern: {canonical}.")
        seen.add(canonical)
        normalized.append(canonical)
    return tuple(normalized)


def _validate_offset_pattern(
    diabatic_offset_pattern: Sequence[str],
    state_names: Sequence[str],
) -> tuple[str, ...]:
    state_set = set(state_names)
    seen: set[str] = set()
    normalized: list[str] = []
    for name in diabatic_offset_pattern:
        if name not in state_set:
            raise ValueError(f"Offset state {name!r} is not one of the EVB state names.")
        if name in seen:
            raise ValueError(f"Duplicate state in diabatic_offset_pattern: {name!r}.")
        seen.add(name)
        normalized.append(name)
    return tuple(normalized)


def _resolve_pair_bounds(
    pair: CouplingPair,
    coupling_bounds: tuple[float, float] | Mapping[CouplingPair, tuple[float, float]],
) -> tuple[float, float]:
    if isinstance(coupling_bounds, Mapping):
        canonical = coupling_bounds.get(pair)
        if canonical is None:
            canonical = coupling_bounds.get((pair[1], pair[0]))
        if canonical is None:
            return (-20.0, 20.0)
        return (float(canonical[0]), float(canonical[1]))
    return (float(coupling_bounds[0]), float(coupling_bounds[1]))


def _resolve_offset_bounds(
    state_name: str,
    offset_bounds: tuple[float, float] | Mapping[str, tuple[float, float]],
) -> tuple[float, float]:
    if isinstance(offset_bounds, Mapping):
        canonical = offset_bounds.get(state_name)
        if canonical is None:
            return (-50.0, 50.0)
        return (float(canonical[0]), float(canonical[1]))
    return (float(offset_bounds[0]), float(offset_bounds[1]))


def fit_evb_parameters_nstate(
    *,
    state_names: Sequence[str],
    coupling_pattern: Sequence[tuple[str, str]],
    diabatic_offset_pattern: Sequence[str],
    analyze_candidate: Callable[
        [Mapping[CouplingPair, float], Mapping[str, float]], AnalysisResult
    ],
    target_dg_dagger_kcal_mol: float | None = None,
    target_dg0_kcal_mol: float | None = None,
    weight_dg_dagger: float = 1.0,
    weight_dg0: float = 1.0,
    initial_couplings_kcal_mol: Mapping[tuple[str, str], float] | None = None,
    initial_offsets_kcal_mol: Mapping[str, float] | None = None,
    coupling_bounds_kcal_mol: tuple[float, float] | Mapping[tuple[str, str], tuple[float, float]] = (-20.0, 20.0),
    offset_bounds_kcal_mol: tuple[float, float] | Mapping[str, tuple[float, float]] = (-50.0, 50.0),
    scipy_method: str = "Powell",
    scipy_tolerance: float | None = None,
    scipy_options: Mapping[str, object] | None = None,
) -> dict[str, Any]:
    """Fit a sparse set of N-state EVB coupling entries + diabatic offsets.

    The caller supplies ``analyze_candidate(couplings, offsets) -> analysis_dict``
    that runs :func:`evbopenmm.analysis.analyze_evb_profile_records` (or an
    EVBSystem method wrapper) with the candidate parameters. The objective
    is weighted squared error against ``dg_dagger`` and/or ``dg0`` targets.
    """
    _require_scipy()

    if target_dg_dagger_kcal_mol is None and target_dg0_kcal_mol is None:
        raise ValueError(
            "At least one of target_dg_dagger_kcal_mol or target_dg0_kcal_mol must be provided."
        )

    normalized_pattern = _validate_coupling_pattern(coupling_pattern, state_names)
    normalized_offsets = _validate_offset_pattern(diabatic_offset_pattern, state_names)

    if not normalized_pattern and not normalized_offsets:
        raise ValueError(
            "fit_evb_parameters_nstate needs at least one varied parameter "
            "(non-empty coupling_pattern or diabatic_offset_pattern)."
        )

    initial_couplings = {
        _canonical_pair(pair): float(value)
        for pair, value in (initial_couplings_kcal_mol or {}).items()
    }
    initial_offsets = dict(initial_offsets_kcal_mol or {})

    n_couplings = len(normalized_pattern)
    n_offsets = len(normalized_offsets)

    x0 = np.zeros(n_couplings + n_offsets, dtype=float)
    for idx, pair in enumerate(normalized_pattern):
        x0[idx] = initial_couplings.get(pair, 0.0)
    for idx, name in enumerate(normalized_offsets):
        x0[n_couplings + idx] = float(initial_offsets.get(name, 0.0))

    bounds: list[tuple[float, float]] = []
    for pair in normalized_pattern:
        bounds.append(_resolve_pair_bounds(pair, coupling_bounds_kcal_mol))
    for name in normalized_offsets:
        bounds.append(_resolve_offset_bounds(name, offset_bounds_kcal_mol))

    history: list[dict[str, Any]] = []

    def _unpack(x_vec: np.ndarray) -> tuple[dict[CouplingPair, float], dict[str, float]]:
        couplings: dict[CouplingPair, float] = {}
        for idx, pair in enumerate(normalized_pattern):
            couplings[pair] = float(x_vec[idx])
        offsets: dict[str, float] = {}
        for idx, name in enumerate(normalized_offsets):
            offsets[name] = float(x_vec[n_couplings + idx])
        return couplings, offsets

    infeasible_penalty = 1.0e12

    def _objective(x_vec: np.ndarray) -> float:
        couplings, offsets = _unpack(np.asarray(x_vec, dtype=float))
        try:
            analysis = analyze_candidate(couplings, offsets)
        except Exception as exc:  # pragma: no cover - candidate failure is rare but not impossible
            history.append(
                {
                    "couplings_kcal_mol": dict(couplings),
                    "offsets_kcal_mol": dict(offsets),
                    "dg_dagger_kcal_mol": None,
                    "dg0_kcal_mol": None,
                    "objective_kcal2_mol2": float("inf"),
                    "error": str(exc),
                }
            )
            return infeasible_penalty

        thermo = analysis["thermochemistry"]
        dg_dagger_raw = thermo.get("dg_dagger_kcal_mol")
        dg0_raw = thermo.get("dg0_kcal_mol")
        dg_dagger = None if dg_dagger_raw is None else float(dg_dagger_raw)
        dg0 = None if dg0_raw is None else float(dg0_raw)

        objective = 0.0
        if target_dg_dagger_kcal_mol is not None:
            if dg_dagger is None:
                objective += float(weight_dg_dagger) * 1.0e6
            else:
                delta = dg_dagger - float(target_dg_dagger_kcal_mol)
                objective += float(weight_dg_dagger) * delta * delta
        if target_dg0_kcal_mol is not None:
            if dg0 is None:
                objective += float(weight_dg0) * 1.0e6
            else:
                delta = dg0 - float(target_dg0_kcal_mol)
                objective += float(weight_dg0) * delta * delta

        history.append(
            {
                "couplings_kcal_mol": dict(couplings),
                "offsets_kcal_mol": dict(offsets),
                "dg_dagger_kcal_mol": dg_dagger,
                "dg0_kcal_mol": dg0,
                "objective_kcal2_mol2": float(objective),
            }
        )
        return float(objective)

    options = dict(scipy_options or {})
    minimize_kwargs: dict[str, Any] = {
        "method": scipy_method,
        "options": options,
    }
    if scipy_tolerance is not None:
        minimize_kwargs["tol"] = float(scipy_tolerance)
    if bounds:
        minimize_kwargs["bounds"] = bounds

    result = _scipy_optimize.minimize(_objective, x0, **minimize_kwargs)

    best_x = np.asarray(result.x, dtype=float).ravel()
    best_couplings, best_offsets = _unpack(best_x)
    final_analysis = analyze_candidate(best_couplings, best_offsets)
    thermo = final_analysis["thermochemistry"]
    best_dg_dagger = thermo.get("dg_dagger_kcal_mol")
    best_dg0 = thermo.get("dg0_kcal_mol")

    warnings_list: list[str] = []
    if not bool(result.success):
        warnings_list.append(f"scipy minimize did not converge: {result.message!s}")
    if target_dg_dagger_kcal_mol is not None and best_dg_dagger is None:
        warnings_list.append(
            "Best-fit profile has no detectable transition state; dG_dagger could not be measured."
        )
    if target_dg0_kcal_mol is not None and best_dg0 is None:
        warnings_list.append(
            "Best-fit profile has no identified product minimum; dG0 could not be measured."
        )

    return {
        "state_names": tuple(state_names),
        "coupling_pattern": tuple(normalized_pattern),
        "diabatic_offset_pattern": tuple(normalized_offsets),
        "best_couplings_kcal_mol": dict(best_couplings),
        "best_offsets_kcal_mol": dict(best_offsets),
        "best_dg_dagger_kcal_mol": None if best_dg_dagger is None else float(best_dg_dagger),
        "best_dg0_kcal_mol": None if best_dg0 is None else float(best_dg0),
        "objective_kcal2_mol2": float(result.fun),
        "target_dg_dagger_kcal_mol": None if target_dg_dagger_kcal_mol is None else float(target_dg_dagger_kcal_mol),
        "target_dg0_kcal_mol": None if target_dg0_kcal_mol is None else float(target_dg0_kcal_mol),
        "weight_dg_dagger": float(weight_dg_dagger),
        "weight_dg0": float(weight_dg0),
        "optimizer": {
            "method": str(scipy_method),
            "success": bool(result.success),
            "message": str(result.message),
            "nit": int(getattr(result, "nit", 0)),
            "nfev": int(getattr(result, "nfev", 0)),
        },
        "history": tuple(history),
        "best_analysis": final_analysis,
        "warnings": tuple(warnings_list),
    }


def fit_evb_parameters_replicas_nstate(
    *,
    state_names: Sequence[str],
    replica_count: int,
    coupling_pattern: Sequence[tuple[str, str]],
    diabatic_offset_pattern: Sequence[str],
    analyze_candidate: Callable[
        [Mapping[CouplingPair, float], Mapping[str, float]], AnalysisResult
    ],
    target_dg_dagger_kcal_mol: float | None = None,
    target_dg0_kcal_mol: float | None = None,
    weight_dg_dagger: float = 1.0,
    weight_dg0: float = 1.0,
    initial_couplings_kcal_mol: Mapping[tuple[str, str], float] | None = None,
    initial_offsets_kcal_mol: Mapping[str, float] | None = None,
    coupling_bounds_kcal_mol: tuple[float, float] | Mapping[tuple[str, str], tuple[float, float]] = (-20.0, 20.0),
    offset_bounds_kcal_mol: tuple[float, float] | Mapping[str, tuple[float, float]] = (-50.0, 50.0),
    scipy_method: str = "Powell",
    scipy_tolerance: float | None = None,
    scipy_options: Mapping[str, object] | None = None,
) -> dict[str, Any]:
    """Replica-ensemble variant: objective uses the aggregate dG_dagger / dG0 means.

    ``analyze_candidate(couplings, offsets)`` must return the replica-ensemble
    result produced by :func:`analyze_evb_profile_replicas_result`; the
    ``aggregate.thermochemistry.dg_dagger_mean_kcal_mol`` and
    ``aggregate.thermochemistry.dg0_mean_kcal_mol`` values drive the fit.
    """

    def _analyze_wrapping(
        couplings: Mapping[CouplingPair, float],
        offsets: Mapping[str, float],
    ) -> AnalysisResult:
        replica_result = analyze_candidate(couplings, offsets)
        aggregate = replica_result["aggregate"]
        thermo = aggregate["thermochemistry"]
        return {
            "thermochemistry": {
                "dg_dagger_kcal_mol": thermo.get("dg_dagger_mean_kcal_mol"),
                "dg0_kcal_mol": thermo.get("dg0_mean_kcal_mol"),
            },
            "replica_result": replica_result,
        }

    result = fit_evb_parameters_nstate(
        state_names=state_names,
        coupling_pattern=coupling_pattern,
        diabatic_offset_pattern=diabatic_offset_pattern,
        analyze_candidate=_analyze_wrapping,
        target_dg_dagger_kcal_mol=target_dg_dagger_kcal_mol,
        target_dg0_kcal_mol=target_dg0_kcal_mol,
        weight_dg_dagger=weight_dg_dagger,
        weight_dg0=weight_dg0,
        initial_couplings_kcal_mol=initial_couplings_kcal_mol,
        initial_offsets_kcal_mol=initial_offsets_kcal_mol,
        coupling_bounds_kcal_mol=coupling_bounds_kcal_mol,
        offset_bounds_kcal_mol=offset_bounds_kcal_mol,
        scipy_method=scipy_method,
        scipy_tolerance=scipy_tolerance,
        scipy_options=scipy_options,
    )

    best_analysis = result["best_analysis"]
    replica_result = best_analysis["replica_result"]
    result["best_replica_result"] = replica_result
    result["best_analysis"] = replica_result  # upstream consumers see the full ensemble dict
    result["replica_count"] = int(replica_count)
    return result


__all__ = [
    "fit_evb_parameters_nstate",
    "fit_evb_parameters_replicas_nstate",
]
