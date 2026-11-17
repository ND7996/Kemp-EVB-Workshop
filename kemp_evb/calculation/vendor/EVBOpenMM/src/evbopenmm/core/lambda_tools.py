from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from numbers import Real


def resolve_lambda_parameters(
    current_mapping: Mapping[str, str],
    *,
    state_names: Sequence[str],
    overrides: Mapping[str, str] | None,
) -> dict[str, str]:
    mapping = dict(current_mapping)
    if overrides:
        mapping.update(overrides)
    for state in state_names:
        mapping.setdefault(state, f"lambda_{state}")
    return mapping


def coerce_lambda_sequence(
    lambda_values: Iterable[float] | float,
    *,
    state_names: Sequence[str],
) -> tuple[float, ...]:
    n_states = len(state_names)
    if n_states < 2:
        raise RuntimeError("Lambda assignment requires at least two EVB states.")

    if isinstance(lambda_values, Real):
        provided = (float(lambda_values),)
    else:
        provided = tuple(float(value) for value in lambda_values)
    if len(provided) != n_states - 1:
        raise ValueError(
            f"Expected {n_states - 1} lambda values for states {tuple(state_names[:-1])}, but received {len(provided)}."
        )

    tolerance = 1e-8
    for idx, value in enumerate(provided):
        if value < -tolerance:
            raise ValueError(
                f"Lambda for state '{state_names[idx]}' is negative ({value}). Provide non-negative weights."
            )

    provided_sum = sum(provided)
    if provided_sum > 1.0 + tolerance:
        raise ValueError(
            f"Sum of provided lambdas ({provided_sum}) exceeds 1. Supply values whose sum is <= 1."
        )

    residual = 1.0 - provided_sum
    if residual < -tolerance:
        raise ValueError(
            f"Computed lambda for state '{state_names[-1]}' is negative ({residual}). Adjust the input values."
        )
    if residual < 0.0:
        residual = 0.0

    return tuple(provided) + (residual,)


__all__ = [
    "resolve_lambda_parameters",
    "coerce_lambda_sequence",
]
