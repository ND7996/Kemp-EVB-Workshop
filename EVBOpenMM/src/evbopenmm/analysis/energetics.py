from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
import math

from openmm import app, unit

from ..constants import KJ_PER_KCAL, REPORTING_ENERGY_UNIT


def sample_diabatic_state_energies_kj_mol(
    *,
    simulation: app.Simulation,
    state_names: Sequence[str],
    context_lambda_snapshot: Callable[[object], Mapping[str, float]],
    apply_indicator_lambdas: Callable[..., None],
    apply_context_lambdas: Callable[[object, Mapping[str, float]], None],
) -> tuple[float, ...]:
    context = simulation.context
    original_lambdas = context_lambda_snapshot(context)
    energies: list[float] = []
    try:
        for state_name in state_names:
            apply_indicator_lambdas(context, active_state=state_name)
            state = context.getState(getEnergy=True)
            energies.append(float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)))
    finally:
        apply_context_lambdas(context, original_lambdas)
    return tuple(energies)


def sample_diabatic_state_energies_via_force_groups_kj_mol(
    *,
    simulation: app.Simulation,
    state_names: Sequence[str],
    state_group_masks: Mapping[str, int],
    shared_group_mask: int,
    context_lambda_snapshot: Callable[[object], Mapping[str, float]],
    apply_context_lambdas: Callable[[object, Mapping[str, float]], None],
) -> tuple[float, ...]:
    """Same return shape as ``sample_diabatic_state_energies_kj_mol`` but
    avoids the per-state setParameter/getEnergy round-trip.

    Sets every ``lambda_<state>`` to 1.0 once so the per-state force groups
    read their unscaled contribution, then issues one
    ``getState(getEnergy=True, groups=mask)`` per diabatic state with
    ``mask = shared_group_mask | state_group_masks[state]``. Bit-exactly
    equivalent to the legacy path on every system tested; ~2× faster on
    CUDA because each state read no longer triggers a full re-evaluation
    of every force.
    """
    context = simulation.context
    original_lambdas = context_lambda_snapshot(context)
    energies: list[float] = []
    try:
        apply_context_lambdas(context, {state: 1.0 for state in state_names})
        for state_name in state_names:
            mask = int(shared_group_mask) | int(state_group_masks[state_name])
            state = context.getState(getEnergy=True, groups=mask)
            energies.append(float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)))
    finally:
        apply_context_lambdas(context, original_lambdas)
    return tuple(energies)


def effective_energy_from_state_energies(
    lambda_values: Sequence[float],
    state_energies_kj_mol: Sequence[float],
) -> float:
    return sum(weight * energy for weight, energy in zip(lambda_values, state_energies_kj_mol))


def reaction_coordinate_from_state_energies(
    coefficients: Sequence[float],
    state_energies_kj_mol: Sequence[float],
) -> float:
    return sum(weight * energy for weight, energy in zip(coefficients, state_energies_kj_mol))


def log_mean_exp(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("Cannot evaluate log-mean-exp for an empty sequence.")
    max_value = max(values)
    shifted_sum = sum(math.exp(value - max_value) for value in values)
    return max_value + math.log(shifted_sum / float(len(values)))


def kj_to_kcal(value_kj_mol: float) -> float:
    return float(value_kj_mol) / KJ_PER_KCAL


def energy_scale_to_kj(input_energy_unit: str) -> float:
    normalized = str(input_energy_unit).strip().lower()
    if normalized in {"kj", "kj/mol", "kilojoule", "kilojoules", "kilojoule/mol", "kilojoules/mol"}:
        return 1.0
    if normalized in {
        "kcal",
        "kcal/mol",
        "kilocalorie",
        "kilocalories",
        "kilocalorie/mol",
        "kilocalories/mol",
    }:
        return KJ_PER_KCAL
    raise ValueError(
        f"Unsupported 'input_energy_unit'. Use either '{REPORTING_ENERGY_UNIT}' or 'kJ/mol'."
    )


def evb_ground_state_energy_kj_mol(
    state1_energy_kj_mol: float,
    state2_energy_kj_mol: float,
    *,
    alpha_kj_mol: float,
    hij_kj_mol: float,
) -> float:
    shifted_gap = state1_energy_kj_mol - state2_energy_kj_mol - alpha_kj_mol
    coupling_term = math.sqrt(shifted_gap * shifted_gap + 4.0 * hij_kj_mol * hij_kj_mol)
    return 0.5 * (state1_energy_kj_mol + state2_energy_kj_mol + alpha_kj_mol - coupling_term)


def evb_ground_state_energy_nstate_kj_mol(
    diabatic_energies_kj_mol: Sequence[float],
    coupling_matrix_kj_mol: Sequence[Sequence[float]],
) -> float:
    """Ground-state EVB energy as the lowest eigenvalue of the N×N Hamiltonian.

    Builds ``H = diag(diabatic_energies) + coupling_matrix`` and returns its
    smallest eigenvalue. The caller is responsible for:

    - folding any diagonal ``alpha`` shifts (relative-state offsets) into the
      ``diabatic_energies`` vector before calling;
    - providing a **symmetric** ``coupling_matrix`` with zeros on the diagonal.

    For N=2 this is mathematically identical to :func:`evb_ground_state_energy_kj_mol`
    (which is retained because a closed-form 2×2 root avoids the numpy dispatch
    cost in hot loops). The matrix form exists as the scaffold for Phase D
    multi-state analysis work — it is not yet wired into the public
    analyze/fit pipelines.
    """

    import numpy as np

    diag = [float(value) for value in diabatic_energies_kj_mol]
    n_states = len(diag)
    if n_states < 2:
        raise ValueError("N-state ground-state energy requires at least 2 diabatic states.")

    rows = [list(map(float, row)) for row in coupling_matrix_kj_mol]
    if len(rows) != n_states or any(len(row) != n_states for row in rows):
        raise ValueError(
            f"'coupling_matrix_kj_mol' must be a {n_states}x{n_states} matrix; "
            f"received shape {len(rows)}x{(len(rows[0]) if rows else 0)}."
        )

    tolerance = 1e-8
    for i in range(n_states):
        if abs(rows[i][i]) > tolerance:
            raise ValueError(
                f"'coupling_matrix_kj_mol' diagonal entry [{i}][{i}] must be zero "
                "(fold diabatic-state offsets into 'diabatic_energies_kj_mol' instead)."
            )
        for j in range(i + 1, n_states):
            if abs(rows[i][j] - rows[j][i]) > tolerance:
                raise ValueError(
                    f"'coupling_matrix_kj_mol' must be symmetric; "
                    f"[{i}][{j}]={rows[i][j]} does not match [{j}][{i}]={rows[j][i]}."
                )

    hamiltonian = np.array(rows, dtype=float)
    for i in range(n_states):
        hamiltonian[i, i] = diag[i]
    eigenvalues = np.linalg.eigvalsh(hamiltonian)
    return float(eigenvalues[0])


def evb_ground_state_energy_gaussian_coupling_nstate_kj_mol(
    diabatic_energies_kj_mol: Sequence[float],
    coupling_amplitude_matrix_kj_mol: Sequence[Sequence[float]],
    coupling_width_matrix_per_kj2_mol2: Sequence[Sequence[float]],
) -> float:
    """Ground-state EVB energy with Gaussian-in-energy-gap off-diagonal coupling.

    The off-diagonal coupling between diabatic states i and j is evaluated per
    sample as ``H_ij(dE) = A_ij * exp(-mu_ij * (E_i - E_j)^2)``. This is the
    classic Warshel/Aqvist form: a small amplitude ``A`` localized to the
    crossing region by a width ``mu``, so the coupling does not artificially
    depress the reactant- and product-state minima.

    Both matrices must be symmetric, N×N, with zero on the diagonal.

    See :func:`evb_ground_state_energy_nstate_kj_mol` for the constant-coupling
    variant.
    """

    import numpy as np

    diag = [float(value) for value in diabatic_energies_kj_mol]
    n_states = len(diag)
    if n_states < 2:
        raise ValueError("N-state ground-state energy requires at least 2 diabatic states.")

    amplitude_rows = [list(map(float, row)) for row in coupling_amplitude_matrix_kj_mol]
    width_rows = [list(map(float, row)) for row in coupling_width_matrix_per_kj2_mol2]
    for label, rows in (
        ("coupling_amplitude_matrix_kj_mol", amplitude_rows),
        ("coupling_width_matrix_per_kj2_mol2", width_rows),
    ):
        if len(rows) != n_states or any(len(row) != n_states for row in rows):
            raise ValueError(f"'{label}' must be a {n_states}x{n_states} matrix.")

    tolerance = 1e-8
    for i in range(n_states):
        for label, rows in (
            ("coupling_amplitude_matrix_kj_mol", amplitude_rows),
            ("coupling_width_matrix_per_kj2_mol2", width_rows),
        ):
            if abs(rows[i][i]) > tolerance:
                raise ValueError(
                    f"'{label}' diagonal entry [{i}][{i}] must be zero."
                )
            for j in range(i + 1, n_states):
                if abs(rows[i][j] - rows[j][i]) > tolerance:
                    raise ValueError(f"'{label}' must be symmetric at [{i}][{j}].")

    hamiltonian = np.zeros((n_states, n_states), dtype=float)
    for i in range(n_states):
        hamiltonian[i, i] = diag[i]
        for j in range(i + 1, n_states):
            amplitude = amplitude_rows[i][j]
            width = width_rows[i][j]
            if width < 0.0:
                raise ValueError(
                    f"'coupling_width_matrix_per_kj2_mol2' entry [{i}][{j}]={width} "
                    "must be >= 0 (nonnegative Gaussian width)."
                )
            delta_e = diag[i] - diag[j]
            coupling = amplitude * math.exp(-width * delta_e * delta_e)
            hamiltonian[i, j] = coupling
            hamiltonian[j, i] = coupling
    eigenvalues = np.linalg.eigvalsh(hamiltonian)
    return float(eigenvalues[0])


def evb_distance_coupling_kj_mol(
    distance_nm: float,
    *,
    amplitude_kj_mol: float,
    decay_per_nm: float = 0.0,
    reference_distance_nm: float = 0.0,
) -> float:
    """Off-diagonal coupling as an exponential in a geometric distance.

    ``H12(r) = A * exp(-mu * (r - r0))`` — the standard Warshel/Q6 distance-dependent
    EVB coupling (Q6 ``qfep``: ``Hij = A*exp(-mu*(rkl - rxy0))``). ``r`` is normally
    the distance between the two reactive atoms that define the bond being made/broken.

    With ``decay_per_nm == 0`` this reduces to the constant coupling ``A``, so the
    constant-H12 workflows are the ``mu = 0`` special case.
    """
    if decay_per_nm == 0.0:
        return float(amplitude_kj_mol)
    return float(amplitude_kj_mol) * math.exp(
        -float(decay_per_nm) * (float(distance_nm) - float(reference_distance_nm))
    )


def evb_ground_state_energy_distance_coupling_kj_mol(
    state1_energy_kj_mol: float,
    state2_energy_kj_mol: float,
    *,
    alpha_kj_mol: float,
    coupling_distance_nm: float,
    amplitude_kj_mol: float,
    decay_per_nm: float = 0.0,
    reference_distance_nm: float = 0.0,
) -> float:
    """2-state EVB ground state with a distance-dependent off-diagonal coupling.

    Computes ``H12 = A*exp(-mu*(r - r0))`` from the supplied reactive-pair distance
    and feeds it into the closed-form 2x2 ground state
    (:func:`evb_ground_state_energy_kj_mol`). With ``decay_per_nm == 0`` it is exactly
    the constant-H12 result, so this is a strict superset of the existing path.
    """
    hij = evb_distance_coupling_kj_mol(
        coupling_distance_nm,
        amplitude_kj_mol=amplitude_kj_mol,
        decay_per_nm=decay_per_nm,
        reference_distance_nm=reference_distance_nm,
    )
    return evb_ground_state_energy_kj_mol(
        state1_energy_kj_mol,
        state2_energy_kj_mol,
        alpha_kj_mol=alpha_kj_mol,
        hij_kj_mol=hij,
    )


def profile_min_normalize(values: Sequence[float | None]) -> tuple[float | None, ...]:
    finite = [value for value in values if value is not None]
    if not finite:
        return tuple(values)
    offset = min(finite)
    normalized: list[float | None] = []
    for value in values:
        if value is None:
            normalized.append(None)
        else:
            normalized.append(value - offset)
    return tuple(normalized)


def mean_std_sem(values: Sequence[float]) -> tuple[float, float, float]:
    if not values:
        raise ValueError("Cannot compute statistics for an empty sequence.")
    n_values = len(values)
    mean = sum(values) / float(n_values)
    if n_values == 1:
        return mean, 0.0, 0.0
    variance = sum((value - mean) ** 2 for value in values) / float(n_values - 1)
    std = math.sqrt(max(variance, 0.0))
    sem = std / math.sqrt(float(n_values))
    return mean, std, sem


def inclusive_scan_values(
    bounds: tuple[float, float],
    step: float,
    *,
    name: str,
) -> tuple[float, ...]:
    if len(bounds) != 2:
        raise ValueError(f"'{name}_bounds' must be a (lower, upper) pair.")
    lower = float(bounds[0])
    upper = float(bounds[1])
    step_value = float(step)
    if step_value <= 0.0:
        raise ValueError(f"'{name}_step' must be > 0.")
    if lower > upper:
        raise ValueError(f"'{name}_bounds' must satisfy lower <= upper.")
    tolerance = max(1e-12, abs(step_value) * 1e-9)
    values: list[float] = []
    cursor = lower
    while cursor <= upper + tolerance:
        values.append(cursor)
        cursor += step_value
    if values:
        if values[-1] > upper + tolerance:
            values.pop()
        elif abs(values[-1] - upper) > tolerance:
            values.append(upper)
    return tuple(values)


__all__ = [
    "sample_diabatic_state_energies_kj_mol",
    "sample_diabatic_state_energies_via_force_groups_kj_mol",
    "effective_energy_from_state_energies",
    "reaction_coordinate_from_state_energies",
    "log_mean_exp",
    "kj_to_kcal",
    "energy_scale_to_kj",
    "evb_ground_state_energy_kj_mol",
    "evb_ground_state_energy_nstate_kj_mol",
    "evb_ground_state_energy_gaussian_coupling_nstate_kj_mol",
    "profile_min_normalize",
    "mean_std_sem",
    "inclusive_scan_values",
]
