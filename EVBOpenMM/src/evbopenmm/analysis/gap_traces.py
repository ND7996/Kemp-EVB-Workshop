"""Per-frame derived traces from diabatic EVB state energies.

Given an already-sampled ``diabatic_energies`` array of shape
``(n_frames, n_states)`` — produced, for example, by looping
:func:`evbopenmm.analysis.sample_diabatic_state_energies_kj_mol` across the
frames of a trajectory — these helpers produce the derived time series most
EVB post-processing needs:

* :func:`diabatic_gap_trace` — the antisymmetric gap ``ΔE_ij(t) = E_j − E_i``
  between two named/indexed states. The standard EVB energy-gap reaction
  coordinate.
* :func:`evb_ground_state_trace_kj_mol` — per-frame EVB ground-state
  eigenvalue for a given N×N diabatic-coupling matrix. Loops
  ``np.linalg.eigvalsh`` per frame; N is small (2–4 for typical EVB), so the
  cost is negligible.
* :func:`mapping_potential_trace` — linear mapping potential
  ``V_m(t; λ) = Σ_i λ_i · E_i(t)`` for any λ vector (e.g. a window's λ in
  mapping FEP).

These utilities are pure-numpy and do not drive an OpenMM context. To
sample diabatic energies online during a simulation use
:func:`~evbopenmm.analysis.sample_diabatic_state_energies_kj_mol`; to apply
it along an existing trajectory, loop it while setting context positions.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def _as_2d_energies(diabatic_energies: np.ndarray) -> np.ndarray:
    """Accept ``(n_frames, n_states)`` or a single-frame ``(n_states,)``."""
    a = np.asarray(diabatic_energies, dtype=float)
    if a.ndim == 1:
        a = a[np.newaxis, :]
    elif a.ndim != 2:
        raise ValueError(
            "diabatic_energies must have shape (n_states,) or (n_frames, n_states); "
            f"got {a.shape}."
        )
    return a


def diabatic_gap_trace(
    diabatic_energies: np.ndarray,
    state_i: int,
    state_j: int,
) -> np.ndarray:
    """Per-frame energy gap ``ΔE_ij(t) = E_j(t) − E_i(t)``.

    Sign convention: positive = state ``j`` higher than ``i``. For the usual
    EVB reactant→product gap, pass ``state_i=reactant_index,
    state_j=product_index`` so the trace is negative in the reactant basin,
    zero at the crossing, positive in the product basin.

    Returns a 1-D array of shape ``(n_frames,)``.
    """
    energies = _as_2d_energies(diabatic_energies)
    n_states = energies.shape[1]
    if not (0 <= int(state_i) < n_states and 0 <= int(state_j) < n_states):
        raise ValueError(
            f"state_i={state_i}, state_j={state_j} out of range [0, {n_states - 1}]."
        )
    if int(state_i) == int(state_j):
        raise ValueError("state_i and state_j must be different.")
    return energies[:, int(state_j)] - energies[:, int(state_i)]


def evb_ground_state_trace_kj_mol(
    diabatic_energies_kj_mol: np.ndarray,
    coupling_matrix_kj_mol: np.ndarray,
) -> np.ndarray:
    """Per-frame EVB ground-state energy in kJ/mol.

    The Hamiltonian for frame ``t`` is built by placing
    ``diabatic_energies_kj_mol[t, :]`` on the diagonal and the constant
    ``coupling_matrix_kj_mol`` on the off-diagonals. Ground state is the
    lowest eigenvalue of that symmetric matrix.

    ``coupling_matrix_kj_mol`` must be ``(n_states, n_states)``, symmetric,
    and have zero diagonal (couplings only modify off-diagonal entries;
    diagonal offsets belong in the diabatic energies themselves).

    Returns a 1-D array of shape ``(n_frames,)``.
    """
    energies = _as_2d_energies(diabatic_energies_kj_mol)
    n_frames, n_states = energies.shape

    coupling = np.asarray(coupling_matrix_kj_mol, dtype=float)
    if coupling.shape != (n_states, n_states):
        raise ValueError(
            f"coupling_matrix_kj_mol must be ({n_states}, {n_states}); got {coupling.shape}."
        )
    if not np.allclose(coupling, coupling.T, atol=1e-12):
        raise ValueError("coupling_matrix_kj_mol must be symmetric.")
    if not np.allclose(np.diag(coupling), 0.0, atol=1e-12):
        raise ValueError("coupling_matrix_kj_mol must have zero diagonal.")

    out = np.empty(n_frames, dtype=float)
    hamiltonian = coupling.copy()
    diag_idx = np.diag_indices(n_states)
    for t in range(n_frames):
        hamiltonian[diag_idx] = energies[t]
        out[t] = float(np.linalg.eigvalsh(hamiltonian)[0])
    return out


def mapping_potential_trace(
    diabatic_energies: np.ndarray,
    lambda_values: Sequence[float],
) -> np.ndarray:
    """Linear mapping potential ``V_m(t; λ) = Σ_i λ_i · E_i(t)``.

    The λ vector is not normalised here — pass it as-is (for an EVB mapping
    FEP window at endpoint λ=(1, 0, …) the trace recovers ``E_0(t)``). Must
    have length matching the state count in ``diabatic_energies``.

    Returns a 1-D array of shape ``(n_frames,)``.
    """
    energies = _as_2d_energies(diabatic_energies)
    n_states = energies.shape[1]
    lam = np.asarray(lambda_values, dtype=float).ravel()
    if lam.shape != (n_states,):
        raise ValueError(
            f"lambda_values must have length {n_states}; got shape {lam.shape}."
        )
    return energies @ lam


__all__ = [
    "diabatic_gap_trace",
    "evb_ground_state_trace_kj_mol",
    "mapping_potential_trace",
]
