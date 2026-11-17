"""Pairwise nonbond-energy traces for EVB trajectories.

Evaluates analytical Coulomb + Lennard-Jones between a user-supplied list of
atom pairs, frame by frame. A single-state monitor (one per-state parameter
set) is the common case — typical uses are:

* Tracking specific Q–Q or Q–solvent interactions across a mapping-FEP window.
* Decomposing a residue's electrostatic vs van-der-Waals contribution to a
  particular bond-breaking event.
* Building control traces next to :mod:`evbopenmm.analysis.group_contributions`
  (which is residue-wise and state-difference-oriented; this module is
  pair-wise and single-state).

Two entry points:

* :func:`compute_nonbond_pair_traces_raw` — pure numpy. Takes positions +
  ``(charges, sigmas, epsilons)`` arrays directly; no OpenMM or EVBSystem
  dependency. Use when you already have parameters in hand.
* :func:`compute_nonbond_pair_traces` — convenience wrapper that reads the
  per-state ``NonbondedForce`` parameters from an ``EVBSystem`` for a named
  state and delegates to the raw version.

Units: positions in nm, charges in e, sigma in nm, epsilon in kJ/mol.
Returns kcal/mol by default (``output_energy_unit="kcal/mol"``) for
consistency with the rest of ``evbopenmm.analysis``; pass ``"kJ/mol"`` for
kJ/mol output.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

# Match group_contributions.py so both modules agree on the Coulomb constant.
COULOMB_CONSTANT_KJ_NM_PER_E2 = 138.935456
KJ_PER_KCAL = 4.184


def _normalize_positions(positions_nm: np.ndarray) -> np.ndarray:
    a = np.asarray(positions_nm, dtype=float)
    if a.ndim == 2:
        a = a[np.newaxis, :, :]
    elif a.ndim != 3:
        raise ValueError(
            "positions_nm must have shape (n_atoms, 3) or (n_frames, n_atoms, 3); "
            f"got {a.shape}."
        )
    return a


def _apply_mixing(sigma_i: np.ndarray, sigma_j: np.ndarray,
                  epsilon_i: np.ndarray, epsilon_j: np.ndarray,
                  rule: str) -> tuple[np.ndarray, np.ndarray]:
    if rule == "lorentz-berthelot":
        return 0.5 * (sigma_i + sigma_j), np.sqrt(epsilon_i * epsilon_j)
    if rule == "geometric":
        return np.sqrt(sigma_i * sigma_j), np.sqrt(epsilon_i * epsilon_j)
    raise ValueError(f"Unknown lj_combining '{rule}'.")


def compute_nonbond_pair_traces_raw(
    positions_nm: np.ndarray,
    atom_pairs: Sequence[tuple[int, int]],
    charges_e: np.ndarray,
    sigmas_nm: np.ndarray,
    epsilons_kj_mol: np.ndarray,
    *,
    lj_combining: str = "lorentz-berthelot",
    decompose: bool = False,
    output_energy_unit: str = "kcal/mol",
) -> dict[str, np.ndarray]:
    """Per-frame Coulomb + LJ energy for each ``(i, j)`` atom pair.

    Returns a dict with key ``"total"`` of shape ``(n_frames, n_pairs)``.
    With ``decompose=True`` also returns ``"coulomb"`` and ``"lj"`` components.
    """
    pos = _normalize_positions(positions_nm)
    n_frames, n_atoms, _ = pos.shape

    q = np.asarray(charges_e, dtype=float)
    sig = np.asarray(sigmas_nm, dtype=float)
    eps = np.asarray(epsilons_kj_mol, dtype=float)
    for name, arr in (("charges_e", q), ("sigmas_nm", sig), ("epsilons_kj_mol", eps)):
        if arr.shape != (n_atoms,):
            raise ValueError(
                f"'{name}' must be shape ({n_atoms},) matching positions; got {arr.shape}."
            )

    pairs = np.asarray(atom_pairs, dtype=int)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError(f"atom_pairs must be (n_pairs, 2); got {pairs.shape}.")
    if pairs.size == 0:
        raise ValueError("atom_pairs must be non-empty.")
    lo = int(pairs.min())
    hi = int(pairs.max())
    if lo < 0 or hi >= n_atoms:
        raise ValueError(
            f"atom_pairs contains index out of range [0, {n_atoms - 1}]; got [{lo}, {hi}]."
        )
    if np.any(pairs[:, 0] == pairs[:, 1]):
        raise ValueError("atom_pairs must not contain self-pairs (i == j).")

    i = pairs[:, 0]
    j = pairs[:, 1]
    q_pair = q[i] * q[j]                                 # (n_pairs,)
    sig_pair, eps_pair = _apply_mixing(sig[i], sig[j], eps[i], eps[j], lj_combining)

    diff = pos[:, i, :] - pos[:, j, :]                   # (n_frames, n_pairs, 3)
    r_nm = np.sqrt((diff ** 2).sum(axis=-1))             # (n_frames, n_pairs)
    if np.any(r_nm == 0.0):
        raise ValueError(
            "Zero interatomic distance detected between a pair; check atom_pairs."
        )

    v_elec_kj = COULOMB_CONSTANT_KJ_NM_PER_E2 * q_pair[None, :] / r_nm
    sr6 = (sig_pair[None, :] / r_nm) ** 6
    v_vdw_kj = 4.0 * eps_pair[None, :] * (sr6 ** 2 - sr6)

    if output_energy_unit == "kcal/mol":
        scale = 1.0 / KJ_PER_KCAL
    elif output_energy_unit == "kJ/mol":
        scale = 1.0
    else:
        raise ValueError(
            f"output_energy_unit must be 'kcal/mol' or 'kJ/mol'; got {output_energy_unit!r}."
        )

    coulomb = v_elec_kj * scale
    lj = v_vdw_kj * scale
    total = coulomb + lj

    out: dict[str, np.ndarray] = {"total": total}
    if decompose:
        out["coulomb"] = coulomb
        out["lj"] = lj
    return out


def compute_nonbond_pair_traces(
    evb_system: Any,
    positions_nm: np.ndarray,
    atom_pairs: Sequence[tuple[int, int]],
    *,
    state_name: str,
    lj_combining: str = "lorentz-berthelot",
    decompose: bool = False,
    output_energy_unit: str = "kcal/mol",
) -> dict[str, np.ndarray]:
    """Pair-energy traces using per-state ``NonbondedForce`` parameters.

    Thin wrapper: pulls ``(charges, sigmas, epsilons)`` from
    ``evb_system._state_systems[state_name]`` and delegates to
    :func:`compute_nonbond_pair_traces_raw`.
    """
    from .group_contributions import _collect_state_nonbonded_params

    q, sig, eps = _collect_state_nonbonded_params(evb_system, state_name)
    return compute_nonbond_pair_traces_raw(
        positions_nm,
        atom_pairs,
        q,
        sig,
        eps,
        lj_combining=lj_combining,
        decompose=decompose,
        output_energy_unit=output_energy_unit,
    )


__all__ = [
    "COULOMB_CONSTANT_KJ_NM_PER_E2",
    "compute_nonbond_pair_traces",
    "compute_nonbond_pair_traces_raw",
]
