"""Linear-response-approximation (LRA) group-contribution analysis for EVB.

Decomposes ⟨ΔE⟩ = ⟨E_state2 − E_state1⟩ into per-residue contributions on a
trajectory sampled at a single λ endpoint. Combining the reactant-basin
(λ≈1) and product-basin (λ≈0) ensembles gives the LRA:

    LRA_r = ½·(⟨ΔE_r⟩_reactant + ⟨ΔE_r⟩_product)
    REORG_r = ⟨ΔE_r⟩_reactant − LRA_r

Ionized residues (ARG, LYS, HIP, GLU, ASP) may be scaled by 0.5 per
Warshel et al., J. Phys. Chem. B 101, 5583 (1997), doi:10.1021/jp962478o —
a heuristic correction for overestimated electrostatic contributions at
charged side chains in a linear-response framework. Scaling is opt-in
(``apply_warshel_scaling=True``), not applied by default.

Implementation: analytical pairwise Coulomb + Lennard-Jones using
state-specific per-atom charges, σ, ε read from each per-state
``NonbondedForce``. Only residue–solute pairs contribute to ΔE; residue–
residue terms vanish because non-solute atom parameters are state-independent.

This module is OpenMM-only (no mdtraj dependency). Inputs are raw position
arrays (shape ``(n_frames, n_atoms, 3)``) in nanometers plus an OpenMM
``Topology`` to enumerate residues.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

# Coulomb constant in OpenMM units: k · e²·nm·(kJ/mol)⁻¹
#   138.9354... kJ/mol × nm / e²   (this is OpenMM's ONE_4PI_EPS0)
COULOMB_CONSTANT_KJ_NM_PER_E2 = 138.935456

KJ_PER_KCAL = 4.184
NM_PER_ANGSTROM = 0.1


def _reaction_field_coeffs(cutoff_nm: float, dielectric: float) -> tuple[float, float]:
    """Barker-Watts reaction-field coefficients (k_rf, c_rf) for cutoff ``cutoff_nm``.

    The shifted Coulomb kernel ``1/r + k_rf r^2 - c_rf`` goes smoothly to zero at the
    cutoff, so the environment sum CONVERGES even around a net-charged solute (a bare
    ``1/r`` truncated sum is only conditionally convergent and oscillates with cutoff).
    ``dielectric`` is the reaction-field continuum permittivity; ``math.inf`` gives the
    conductor limit (k_rf = 1/2rc^3, c_rf = 3/2rc)."""
    if cutoff_nm is None or cutoff_nm <= 0:
        raise ValueError("reaction_field=True requires a positive cutoff_nm.")
    rc = float(cutoff_nm)
    if dielectric == math.inf:
        return 1.0 / (2.0 * rc**3), 3.0 / (2.0 * rc)
    e = float(dielectric)
    if e < 1.0:
        raise ValueError(f"rf_dielectric must be >= 1 (or math.inf); got {dielectric!r}.")
    k_rf = (1.0 / rc**3) * (e - 1.0) / (2.0 * e + 1.0)
    c_rf = (1.0 / rc) * (3.0 * e) / (2.0 * e + 1.0)
    return k_rf, c_rf

# Warshel 1997 ionized-residue scaling factor (half-charge correction).
WARSHEL_SCALING = 0.5
IONIZED_RESIDUE_NAMES = frozenset({"ARG", "LYS", "HIP", "GLU", "ASP"})


@dataclass
class ResidueContribution:
    """Per-residue ΔE_r decomposition across a trajectory.

    Fields use kcal/mol units. ``delta_e_*`` arrays have shape ``(n_frames,)``.
    """

    residue_index: int
    residue_name: str
    chain_id: str
    atom_indices: tuple[int, ...]
    delta_e_elec_kcal_mol: np.ndarray
    delta_e_vdw_kcal_mol: np.ndarray
    warshel_scaled: bool = False

    @property
    def mean_elec(self) -> float:
        return float(self.delta_e_elec_kcal_mol.mean())

    @property
    def mean_vdw(self) -> float:
        return float(self.delta_e_vdw_kcal_mol.mean())

    @property
    def mean_total(self) -> float:
        return float((self.delta_e_elec_kcal_mol + self.delta_e_vdw_kcal_mol).mean())

    @property
    def std_total(self) -> float:
        return float((self.delta_e_elec_kcal_mol + self.delta_e_vdw_kcal_mol).std())

    @property
    def sem_total(self) -> float:
        n = self.delta_e_elec_kcal_mol.size
        if n < 2:
            return float("nan")
        return self.std_total / float(np.sqrt(n))


def _collect_state_nonbonded_params(
    evb_system: Any, state_name: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract (charge_e, sigma_nm, epsilon_kj_mol) per atom from the given
    state's ``NonbondedForce``. Returns three 1-D arrays of length n_atoms.
    """
    import openmm  # local import to avoid top-level dependency cycle
    from openmm import unit

    try:
        omm_sys = evb_system._state_systems[state_name]
    except (AttributeError, KeyError) as exc:
        raise ValueError(
            f"EVBSystem has no state '{state_name}'. Call build_state_systems() first."
        ) from exc

    nb = None
    for fi in range(omm_sys.getNumForces()):
        f = omm_sys.getForce(fi)
        if isinstance(f, openmm.NonbondedForce):
            nb = f
            break
    if nb is None:
        raise ValueError(f"State '{state_name}' has no NonbondedForce.")

    n = nb.getNumParticles()
    q = np.empty(n, dtype=float)
    sigma = np.empty(n, dtype=float)
    epsilon = np.empty(n, dtype=float)
    for i in range(n):
        q_i, sig_i, eps_i = nb.getParticleParameters(i)
        q[i] = float(q_i.value_in_unit(unit.elementary_charge))
        sigma[i] = float(sig_i.value_in_unit(unit.nanometer))
        epsilon[i] = float(eps_i.value_in_unit(unit.kilojoules_per_mole))
    return q, sigma, epsilon


def _residue_map_from_topology(topology: Any) -> list[tuple[int, str, str, tuple[int, ...]]]:
    """Walk an OpenMM ``Topology`` and return a list of
    ``(residue_index, residue_name, chain_id, (atom_indices...))`` tuples.

    Residue index follows ``residue.index`` which is unique across the topology.
    Chain id falls back to ``str(chain.index)`` if the chain has no id.
    """
    out: list[tuple[int, str, str, tuple[int, ...]]] = []
    for residue in topology.residues():
        atoms = tuple(atom.index for atom in residue.atoms())
        chain_id = residue.chain.id if residue.chain.id else str(residue.chain.index)
        out.append((int(residue.index), str(residue.name).upper(), str(chain_id), atoms))
    return out


def _normalize_positions(positions_nm: np.ndarray) -> np.ndarray:
    """Accept ``(n_atoms, 3)`` or ``(n_frames, n_atoms, 3)``, return 3-D."""
    a = np.asarray(positions_nm, dtype=float)
    if a.ndim == 2:
        a = a[np.newaxis, :, :]
    elif a.ndim != 3:
        raise ValueError(
            f"positions_nm must have shape (n_atoms, 3) or (n_frames, n_atoms, 3); got {a.shape}."
        )
    return a


def _apply_mixing_lb(sigma_i: np.ndarray, sigma_j: np.ndarray,
                    epsilon_i: np.ndarray, epsilon_j: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Lorentz-Berthelot: σ arithmetic, ε geometric."""
    return 0.5 * (sigma_i + sigma_j), np.sqrt(epsilon_i * epsilon_j)


def _apply_mixing_geometric(sigma_i: np.ndarray, sigma_j: np.ndarray,
                            epsilon_i: np.ndarray, epsilon_j: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """OPLS geometric: both σ and ε geometric."""
    return np.sqrt(sigma_i * sigma_j), np.sqrt(epsilon_i * epsilon_j)


def compute_per_residue_delta_e(
    evb_system: Any,
    topology: Any,
    positions_nm: np.ndarray,
    *,
    solute_atom_indices: Sequence[int],
    state_A: str = "reactant",
    state_B: str = "product",
    cutoff_nm: float | None = None,
    lj_combining: str = "lorentz-berthelot",
) -> list[ResidueContribution]:
    """Compute per-residue ``ΔE = E_B − E_A`` on a trajectory.

    Parameters
    ----------
    evb_system :
        The ``EVBSystem`` whose per-state ``NonbondedForce`` parameters define
        state A and state B.
    topology :
        OpenMM ``Topology`` for the solvated EVB system (scaffold topology).
        Used to enumerate residues.
    positions_nm :
        Positions in nanometers. Either a single frame ``(n_atoms, 3)`` or a
        trajectory ``(n_frames, n_atoms, 3)``. Must match the scaffold atom
        ordering.
    solute_atom_indices :
        Indices of Q-atoms (solute). Residue-wise contributions are computed
        as residue ↔ solute interactions; residue–residue terms drop because
        non-solute parameters are state-independent.
    state_A, state_B :
        Names of the two EVB states; ΔE = E_B − E_A.
    cutoff_nm :
        If given, pairs with r > cutoff are zeroed. Default None: no cutoff
        (appropriate for small spherical SCAAS systems with cutoff ≥ diameter).
    lj_combining :
        ``"lorentz-berthelot"`` (OpenMM default; σ arithmetic, ε geometric)
        or ``"geometric"`` (OPLS style; both σ and ε geometric).

    Returns
    -------
    list[ResidueContribution]
        One entry per residue in the topology (excluding residues entirely
        contained in ``solute_atom_indices``).
    """
    pos = _normalize_positions(positions_nm)  # (n_frames, n_atoms, 3)
    n_frames, n_atoms, _ = pos.shape

    # Per-state nonbonded parameters (same n_atoms as positions).
    qA, sigA, epsA = _collect_state_nonbonded_params(evb_system, state_A)
    qB, sigB, epsB = _collect_state_nonbonded_params(evb_system, state_B)
    if qA.size != n_atoms or qB.size != n_atoms:
        raise ValueError(
            f"Parameter count ({qA.size}) != n_atoms in positions ({n_atoms}). "
            "Check that the topology and position array match the scaffold."
        )

    solute_set = set(int(i) for i in solute_atom_indices)
    if not solute_set:
        raise ValueError("solute_atom_indices must be non-empty.")
    solute_idx = np.fromiter(sorted(solute_set), dtype=int)

    if lj_combining == "lorentz-berthelot":
        mix_fn = _apply_mixing_lb
    elif lj_combining == "geometric":
        mix_fn = _apply_mixing_geometric
    else:
        raise ValueError(f"Unknown lj_combining '{lj_combining}'.")

    # Per-Q-atom state-specific params (vectors over len(solute_idx)):
    qA_q = qA[solute_idx]
    qB_q = qB[solute_idx]
    sigA_q = sigA[solute_idx]
    sigB_q = sigB[solute_idx]
    epsA_q = epsA[solute_idx]
    epsB_q = epsB[solute_idx]

    residues = _residue_map_from_topology(topology)
    results: list[ResidueContribution] = []

    # Precompute solute positions for all frames: (n_frames, n_solute, 3).
    pos_sol = pos[:, solute_idx, :]

    for resi_idx, resi_name, chain_id, atom_idxs in residues:
        atom_arr = np.asarray(atom_idxs, dtype=int)
        # Skip residues that are entirely solute.
        if set(atom_arr.tolist()) <= solute_set:
            continue
        # Non-solute atoms of this residue.
        non_sol = np.array([a for a in atom_arr if a not in solute_set], dtype=int)
        if non_sol.size == 0:
            continue

        # Residue atom params (state-independent for solvent/non-Q):
        # We use state A params for non-Q atoms; they equal state B params
        # because non-Q atoms are not in the apply_q6 override or reactive set.
        q_r = qA[non_sol]
        sig_r = sigA[non_sol]
        eps_r = epsA[non_sol]

        # Pair expansion: (n_r, n_q) broadcasting.
        # q_pair_A = q_r[:, None] * qA_q[None, :]  (A charges for Q atoms)
        q_pair_A = q_r[:, None] * qA_q[None, :]
        q_pair_B = q_r[:, None] * qB_q[None, :]

        sig_ij_A, eps_ij_A = mix_fn(
            sig_r[:, None], sigA_q[None, :], eps_r[:, None], epsA_q[None, :]
        )
        sig_ij_B, eps_ij_B = mix_fn(
            sig_r[:, None], sigB_q[None, :], eps_r[:, None], epsB_q[None, :]
        )

        # Positions: (n_frames, n_r, 3) and (n_frames, n_q, 3).
        pos_r = pos[:, non_sol, :]  # (n_frames, n_r, 3)
        # Pair distances: (n_frames, n_r, n_q)
        diff = pos_r[:, :, None, :] - pos_sol[:, None, :, :]
        r_nm = np.sqrt((diff ** 2).sum(axis=-1))

        if cutoff_nm is not None:
            mask = r_nm <= cutoff_nm
        else:
            mask = None

        # Coulomb: V = k · q1·q2 / r  (kJ/mol); convert to kcal/mol at the end.
        with np.errstate(divide="ignore", invalid="ignore"):
            v_elec_A = COULOMB_CONSTANT_KJ_NM_PER_E2 * q_pair_A[None, :, :] / r_nm
            v_elec_B = COULOMB_CONSTANT_KJ_NM_PER_E2 * q_pair_B[None, :, :] / r_nm

            # LJ: V = 4ε·[(σ/r)^12 − (σ/r)^6]
            sr6_A = (sig_ij_A[None, :, :] / r_nm) ** 6
            sr6_B = (sig_ij_B[None, :, :] / r_nm) ** 6
            v_vdw_A = 4.0 * eps_ij_A[None, :, :] * (sr6_A ** 2 - sr6_A)
            v_vdw_B = 4.0 * eps_ij_B[None, :, :] * (sr6_B ** 2 - sr6_B)

        # Mask out any zero-distance or cutoff-excluded pairs.
        v_elec_A = np.where(np.isfinite(v_elec_A), v_elec_A, 0.0)
        v_elec_B = np.where(np.isfinite(v_elec_B), v_elec_B, 0.0)
        v_vdw_A = np.where(np.isfinite(v_vdw_A), v_vdw_A, 0.0)
        v_vdw_B = np.where(np.isfinite(v_vdw_B), v_vdw_B, 0.0)
        if mask is not None:
            v_elec_A = v_elec_A * mask
            v_elec_B = v_elec_B * mask
            v_vdw_A = v_vdw_A * mask
            v_vdw_B = v_vdw_B * mask

        # Sum over residue-atom × Q-atom pairs → per-frame ΔE_r in kJ/mol.
        delta_elec_kj = (v_elec_B - v_elec_A).sum(axis=(1, 2))
        delta_vdw_kj = (v_vdw_B - v_vdw_A).sum(axis=(1, 2))

        results.append(
            ResidueContribution(
                residue_index=int(resi_idx),
                residue_name=resi_name,
                chain_id=chain_id,
                atom_indices=tuple(int(a) for a in atom_arr),
                delta_e_elec_kcal_mol=delta_elec_kj / KJ_PER_KCAL,
                delta_e_vdw_kcal_mol=delta_vdw_kj / KJ_PER_KCAL,
                warshel_scaled=False,
            )
        )

    return results


@dataclass
class SoluteAtomContribution:
    """Per-solute-atom decomposition of the solute↔environment ΔE across a
    trajectory (the transpose marginalization of :class:`ResidueContribution`:
    summed over ALL environment atoms, resolved per Q-atom).

    ``phi_env_elec_kcal_per_e`` is the environment Coulomb potential at the
    solute atom, Φ_a = Σ_env k·q_env/r, in kcal/mol per unit charge. It is
    **independent of the solute charge model**, so the electrostatic ΔE of any
    charge set follows by ΔE_elec,a = Δq_a · Φ_a — letting one trajectory score
    the catalytic contribution of arbitrary charge models. Arrays are
    ``(n_frames,)`` in kcal/mol (or kcal/mol/e for the potential).
    """

    solute_atom_index: int
    delta_q_e: float
    phi_env_elec_kcal_per_e: np.ndarray
    delta_e_elec_kcal_mol: np.ndarray
    delta_e_vdw_kcal_mol: np.ndarray

    @property
    def mean_phi(self) -> float:
        return float(self.phi_env_elec_kcal_per_e.mean())

    @property
    def mean_elec(self) -> float:
        return float(self.delta_e_elec_kcal_mol.mean())

    @property
    def mean_vdw(self) -> float:
        return float(self.delta_e_vdw_kcal_mol.mean())


def compute_per_solute_atom_delta_e(
    evb_system: Any,
    topology: Any,
    positions_nm: np.ndarray,
    *,
    solute_atom_indices: Sequence[int],
    state_A: str = "reactant",
    state_B: str = "product",
    cutoff_nm: float | None = None,
    lj_combining: str = "lorentz-berthelot",
    env_atom_indices: Sequence[int] | None = None,
    reaction_field: bool = False,
    rf_dielectric: float = math.inf,
) -> list[SoluteAtomContribution]:
    """Per-solute-atom ΔE and environment potential on a trajectory.

    ``env_atom_indices`` optionally restricts the environment sum to a subset of
    atoms (any that also fall in the solute set are dropped). Use it to resolve
    how a *specific* environment group (e.g. one catalytic residue) couples to
    each solute atom; default ``None`` uses all non-solute atoms.

    Same pairwise Coulomb/LJ physics as :func:`compute_per_residue_delta_e`, but
    summed over every environment (non-solute) atom and resolved **per Q-atom**
    instead of per residue. Returns one entry per solute atom, ordered by sorted
    ``solute_atom_indices``.

    The per-atom electrostatic potential Φ_a is charge-model-independent, so the
    catalytic contribution of any candidate charge set can be scored post hoc as
    Δq_a·Φ_a without re-running dynamics. Summing ``delta_e_elec_kcal_mol`` over
    all returned atoms equals the total solute↔environment ΔE_elec (i.e. the sum
    of the per-residue elec contributions) — the two decompositions are
    consistent marginalizations of the same pair matrix.

    For a **net-charged solute** the bare truncated ``1/r`` sum is only conditionally
    convergent -- its total oscillates with ``cutoff_nm`` and is meaningless. Pass
    ``reaction_field=True`` (with a ``cutoff_nm``) to use the Barker-Watts shifted
    kernel ``1/r + k_rf r^2 - c_rf``, which goes smoothly to zero at the cutoff and
    gives a convergent, decomposable potential. ``rf_dielectric`` sets the continuum
    permittivity beyond the cutoff (``math.inf`` = conductor limit, the default).
    """
    pos = _normalize_positions(positions_nm)
    n_frames, n_atoms, _ = pos.shape

    qA, sigA, epsA = _collect_state_nonbonded_params(evb_system, state_A)
    qB, sigB, epsB = _collect_state_nonbonded_params(evb_system, state_B)
    if qA.size != n_atoms or qB.size != n_atoms:
        raise ValueError(
            f"Parameter count ({qA.size}) != n_atoms in positions ({n_atoms}). "
            "Check that the topology and position array match the scaffold."
        )

    solute_set = set(int(i) for i in solute_atom_indices)
    if not solute_set:
        raise ValueError("solute_atom_indices must be non-empty.")
    solute_idx = np.fromiter(sorted(solute_set), dtype=int)
    n_q = solute_idx.size

    if lj_combining == "lorentz-berthelot":
        mix_fn = _apply_mixing_lb
    elif lj_combining == "geometric":
        mix_fn = _apply_mixing_geometric
    else:
        raise ValueError(f"Unknown lj_combining '{lj_combining}'.")

    rf_krf = rf_crf = None
    if reaction_field:
        if cutoff_nm is None:
            raise ValueError("reaction_field=True requires cutoff_nm (the reaction-field cutoff).")
        rf_krf, rf_crf = _reaction_field_coeffs(cutoff_nm, rf_dielectric)

    qA_q = qA[solute_idx]
    qB_q = qB[solute_idx]
    dq_q = qB_q - qA_q
    sigA_q = sigA[solute_idx]
    sigB_q = sigB[solute_idx]
    epsA_q = epsA[solute_idx]
    epsB_q = epsB[solute_idx]
    pos_sol = pos[:, solute_idx, :]

    # Environment = every atom not in the solute set (optionally restricted).
    if env_atom_indices is None:
        env_pool = range(n_atoms)
    else:
        env_pool = sorted(set(int(i) for i in env_atom_indices))
    env = np.array([a for a in env_pool if a not in solute_set], dtype=int)
    if env.size == 0:
        raise ValueError("no environment atoms left after applying env_atom_indices/solute filter.")
    q_env = qA[env]           # non-solute params are state-independent
    sig_env = sigA[env]
    eps_env = epsA[env]

    # Accumulators over Q-atoms, per frame.
    phi = np.zeros((n_frames, n_q))          # Σ_env k·q_env/r  (kJ/mol/e)
    dvdw = np.zeros((n_frames, n_q))         # Σ_env (vdw_B − vdw_A)  (kJ/mol)

    # Chunk over environment atoms to bound memory (n_frames × chunk × n_q).
    chunk = max(1, int(4_000_000 / max(1, n_frames * n_q)))
    for s in range(0, env.size, chunk):
        e_sl = slice(s, min(s + chunk, env.size))
        pos_e = pos[:, env[e_sl], :]                       # (nf, ne, 3)
        diff = pos_e[:, :, None, :] - pos_sol[:, None, :, :]  # (nf, ne, n_q, 3)
        r_nm = np.sqrt((diff ** 2).sum(axis=-1))           # (nf, ne, n_q)
        mask = (r_nm <= cutoff_nm) if cutoff_nm is not None else None
        with np.errstate(divide="ignore", invalid="ignore"):
            if reaction_field:
                # Barker-Watts shifted kernel -> smoothly 0 at cutoff -> convergent for net-charged solute
                kernel = COULOMB_CONSTANT_KJ_NM_PER_E2 * (1.0 / r_nm + rf_krf * r_nm ** 2 - rf_crf)
            else:
                kernel = COULOMB_CONSTANT_KJ_NM_PER_E2 / r_nm
            phi_c = q_env[None, e_sl, None] * kernel        # (nf, ne, n_q) kJ/mol/e
            sig_ij_A, eps_ij_A = mix_fn(
                sig_env[e_sl][None, :, None], sigA_q[None, None, :],
                eps_env[e_sl][None, :, None], epsA_q[None, None, :])
            sig_ij_B, eps_ij_B = mix_fn(
                sig_env[e_sl][None, :, None], sigB_q[None, None, :],
                eps_env[e_sl][None, :, None], epsB_q[None, None, :])
            sr6_A = (sig_ij_A / r_nm) ** 6
            sr6_B = (sig_ij_B / r_nm) ** 6
            vdw_A = 4.0 * eps_ij_A * (sr6_A ** 2 - sr6_A)
            vdw_B = 4.0 * eps_ij_B * (sr6_B ** 2 - sr6_B)
            dvdw_c = vdw_B - vdw_A
        phi_c = np.where(np.isfinite(phi_c), phi_c, 0.0)
        dvdw_c = np.where(np.isfinite(dvdw_c), dvdw_c, 0.0)
        if mask is not None:
            phi_c = phi_c * mask
            dvdw_c = dvdw_c * mask
        phi += phi_c.sum(axis=1)
        dvdw += dvdw_c.sum(axis=1)

    phi_kcal = phi / KJ_PER_KCAL          # kcal/mol/e
    dvdw_kcal = dvdw / KJ_PER_KCAL
    delec_kcal = dq_q[None, :] * phi_kcal  # ΔE_elec,a = Δq_a · Φ_a

    return [
        SoluteAtomContribution(
            solute_atom_index=int(solute_idx[j]),
            delta_q_e=float(dq_q[j]),
            phi_env_elec_kcal_per_e=phi_kcal[:, j],
            delta_e_elec_kcal_mol=delec_kcal[:, j],
            delta_e_vdw_kcal_mol=dvdw_kcal[:, j],
        )
        for j in range(n_q)
    ]


@dataclass
class LRAResult:
    """LRA and reorganization energy decomposition per residue.

    All values are kcal/mol.
    """

    residue_index: int
    residue_name: str
    chain_id: str
    lra_elec: float
    lra_vdw: float
    reorg_elec: float
    reorg_vdw: float
    warshel_scaled: bool
    n_frames_A: int
    n_frames_B: int

    @property
    def lra_total(self) -> float:
        return self.lra_elec + self.lra_vdw

    @property
    def reorg_total(self) -> float:
        return self.reorg_elec + self.reorg_vdw


def lra_per_residue(
    contributions_state_A: Sequence[ResidueContribution],
    contributions_state_B: Sequence[ResidueContribution],
    *,
    apply_warshel_scaling: bool = False,
) -> list[LRAResult]:
    """Combine two-ensemble per-residue ΔE into LRA and reorganization.

    ``contributions_state_A`` comes from sampling at state A (reactant basin),
    ``contributions_state_B`` from state B (product basin). Residues are
    matched by ``(chain_id, residue_index)``.

    Parameters
    ----------
    apply_warshel_scaling :
        If True, multiply electrostatic LRA for ionized residues (ARG, LYS,
        HIP, GLU, ASP) by 0.5 following Warshel 1997 doi:10.1021/jp962478o.
        Off by default because the scaling is a heuristic correction that
        not every workflow applies.

    Notes
    -----
    For identity reactions ΔG₀ = 0 and the LRA estimate of ΔG‡ is
    directly sum_r(LRA_r). For asymmetric reactions both LRA and the
    reorganization term play a role; see Hwang & Warshel 1987,
    J. Am. Chem. Soc. 109, 715.
    """
    by_key_A = {(c.chain_id, c.residue_index): c for c in contributions_state_A}
    by_key_B = {(c.chain_id, c.residue_index): c for c in contributions_state_B}
    keys = sorted(set(by_key_A) & set(by_key_B), key=lambda k: (k[0], k[1]))

    results: list[LRAResult] = []
    for key in keys:
        cA = by_key_A[key]
        cB = by_key_B[key]
        assert cA.residue_name == cB.residue_name, (
            f"Residue name mismatch at {key}: {cA.residue_name} vs {cB.residue_name}"
        )
        mean_A_elec = cA.mean_elec
        mean_A_vdw = cA.mean_vdw
        mean_B_elec = cB.mean_elec
        mean_B_vdw = cB.mean_vdw
        lra_elec = 0.5 * (mean_A_elec + mean_B_elec)
        lra_vdw = 0.5 * (mean_A_vdw + mean_B_vdw)
        reorg_elec = mean_A_elec - lra_elec
        reorg_vdw = mean_A_vdw - lra_vdw

        scaled = False
        if apply_warshel_scaling and cA.residue_name in IONIZED_RESIDUE_NAMES:
            lra_elec *= WARSHEL_SCALING
            reorg_elec *= WARSHEL_SCALING
            scaled = True

        results.append(
            LRAResult(
                residue_index=cA.residue_index,
                residue_name=cA.residue_name,
                chain_id=cA.chain_id,
                lra_elec=lra_elec,
                lra_vdw=lra_vdw,
                reorg_elec=reorg_elec,
                reorg_vdw=reorg_vdw,
                warshel_scaled=scaled,
                n_frames_A=cA.delta_e_elec_kcal_mol.size,
                n_frames_B=cB.delta_e_elec_kcal_mol.size,
            )
        )
    return results


def summarize_lra(results: Iterable[LRAResult], *, top_n: int = 10) -> str:
    """Human-readable summary of the top-N largest-|LRA| contributions."""
    items = list(results)
    items.sort(key=lambda r: -abs(r.lra_total))
    lines = ["# Per-residue LRA contributions (kcal/mol), sorted by |LRA_total|:"]
    lines.append(f"# {'chain':>5} {'resi':>5} {'name':>4} {'LRA_el':>8} "
                 f"{'LRA_vdw':>8} {'LRA_tot':>8} {'REORG_tot':>9}")
    for r in items[:top_n]:
        mark = "*" if r.warshel_scaled else " "
        lines.append(
            f"{mark} {r.chain_id:>5} {r.residue_index:>5} {r.residue_name:>4} "
            f"{r.lra_elec:>+8.3f} {r.lra_vdw:>+8.3f} {r.lra_total:>+8.3f} "
            f"{r.reorg_total:>+9.3f}"
        )
    lines.append("# (* = Warshel-scaled ionized residue)")
    return "\n".join(lines)


__all__ = [
    "ResidueContribution",
    "LRAResult",
    "IONIZED_RESIDUE_NAMES",
    "WARSHEL_SCALING",
    "compute_per_residue_delta_e",
    "lra_per_residue",
    "summarize_lra",
]
