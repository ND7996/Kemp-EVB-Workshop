"""Structural / geometric analysis utilities for EVB trajectories.

A small pure-numpy replacement for the subset of ``Qpyl.core.qcalc``
functionality that's useful for EVB post-processing:

* ``compute_distances(positions_nm, atom_pairs)``
* ``compute_angles(positions_nm, atom_triplets)``
* ``compute_torsions(positions_nm, atom_quartets)``
* ``compute_rmsd(positions_nm, reference_positions_nm, atom_indices=...)``
* ``compute_rmsf(positions_nm, reference_positions_nm=..., atom_indices=...)``
* ``compute_rdf(positions_nm, group_a_indices, group_b_indices, ...)``
* ``evb_reaction_coordinate(positions_nm, central, leaving, attacking)``
  — common antisymmetric reaction coordinate
  ``r(central−leaving) − r(central−attacking)`` for SN2-like reactions.

All inputs are position arrays in nanometers with shape ``(n_frames, n_atoms, 3)``
(or ``(n_atoms, 3)`` for a single frame) — matching what mdtraj or OpenMM
``State.getPositions()`` produces. No periodic wrapping: assumes a non-
periodic spherical (SCAAS) system, which is the EVB default. If your system
is periodic, pre-wrap positions with mdtraj before calling.

RMSD alignment uses the standard Kabsch algorithm (SVD-based closed-form
solution). No mdtraj dependency is imposed on the library.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def _normalize_positions(positions_nm: np.ndarray) -> np.ndarray:
    """Accept ``(n_atoms, 3)`` or ``(n_frames, n_atoms, 3)`` -> always 3-D."""
    a = np.asarray(positions_nm, dtype=float)
    if a.ndim == 2:
        a = a[np.newaxis, :, :]
    elif a.ndim != 3:
        raise ValueError(
            "positions_nm must have shape (n_atoms, 3) or (n_frames, n_atoms, 3); "
            f"got {a.shape}."
        )
    return a


def compute_distances(positions_nm: np.ndarray,
                      atom_pairs: Sequence[tuple[int, int]]) -> np.ndarray:
    """Per-frame distances for a list of atom-index pairs, in nanometers.

    Returns shape ``(n_frames, n_pairs)``.
    """
    pos = _normalize_positions(positions_nm)
    pairs = np.asarray(atom_pairs, dtype=int)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError(f"atom_pairs must be (n_pairs, 2); got {pairs.shape}.")
    diff = pos[:, pairs[:, 0], :] - pos[:, pairs[:, 1], :]
    return np.sqrt((diff ** 2).sum(axis=-1))


def compute_angles(positions_nm: np.ndarray,
                   atom_triplets: Sequence[tuple[int, int, int]]) -> np.ndarray:
    """Per-frame angles (i-j-k with j as vertex), in radians.

    Returns shape ``(n_frames, n_triplets)``.
    """
    pos = _normalize_positions(positions_nm)
    trips = np.asarray(atom_triplets, dtype=int)
    if trips.ndim != 2 or trips.shape[1] != 3:
        raise ValueError(f"atom_triplets must be (n_triplets, 3); got {trips.shape}.")
    v1 = pos[:, trips[:, 0], :] - pos[:, trips[:, 1], :]
    v2 = pos[:, trips[:, 2], :] - pos[:, trips[:, 1], :]
    cos = (v1 * v2).sum(axis=-1) / (
        np.linalg.norm(v1, axis=-1) * np.linalg.norm(v2, axis=-1)
    )
    return np.arccos(np.clip(cos, -1.0, 1.0))


def compute_torsions(positions_nm: np.ndarray,
                     atom_quartets: Sequence[tuple[int, int, int, int]]) -> np.ndarray:
    """Per-frame torsions / proper dihedrals (i-j-k-l), in radians on (−π, π].

    Follows the IUPAC convention; sign matches OpenMM's
    ``PeriodicTorsionForce`` and mdtraj's ``compute_dihedrals``.
    Returns shape ``(n_frames, n_quartets)``.
    """
    pos = _normalize_positions(positions_nm)
    quads = np.asarray(atom_quartets, dtype=int)
    if quads.ndim != 2 or quads.shape[1] != 4:
        raise ValueError(f"atom_quartets must be (n_quartets, 4); got {quads.shape}.")
    p0 = pos[:, quads[:, 0], :]
    p1 = pos[:, quads[:, 1], :]
    p2 = pos[:, quads[:, 2], :]
    p3 = pos[:, quads[:, 3], :]

    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2

    # Normal to plane (p0, p1, p2) and (p1, p2, p3).
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    # m1 = n1 × (b2 normalised) — used for signed angle.
    b2_norm = b2 / np.linalg.norm(b2, axis=-1, keepdims=True)
    m1 = np.cross(n1, b2_norm)
    x = (n1 * n2).sum(axis=-1)
    y = (m1 * n2).sum(axis=-1)
    return np.arctan2(y, x)


def _kabsch_align(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Kabsch rotation + translation to best-fit ``mobile`` onto ``reference``.

    Both inputs have shape ``(n_atoms, 3)``. Returns the aligned mobile
    coordinates and the applied 3x3 rotation matrix.
    """
    mob_c = mobile - mobile.mean(axis=0)
    ref_c = reference - reference.mean(axis=0)
    H = mob_c.T @ ref_c
    U, _s, Vt = np.linalg.svd(H)
    # Correct for reflection (ensure proper rotation, det(R) = +1).
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    aligned = mob_c @ R.T + reference.mean(axis=0)
    return aligned, R


def compute_rmsd(positions_nm: np.ndarray,
                 reference_positions_nm: np.ndarray,
                 atom_indices: Sequence[int] | None = None) -> np.ndarray:
    """Per-frame RMSD vs a reference after Kabsch alignment, in nanometers.

    If ``atom_indices`` is given, alignment and RMSD are computed on that
    subset (e.g. solute only). Returns shape ``(n_frames,)``.
    """
    pos = _normalize_positions(positions_nm)
    ref = np.asarray(reference_positions_nm, dtype=float)
    if ref.ndim == 3:
        if ref.shape[0] != 1:
            raise ValueError("reference_positions_nm must be a single frame.")
        ref = ref[0]
    if ref.ndim != 2 or ref.shape[1] != 3:
        raise ValueError(f"reference must be (n_atoms, 3); got {ref.shape}.")
    if atom_indices is not None:
        idx = np.asarray(atom_indices, dtype=int)
        pos = pos[:, idx, :]
        ref = ref[idx, :]
    n_frames = pos.shape[0]
    rmsd = np.empty(n_frames, dtype=float)
    for f in range(n_frames):
        aligned, _ = _kabsch_align(pos[f], ref)
        rmsd[f] = float(np.sqrt(((aligned - ref) ** 2).sum(axis=-1).mean()))
    return rmsd


def compute_rmsf(positions_nm: np.ndarray,
                 reference_positions_nm: np.ndarray | None = None,
                 atom_indices: Sequence[int] | None = None) -> np.ndarray:
    """Per-atom root-mean-square fluctuation over the trajectory.

    Aligns every frame to ``reference_positions_nm`` (default: the trajectory
    mean structure) and returns ``sqrt(⟨(x − x̄)²⟩)`` per atom, in nm.
    Shape: ``(n_atoms,)`` (or ``(len(atom_indices),)``).
    """
    pos = _normalize_positions(positions_nm)
    if atom_indices is not None:
        pos = pos[:, np.asarray(atom_indices, dtype=int), :]
    if reference_positions_nm is None:
        ref = pos.mean(axis=0)
    else:
        ref = np.asarray(reference_positions_nm, dtype=float)
        if ref.ndim == 3:
            ref = ref[0]
        if atom_indices is not None:
            ref = ref[np.asarray(atom_indices, dtype=int), :]
    aligned = np.empty_like(pos)
    for f in range(pos.shape[0]):
        aligned[f], _ = _kabsch_align(pos[f], ref)
    fluct = aligned - aligned.mean(axis=0)  # (n_frames, n_atoms, 3)
    return np.sqrt((fluct ** 2).sum(axis=-1).mean(axis=0))


def compute_rdf(positions_nm: np.ndarray,
                group_a_indices: Sequence[int],
                group_b_indices: Sequence[int],
                *,
                r_max_nm: float = 1.0,
                n_bins: int = 100,
                volume_nm3: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Radial distribution function g(r) between two atom groups.

    Sums pair-distance histograms across all frames, normalises by shell
    volume ``4·π·r²·Δr`` and by ``N_a·ρ_b`` where ``ρ_b = len(group_b)/V``.
    If ``volume_nm3`` is None, uses the minimum bounding sphere of all
    positions across the trajectory (appropriate for SCAAS spheres).

    Returns ``(r_centers, g_r)`` arrays of length ``n_bins``.
    """
    pos = _normalize_positions(positions_nm)
    A = np.asarray(group_a_indices, dtype=int)
    B = np.asarray(group_b_indices, dtype=int)
    if A.size == 0 or B.size == 0:
        raise ValueError("group_a_indices and group_b_indices must be non-empty.")

    edges = np.linspace(0.0, r_max_nm, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    counts = np.zeros(n_bins, dtype=float)

    for frame in pos:
        diff = frame[A, np.newaxis, :] - frame[np.newaxis, B, :]
        r = np.sqrt((diff ** 2).sum(axis=-1))
        # Exclude self-pairs (same atom in both groups) by zeroing r=0 entries.
        r = r[r > 0]
        h, _ = np.histogram(r, bins=edges)
        counts += h

    n_frames = pos.shape[0]
    if volume_nm3 is None:
        # Bounding-sphere volume of all positions across trajectory.
        centroid = pos.mean(axis=(0, 1))
        r_max_atom = np.linalg.norm(pos - centroid, axis=-1).max()
        volume_nm3 = (4.0 / 3.0) * np.pi * r_max_atom ** 3

    shell_vol = 4.0 * np.pi * centers ** 2 * (edges[1:] - edges[:-1])
    rho_b = float(B.size) / float(volume_nm3)
    ideal_counts_per_frame = float(A.size) * rho_b * shell_vol
    with np.errstate(divide="ignore", invalid="ignore"):
        g_r = counts / (ideal_counts_per_frame * n_frames)
    g_r = np.where(np.isfinite(g_r), g_r, 0.0)
    return centers, g_r


def evb_reaction_coordinate(positions_nm: np.ndarray,
                            central_index: int,
                            leaving_index: int,
                            attacking_index: int) -> np.ndarray:
    """Antisymmetric SN2 reaction coordinate:
    ``RC(t) = r(central−leaving) − r(central−attacking)``.

    Negative in reactant basin, 0 at the symmetric TS, positive in product
    basin. Returns shape ``(n_frames,)`` in nanometers.
    """
    d = compute_distances(positions_nm,
                          [(central_index, leaving_index),
                           (central_index, attacking_index)])
    return d[:, 0] - d[:, 1]


__all__ = [
    "compute_distances",
    "compute_angles",
    "compute_torsions",
    "compute_rmsd",
    "compute_rmsf",
    "compute_rdf",
    "evb_reaction_coordinate",
]
