"""Utilities for loading PDB files needed to configure EVB/OpenMM simulations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import FrozenSet, Tuple
import os

from openmm import app

PathLike = os.PathLike[str] | str


@dataclass(frozen=True)
class EVBState:
    """Holds the topology/positions of a single EVB state sourced from a PDB."""

    path: Path
    topology: app.Topology
    positions: Tuple[app.Vec3, ...]

    @property
    def n_atoms(self) -> int:
        return len(self.positions)

    @property
    def residue_names(self) -> FrozenSet[str]:
        """Unique residue identifiers present in the state."""
        names: set[str] = set()
        for residue in self.topology.residues():
            name = (residue.name or "").strip()
            if not name:
                name = f"RES{residue.index}"
            names.add(name)
        return frozenset(names)


def read_pdb_state(pdb_path: PathLike) -> EVBState:
    """Read a single structure file into an EVBState structure.

    Accepts PDB and PDBx/mmCIF (``.cif``). mmCIF is used by the prepared-bundle
    serializer for large solvated systems, whose residue ids exceed the PDB 4-char
    limit (>9999) and would not survive a PDB round-trip.

    OpenMM's ``PDBFile`` reader normalizes alternate protonation residue names on
    load (e.g. ``ASH``→``ASP``, ``HIE``/``HIP``→``HIS``). For ordinary use that is
    harmless, but in proton-transfer EVB the *per-state protonation is the reactive
    change itself* -- the reactant and product diabats deliberately differ in which
    residue carries the proton (``ASH`` vs ``ASP``, ``HIE`` vs ``HIP``). Normalization
    would collapse both states to the same template and erase that difference. We
    therefore restore the on-disk residue names for ``.pdb`` inputs so the force field
    assigns the correct per-state templates.
    """
    path = _validate_path(pdb_path)
    is_cif = str(path).lower().endswith((".cif", ".pdbx"))
    try:
        pdb_file = app.PDBxFile(str(path)) if is_cif else app.PDBFile(str(path))
    except (OSError, ValueError) as exc:  # pragma: no cover
        raise ValueError(f"Failed to parse structure file: {path}") from exc

    topology = pdb_file.topology
    if not is_cif:
        _restore_on_disk_residue_names(path, topology)
    else:
        # mmCIF (prepared-bundle) loses the on-disk protonation distinction OpenMM normalized
        # (a protonated Asp -> "ASP" + HD2). Re-derive ASH/GLH from atom content so per-state
        # protonation survives the bundle round-trip (the .pdb path uses on-disk names instead).
        _apply_protonation_variants(topology)

    positions = tuple(pdb_file.positions)
    return EVBState(path=path, topology=topology, positions=positions)


def _apply_protonation_variants(topology: app.Topology) -> None:
    """Rename protonated acids (ASP+HD2 -> ASH, GLU+HE2 -> GLH) inferred from atom content."""
    from ..core.naming import protonation_variant_name

    for residue in topology.residues():
        new_name = protonation_variant_name(residue.name, [atom.name for atom in residue.atoms()])
        if new_name != residue.name:
            residue.name = new_name


def _raw_pdb_residue_names(path: Path) -> dict[tuple[str, str], str]:
    """Map ``(chain_id, resSeq+iCode) -> residue name`` straight from the ATOM/HETATM
    records, before any OpenMM normalization.

    Rejects four-character residue names rather than truncating them. The PDB resName field is
    columns 18-20, so a name written into 18-21 (e.g. ``GLUP``) is read here as its three-char
    prefix (``GLU``) -- which would then be *restored* onto the residue, silently renaming a
    protonated glutamate into an ordinary one and colliding it with every other GLU in the chain.
    OpenMM is no kinder: a four-character name has no entry in ``residues.xml``, so its residue
    is built with no bonds at all. Use the three-character variant name instead (``GLH``/``ASH``,
    which OpenMM already knows as alternates).
    """
    out: dict[tuple[str, str], str] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21:22].strip()
                res_id = line[22:26].strip() + line[26:27].strip()  # resSeq + insertion code
                if len(line) > 20 and line[20:21] not in ("", " ", "\n"):
                    raise ValueError(
                        f"{path}: line {line_number} carries a four-character residue name "
                        f"'{line[17:21].strip()}'. The PDB resName field is columns 18-20; a "
                        f"longer name is silently truncated on read (to "
                        f"'{line[17:20].strip()}') and loses all of its standard bonds. Write "
                        f"the three-character variant name instead."
                    )
                out[(chain, res_id)] = line[17:20].strip()
    return out


def _restore_on_disk_residue_names(path: Path, topology: app.Topology) -> None:
    """Re-apply the residue names exactly as written on disk, undoing OpenMM's
    protonation-variant normalization. Matches residues by (chain id, residue id)."""
    raw = _raw_pdb_residue_names(path)
    if not raw:
        return
    for residue in topology.residues():
        key = (residue.chain.id, str(residue.id).strip())
        disk_name = raw.get(key)
        if disk_name and disk_name != residue.name:
            residue.name = disk_name


def read_evb_pdbs(state_a_path: PathLike, state_b_path: PathLike) -> Tuple[EVBState, EVBState]:
    """Read two PDB files (e.g., initial/final states) and return their EVB states."""
    return (
        read_pdb_state(state_a_path),
        read_pdb_state(state_b_path),
    )


def _validate_path(pdb_path: PathLike) -> Path:
    path = Path(pdb_path).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"PDB file not found: {path}")
    return path


__all__ = [
    "EVBState",
    "read_pdb_state",
    "read_evb_pdbs",
]
