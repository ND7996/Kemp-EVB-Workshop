from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Dict, Tuple

from openmm import app

from ..models import ResidueTarget
from ..types import CanonicalPairKey

ResidueSignature = Tuple[Tuple[Tuple[str, str | None], ...], Tuple[Tuple[str, str], ...]]


def ensure_residue_name_consistency(
    *,
    residue_targets: Mapping[str, list[ResidueTarget]],
    residue_signature: Callable[[ResidueTarget, bool], ResidueSignature],
    format_signature_difference: Callable[[ResidueSignature, ResidueSignature], str],
    include_bonds: bool = True,
    bonds_optional_residues: "frozenset[str] | set[str]" = frozenset(),
    skip_residues: "frozenset[str] | set[str]" = frozenset(),
) -> None:
    """Validate that same-named residues are chemically identical across states.

    ``bonds_optional_residues`` lists residues whose bonds may legitimately differ
    between diabatic states (per-state EVB ligand parameterization, where the
    reactive bonds are made/broken). For those, only atom identity is checked.

    ``skip_residues`` lists residues exempted from the check entirely. Use this for
    force-field-standard residues (e.g. standard amino acids), whose same-named
    instances legitimately differ in atom content across positions -- N-/C-terminal
    vs internal variants -- and are parameterized correctly by the force field's own
    terminal templates. The check is meant for custom (e.g. OpenFF) reactive species.
    """

    signatures: Dict[str, ResidueSignature] = {}
    reference_states: Dict[str, str] = {}
    for normalized_name, contexts in residue_targets.items():
        if normalized_name in skip_residues:
            continue
        check_bonds = include_bonds and normalized_name not in bonds_optional_residues
        for target in contexts:
            signature = residue_signature(target, check_bonds)
            if normalized_name not in signatures:
                signatures[normalized_name] = signature
                reference_states[normalized_name] = target.state_name
                continue
            if signature != signatures[normalized_name]:
                display_name = target.residue.name.strip() or normalized_name
                ref_state = reference_states[normalized_name]
                diff_detail = format_signature_difference(signatures[normalized_name], signature)
                raise ValueError(
                    f"Residue '{display_name}' differs between states '{ref_state}' and '{target.state_name}'. "
                    "Please assign unique residue names for chemically distinct species. "
                    f"Differences: {diff_detail}"
                )


def residue_signature(
    *,
    target: ResidueTarget,
    topology_bonds: Mapping[str, Tuple[Tuple[app.Topology.Atom, app.Topology.Atom], ...]],
    include_bonds: bool,
) -> ResidueSignature:
    atoms = tuple(
        sorted(
            (atom.name.strip(), atom.element.symbol if atom.element else None) for atom in target.residue.atoms()
        )
    )
    bonds: Tuple[Tuple[str, str], ...]
    if include_bonds:
        bonds_entries = []
        for atom1, atom2 in topology_bonds[target.state_name]:
            if atom1.residue == target.residue and atom2.residue == target.residue:
                bonds_entries.append(tuple(sorted((atom1.name.strip(), atom2.name.strip()))))
        bonds = tuple(sorted(bonds_entries))
    else:
        bonds = tuple()
    return atoms, bonds


def residue_bond_pairs(
    *,
    target: ResidueTarget,
    topology_bonds: Mapping[str, Tuple[Tuple[app.Topology.Atom, app.Topology.Atom], ...]],
    atom_reference_from_topology_atom: Callable[[str, app.Topology.Atom], tuple[str, str, int, str]],
    canonical_pair_key: Callable[[tuple[str, str, int, str], tuple[str, str, int, str]], CanonicalPairKey],
) -> Tuple[CanonicalPairKey, ...]:
    state_name = target.state_name
    residue = target.residue
    bonds: list[CanonicalPairKey] = []
    for atom1, atom2 in topology_bonds[state_name]:
        if atom1.residue is residue and atom2.residue is residue:
            ref1 = atom_reference_from_topology_atom(state_name, atom1)
            ref2 = atom_reference_from_topology_atom(state_name, atom2)
            bonds.append(canonical_pair_key(ref1, ref2))
    bonds.sort()
    return tuple(bonds)


def format_signature_difference(reference: ResidueSignature, candidate: ResidueSignature) -> str:
    ref_atoms, ref_bonds = reference
    cand_atoms, cand_bonds = candidate

    def _fmt_atom(atom: Tuple[str, str | None]) -> str:
        name, element = atom
        return f"{name}({element or '?'})"

    def _fmt_bond(bond: Tuple[str, str]) -> str:
        return f"{bond[0]}-{bond[1]}"

    missing_atoms = sorted(set(ref_atoms) - set(cand_atoms))
    extra_atoms = sorted(set(cand_atoms) - set(ref_atoms))
    missing_bonds = sorted(set(ref_bonds) - set(cand_bonds))
    extra_bonds = sorted(set(cand_bonds) - set(ref_bonds))

    def _summarize(entries, formatter):
        if not entries:
            return ""
        formatted = [formatter(entry) for entry in entries[:5]]
        if len(entries) > 5:
            formatted.append("...")
        return ", ".join(formatted)

    parts = []
    if missing_atoms:
        parts.append(f"missing atoms [{_summarize(missing_atoms, _fmt_atom)}]")
    if extra_atoms:
        parts.append(f"extra atoms [{_summarize(extra_atoms, _fmt_atom)}]")
    if missing_bonds:
        parts.append(f"missing bonds [{_summarize(missing_bonds, _fmt_bond)}]")
    if extra_bonds:
        parts.append(f"extra bonds [{_summarize(extra_bonds, _fmt_bond)}]")

    if not parts:
        return "unknown difference (signatures differ but no atom/bond mismatch detected)."
    return "; ".join(parts)


__all__ = [
    "ResidueSignature",
    "ensure_residue_name_consistency",
    "residue_signature",
    "residue_bond_pairs",
    "format_signature_difference",
]
