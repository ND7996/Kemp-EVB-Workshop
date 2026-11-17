from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Dict, Set, Tuple

from ..types import AtomReference


def find_unmapped_atoms(
    *,
    state_names: Sequence[str],
    state_atom_references: Callable[[str], Tuple[AtomReference, ...]],
    normalize_atom_reference: Callable[[AtomReference], AtomReference],
    mapping: Mapping[AtomReference, AtomReference],
) -> Dict[str, Tuple[AtomReference, ...]]:
    all_atoms: Dict[str, Tuple[AtomReference, ...]] = {}
    for state_name in state_names:
        all_atoms[state_name] = tuple(state_atom_references(state_name))

    referenced: Dict[str, Set[AtomReference]] = {state: set() for state in state_names}
    for left_ref, right_ref in mapping.items():
        left_norm = normalize_atom_reference(left_ref)
        right_norm = normalize_atom_reference(right_ref)
        referenced[left_norm[0]].add(left_norm)
        referenced[right_norm[0]].add(right_norm)

    unmapped: Dict[str, Tuple[AtomReference, ...]] = {}
    for state_name, identifiers in all_atoms.items():
        missing = tuple(sorted(set(identifiers) - referenced[state_name]))
        if missing:
            unmapped[state_name] = missing
    return unmapped


def validate_atom_mapping_atoms(
    mapping: Mapping[AtomReference, AtomReference],
    *,
    normalize_atom_reference: Callable[[AtomReference], AtomReference],
) -> None:
    if not mapping:
        return

    checked: Set[AtomReference] = set()
    for left_ref, right_ref in mapping.items():
        for ref in (left_ref, right_ref):
            if ref in checked:
                continue
            normalize_atom_reference(ref)
            checked.add(ref)


def normalize_atom_mapping(
    mapping: Mapping[AtomReference, AtomReference],
    *,
    state_names: Sequence[str],
    state_atom_references: Callable[[str], Tuple[AtomReference, ...]],
    normalize_atom_reference: Callable[[AtomReference], AtomReference],
) -> Dict[AtomReference, AtomReference]:
    normalized: Dict[AtomReference, AtomReference] = {}

    def _record_pair(left: AtomReference, right: AtomReference) -> None:
        if left in normalized and normalized[left] != right:
            raise ValueError(f"Atom {left!r} is mapped to multiple targets ({normalized[left]!r}, {right!r}).")
        normalized[left] = right

    for left_ref, right_ref in mapping.items():
        left_norm = normalize_atom_reference(left_ref)
        right_norm = normalize_atom_reference(right_ref)
        _record_pair(left_norm, right_norm)
        _record_pair(right_norm, left_norm)

    for state_name in state_names:
        for atom_ref in state_atom_references(state_name):
            if atom_ref not in normalized:
                normalized[atom_ref] = atom_ref

    return normalized


def build_identity_atom_mapping(
    *,
    state_names: Sequence[str],
    state_atom_references: Callable[[str], Tuple[AtomReference, ...]],
) -> Dict[AtomReference, AtomReference]:
    mapping: Dict[AtomReference, AtomReference] = {}
    for state_name in state_names:
        for atom_ref in state_atom_references(state_name):
            mapping[atom_ref] = atom_ref
    return mapping


def build_atom_canonical_labels(
    atom_mapping: Mapping[AtomReference, AtomReference],
) -> Dict[AtomReference, AtomReference]:
    adjacency: Dict[AtomReference, Set[AtomReference]] = {}
    for source, target in atom_mapping.items():
        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set()).add(source)

    visited: Set[AtomReference] = set()
    labels: Dict[AtomReference, AtomReference] = {}

    for atom in adjacency:
        if atom in visited:
            continue
        stack = [atom]
        component: list[AtomReference] = []
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            stack.extend(adjacency.get(current, ()))
        representative = min(component)
        for member in component:
            labels[member] = representative
    return labels


__all__ = [
    "find_unmapped_atoms",
    "validate_atom_mapping_atoms",
    "normalize_atom_mapping",
    "build_identity_atom_mapping",
    "build_atom_canonical_labels",
]
