from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Dict, Set, Tuple

from openmm import app, unit

from ..input_reader.system_reader import EVBState
from ..types import AtomReference, CanonicalAtomKey, CanonicalPairKey


def normalize_atom_reference(
    ref: AtomReference,
    *,
    states: Mapping[str, EVBState],
    state_atom_identifier_cache: Dict[str, Set[AtomReference]],
    state_atom_references: Callable[[str], Tuple[AtomReference, ...]],
) -> AtomReference:
    if not isinstance(ref, tuple) or len(ref) != 4:
        raise TypeError(
            "Atom references must be 4-tuples of (state_name, chain_id, residue_id, atom_name)."
        )
    state_name, chain_id, residue_id, atom_name = ref
    if state_name not in states:
        raise ValueError(f"Atom mapping refers to unknown state '{state_name}'.")
    if not isinstance(chain_id, str) or not chain_id.strip():
        raise ValueError("Atom mapping chain identifiers must be non-empty strings.")
    cleaned_chain = chain_id.strip()
    try:
        residue_number = int(residue_id)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Residue id '{residue_id}' for state '{state_name}' is not a valid integer."
        ) from exc
    if not isinstance(atom_name, str) or not atom_name.strip():
        raise ValueError("Atom mapping atom names must be non-empty strings.")
    cleaned_atom = atom_name.strip()

    normalized = (state_name, cleaned_chain, residue_number, cleaned_atom)
    if state_name not in state_atom_identifier_cache:
        state_atom_identifier_cache[state_name] = set(state_atom_references(state_name))
    if normalized not in state_atom_identifier_cache[state_name]:
        raise ValueError(f"Atom reference {normalized!r} does not exist in state '{state_name}'.")
    return normalized


def canonical_pair_key(first: AtomReference, second: AtomReference) -> CanonicalPairKey:
    return (first, second) if first <= second else (second, first)


def atom_reference_from_topology_atom(
    *,
    state_name: str,
    atom: app.Topology.Atom,
    chain_identifier: Callable[[app.Topology.Chain], str],
    residue_sequence_number: Callable[[app.Topology.Residue, str], int],
) -> AtomReference:
    residue = atom.residue
    chain = residue.chain
    chain_id = chain_identifier(chain)
    residue_code = residue_sequence_number(residue, state_name)
    atom_name = (atom.name or "").strip()
    if not atom_name:
        raise ValueError(
            f"Atom without a name encountered in residue '{residue.name}' (state '{state_name}')."
        )
    return (state_name, chain_id, residue_code, atom_name)


def state_atom_references(
    *,
    state_name: str,
    states: Mapping[str, EVBState],
    state_atom_reference_lists: Dict[str, Tuple[AtomReference, ...]],
    state_atom_identifier_cache: Dict[str, Set[AtomReference]],
    atom_reference_from_topology_atom: Callable[[str, app.Topology.Atom], AtomReference],
) -> Tuple[AtomReference, ...]:
    cached = state_atom_reference_lists.get(state_name)
    if cached is not None:
        return cached

    state = states[state_name]
    atoms: list[AtomReference] = []
    for atom in state.topology.atoms():
        atoms.append(atom_reference_from_topology_atom(state_name, atom))
    references = tuple(atoms)
    state_atom_reference_lists[state_name] = references
    state_atom_identifier_cache[state_name] = set(references)
    return references


def state_particle_masses(
    *,
    state_name: str,
    modellers: Mapping[str, app.Modeller],
    state_particle_mass_lists: Dict[str, Tuple[unit.Quantity, ...]],
) -> Tuple[unit.Quantity, ...]:
    cached = state_particle_mass_lists.get(state_name)
    if cached is not None:
        return cached

    modeller = modellers[state_name]
    masses: list[unit.Quantity] = []
    for atom in modeller.topology.atoms():
        element = atom.element
        if element is None or element.mass is None:
            masses.append(0.0 * unit.dalton)
        else:
            masses.append(element.mass)
    stored = tuple(masses)
    state_particle_mass_lists[state_name] = stored
    return stored


def atom_reference_from_index(
    *,
    state_name: str,
    atom_index: int,
    state_atom_references: Callable[[str], Tuple[AtomReference, ...]],
) -> AtomReference:
    references = state_atom_references(state_name)
    try:
        return references[atom_index]
    except IndexError as exc:  # pragma: no cover - defensive
        raise IndexError(
            f"Atom index {atom_index} out of range for state '{state_name}' (n_atoms={len(references)})."
        ) from exc


def canonical_atom_label(
    *,
    state_name: str,
    atom_index: int,
    atom_reference_from_index: Callable[[str, int], AtomReference],
    atom_canonical_labels: Mapping[AtomReference, AtomReference],
) -> AtomReference:
    reference = atom_reference_from_index(state_name, atom_index)
    return atom_canonical_labels.get(reference, reference)


def canonical_bond_key_from_indices(
    *,
    state_name: str,
    atom1: int,
    atom2: int,
    canonical_atom_label: Callable[[str, int], AtomReference],
) -> CanonicalAtomKey:
    labels = tuple(canonical_atom_label(state_name, idx) for idx in (atom1, atom2))
    return tuple(sorted(labels))


def canonical_angle_key_from_indices(
    *,
    state_name: str,
    atom1: int,
    atom2: int,
    atom3: int,
    canonical_atom_label: Callable[[str, int], AtomReference],
) -> CanonicalAtomKey:
    labels = tuple(canonical_atom_label(state_name, idx) for idx in (atom1, atom2, atom3))
    if labels[0] <= labels[2]:
        return labels
    return (labels[2], labels[1], labels[0])


def canonical_torsion_key_from_indices(
    *,
    state_name: str,
    atom1: int,
    atom2: int,
    atom3: int,
    atom4: int,
    canonical_atom_label: Callable[[str, int], AtomReference],
) -> CanonicalAtomKey:
    labels = tuple(canonical_atom_label(state_name, idx) for idx in (atom1, atom2, atom3, atom4))
    reverse = tuple(reversed(labels))
    return labels if labels <= reverse else reverse


def state_canonical_label_index_map(
    *,
    state_name: str,
    state_canonical_index_maps: Dict[str, Dict[AtomReference, int]],
    state_atom_references: Callable[[str], Tuple[AtomReference, ...]],
    canonical_atom_label: Callable[[str, int], AtomReference],
) -> Mapping[AtomReference, int]:
    cached = state_canonical_index_maps.get(state_name)
    if cached is not None:
        return cached

    references = state_atom_references(state_name)
    mapping: Dict[AtomReference, int] = {}
    for index, _ in enumerate(references):
        label = canonical_atom_label(state_name, index)
        mapping[label] = index
    state_canonical_index_maps[state_name] = mapping
    return mapping


def particle_indices_from_canonical_key(
    *,
    state_name: str,
    canonical_key: CanonicalAtomKey,
    state_canonical_label_index_map: Callable[[str], Mapping[AtomReference, int]],
) -> Tuple[int, ...]:
    lookup = state_canonical_label_index_map(state_name)
    try:
        return tuple(lookup[label] for label in canonical_key)
    except KeyError as exc:  # pragma: no cover - defensive
        missing_label = exc.args[0]
        raise ValueError(
            f"Canonical atom reference {missing_label!r} not present in state '{state_name}'."
        ) from exc


def particle_indices_from_canonical_pair(
    *,
    state_name: str,
    canonical_pair: CanonicalPairKey,
    state_canonical_label_index_map: Callable[[str], Mapping[AtomReference, int]],
) -> Tuple[int, int]:
    lookup = state_canonical_label_index_map(state_name)
    try:
        return lookup[canonical_pair[0]], lookup[canonical_pair[1]]
    except KeyError as exc:  # pragma: no cover - defensive
        missing_label = exc.args[0]
        raise ValueError(
            f"Canonical atom reference {missing_label!r} not present in state '{state_name}'."
        ) from exc


__all__ = [
    "normalize_atom_reference",
    "canonical_pair_key",
    "atom_reference_from_topology_atom",
    "state_atom_references",
    "state_particle_masses",
    "atom_reference_from_index",
    "canonical_atom_label",
    "canonical_bond_key_from_indices",
    "canonical_angle_key_from_indices",
    "canonical_torsion_key_from_indices",
    "state_canonical_label_index_map",
    "particle_indices_from_canonical_key",
    "particle_indices_from_canonical_pair",
]
