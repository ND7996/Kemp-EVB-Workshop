from __future__ import annotations

from collections.abc import Mapping
from typing import Dict, Tuple

from openmm import app

from ..core import normalize_residue_name
from ..input_reader.system_reader import EVBState
from ..models import ResidueTarget, ResidueTargetIndex


def compute_missing_force_field_residues(
    states: Mapping[str, EVBState],
    normalized_force_field_residues: set[str],
) -> Dict[str, Tuple[str, ...]]:
    missing_by_state: Dict[str, Tuple[str, ...]] = {}
    for name, state in states.items():
        missing = tuple(
            sorted(
                residue
                for residue in state.residue_names
                if normalize_residue_name(residue) not in normalized_force_field_residues
            )
        )
        if missing:
            missing_by_state[name] = missing
    return missing_by_state


def build_residue_targets(states: Mapping[str, EVBState]) -> ResidueTargetIndex:
    modellers: Dict[str, app.Modeller] = {}
    topology_bonds: Dict[str, Tuple[Tuple[app.Topology.Atom, app.Topology.Atom], ...]] = {}
    residue_targets: Dict[str, list[ResidueTarget]] = {}

    for state_name, state in states.items():
        modeller = app.Modeller(state.topology, state.positions)
        modellers[state_name] = modeller
        topology_bonds[state_name] = tuple(tuple(bond) for bond in modeller.topology.bonds())
        for residue in modeller.topology.residues():
            residue_name = (residue.name or "").strip()
            if not residue_name:
                continue
            normalized = normalize_residue_name(residue_name)
            target = ResidueTarget(state_name=state_name, modeller=modeller, residue=residue)
            residue_targets.setdefault(normalized, []).append(target)

    return ResidueTargetIndex(
        modellers=modellers,
        topology_bonds=topology_bonds,
        residue_targets=residue_targets,
    )


__all__ = [
    "compute_missing_force_field_residues",
    "build_residue_targets",
]
