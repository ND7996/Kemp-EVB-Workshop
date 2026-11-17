from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Tuple

from openmm import app

from .types import AtomReference, CanonicalAtomKey, CanonicalPairKey


@dataclass(frozen=True)
class BondedTermEntry:
    """Holds bonded parameters extracted from an OpenMM force."""

    parameters: Tuple
    original_indices: Tuple[int, ...]


@dataclass(frozen=True)
class BondedForceComparison:
    """Summary of bonded-term differences for a single OpenMM force type."""

    force_name: str
    universal_identical: Mapping[CanonicalAtomKey, BondedTermEntry]
    universal_conflicting: Mapping[CanonicalAtomKey, Mapping[str, BondedTermEntry]]
    partial: Mapping[CanonicalAtomKey, Mapping[str, BondedTermEntry]]


@dataclass(frozen=True)
class VanDerWaalsPartition:
    """Bookkeeping container describing the EVB VdW partition."""

    universal_atoms: Mapping[AtomReference, Tuple[float, float]]
    state_specific_atoms: Mapping[str, Mapping[AtomReference, Tuple[float, float]]]
    universal_exceptions: Mapping[CanonicalPairKey, Tuple[float, float]]
    state_specific_exceptions: Mapping[str, Mapping[CanonicalPairKey, Tuple[float, float]]]
    universal_zero_exclusions: Tuple[CanonicalPairKey, ...]
    state_specific_zero_exclusions: Mapping[str, Tuple[CanonicalPairKey, ...]]


@dataclass(frozen=True)
class ElectrostaticsPartition:
    """Bookkeeping container describing the EVB electrostatics partition."""

    universal_atoms: Mapping[AtomReference, float]
    state_specific_atoms: Mapping[str, Mapping[AtomReference, float]]
    universal_exceptions: Mapping[CanonicalPairKey, float]
    state_specific_exceptions: Mapping[str, Mapping[CanonicalPairKey, float]]
    universal_zero_exclusions: Tuple[CanonicalPairKey, ...]
    state_specific_zero_exclusions: Mapping[str, Tuple[CanonicalPairKey, ...]]


@dataclass(frozen=True)
class ResidueTarget:
    """Reference to a residue occurrence within one EVB state."""

    state_name: str
    modeller: app.Modeller
    residue: app.topology.Residue


@dataclass(frozen=True)
class ResidueTargetIndex:
    """Prepared topology views used to resolve residues consistently across states."""

    modellers: Mapping[str, app.Modeller]
    topology_bonds: Mapping[str, Tuple[Tuple[app.Topology.Atom, app.Topology.Atom], ...]]
    residue_targets: Mapping[str, list[ResidueTarget]]


@dataclass(frozen=True)
class FEPSampleRecord:
    """Single sampled configuration used for FEP/umbrella post-processing."""

    md_step: int
    effective_energy_kj_mol: float
    reaction_coordinate_kj_mol: float
    state_energies_kj_mol: Tuple[float, ...]
    # Geometric distance (nm) of the reactive pair that defines the bond being
    # made/broken, captured per frame when a distance-dependent off-diagonal
    # coupling H12(r)=A*exp(-mu*(r-r0)) is requested. ``None`` for the default
    # (constant / gap-Gaussian) workflows, which do not need a geometry.
    coupling_distance_nm: float | None = None


@dataclass(frozen=True)
class MDEnergySampleRecord:
    """Single MD sample containing instantaneous energies."""

    md_step: int
    potential_energy_kj_mol: float
    kinetic_energy_kj_mol: float
    total_energy_kj_mol: float


@dataclass(frozen=True)
class FEPWindowRecord:
    """Container storing all sampled frames for one lambda window."""

    index: int
    lambda_values: Tuple[float, ...]
    samples: Tuple[FEPSampleRecord, ...]


@dataclass(frozen=True)
class MorseBondParameter:
    """Morse-potential parameters for a forming/breaking reactive bond.

    The parameters describe the dissociation curve of the bond as

        ``U(r) = de * (1 - exp(-alpha * (r - r0)))^2``

    in OpenMM's native units (kJ/mol, nm, 1/nm). EVBOpenMM **requires** the
    user to provide one of these for every reactive (forming/breaking)
    pair: the Morse parameters cannot be derived from a force-field
    parameterisation step (e.g. OpenFF / antechamber) because a
    force-field bond potential is harmonic by construction and does not
    capture dissociation. The expectation is that callers fit Morse
    parameters to a QM dissociation curve (or import from a Q6
    ``.fep`` file).
    """

    r0_nm: float
    de_kj_mol: float
    alpha_per_nm: float


__all__ = [
    "BondedTermEntry",
    "BondedForceComparison",
    "VanDerWaalsPartition",
    "ElectrostaticsPartition",
    "ResidueTarget",
    "ResidueTargetIndex",
    "FEPSampleRecord",
    "MDEnergySampleRecord",
    "FEPWindowRecord",
    "MorseBondParameter",
]
