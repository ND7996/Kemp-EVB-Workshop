"""Audit a built ``openmm.System`` for bonded terms silently dropped by the force field.

OpenMM's bonded generators skip any bond, angle or torsion whose atom-class combination
has no parameter in the force field -- **without raising**. ``createSystem`` returns
successfully and the term simply does not exist. That is a dangerous default for EVB,
where diabats routinely retype atoms (an sp2 carbonyl carbon becoming sp3 in a
tetrahedral intermediate, say) and thereby orphan parameters that existed in the other
state. The resulting system looks healthy, runs, and reports finite energies.

This module compares a System's bonded coverage against its Topology's bond graph and
reports what is missing.

Constraints count as coverage. Under ``constraints=HBonds`` an X-H bond is removed from
the ``HarmonicBondForce`` and becomes an ``openmm`` constraint; under ``HAngles`` (and
``rigidWater``) an angle becomes a 1-3 constraint between its outer atoms. Ignoring
these would report a false positive for every hydrogen in the system.

Example
-------
>>> from evbopenmm.forces import audit_dropped_bonded_terms
>>> report = audit_dropped_bonded_terms(system, topology)
>>> if not report.ok:
...     print(report.describe(topology))
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import openmm
from openmm import unit

__all__ = [
    "bond_graph_from_system",
    "DroppedTermReport",
    "MissingExclusionReport",
    "audit_dropped_bonded_terms",
    "audit_missing_exclusions",
]

BondKey = tuple[int, int]
AngleKey = tuple[int, int, int]


@dataclass(frozen=True)
class DroppedTermReport:
    """Bonded terms present in the topology but absent from the System's forces."""

    dropped_bonds: tuple[BondKey, ...]
    dropped_angles: tuple[AngleKey, ...]
    checked_bonds: int
    checked_angles: int

    @property
    def ok(self) -> bool:
        """True when every audited term has a corresponding force entry."""
        return not self.dropped_bonds and not self.dropped_angles

    def describe(self, topology=None, *, max_items: int = 15) -> str:
        """Human-readable summary. Pass the topology to get atom names instead of indices."""
        if self.ok:
            return (f"No dropped bonded terms ({self.checked_bonds} bonds, "
                    f"{self.checked_angles} angles checked).")

        names: Mapping[int, str] = {}
        if topology is not None:
            names = {a.index: f"{a.residue.name}{a.residue.id}:{a.name}"
                     for a in topology.atoms()}

        def label(i: int) -> str:
            return names.get(i, str(i))

        lines: list[str] = []
        if self.dropped_bonds:
            lines.append(f"{len(self.dropped_bonds)} bond(s) with no force term "
                         f"(of {self.checked_bonds} checked):")
            for i, j in self.dropped_bonds[:max_items]:
                lines.append(f"    {label(i)} -- {label(j)}")
            if len(self.dropped_bonds) > max_items:
                lines.append(f"    ... and {len(self.dropped_bonds) - max_items} more")
        if self.dropped_angles:
            lines.append(f"{len(self.dropped_angles)} angle(s) with no force term "
                         f"(of {self.checked_angles} checked):")
            for i, j, k in self.dropped_angles[:max_items]:
                lines.append(f"    {label(i)} - {label(j)} - {label(k)}")
            if len(self.dropped_angles) > max_items:
                lines.append(f"    ... and {len(self.dropped_angles) - max_items} more")
        return "\n".join(lines)


def _covered_bonds(system: openmm.System) -> set[BondKey]:
    """Atom pairs carrying a two-body bonded term or a distance constraint."""
    covered: set[BondKey] = set()
    for index in range(system.getNumConstraints()):
        i, j, _ = system.getConstraintParameters(index)
        covered.add((min(i, j), max(i, j)))
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            for index in range(force.getNumBonds()):
                i, j, *_ = force.getBondParameters(index)
                covered.add((min(i, j), max(i, j)))
        elif isinstance(force, openmm.CustomBondForce):
            # EVB per-state bonded terms, Morse replacements and soft-pairs all live here.
            for index in range(force.getNumBonds()):
                i, j, *_ = force.getBondParameters(index)
                covered.add((min(i, j), max(i, j)))
    return covered




def _looks_like_evb_scaffold(system: openmm.System) -> bool:
    """True when a System carries ``lambda_*`` globals, i.e. it is an EVB scaffold.

    Per-state diabat Systems never have them, so this cleanly separates the two.
    """
    for force in system.getForces():
        if not hasattr(force, "getNumGlobalParameters"):
            continue
        for index in range(force.getNumGlobalParameters()):
            if force.getGlobalParameterName(index).startswith("lambda_"):
                return True
    return False


_SCAFFOLD_MISUSE = (
    "This looks like an EVB scaffold (it carries lambda_* global parameters), not a "
    "per-state diabat System. These audits describe ONE state's force field: on a "
    "scaffold the bonded terms are split between shared and per-state forces and the "
    "exclusions are spread across several nonbonded forces, so the result is not a "
    "statement about any diabat and can come back clean while a state is broken. Pass "
    "each per-state System instead, or use EVBSystem.validate_diabat_topology()."
)


def bond_graph_from_system(system: openmm.System) -> set[BondKey]:
    """The bonds a built ``System`` actually contains.

    For a hand-assembled diabatic state the System -- not the Topology -- is the
    authority on connectivity. A Topology read from a PDB without CONECT records gets
    standard residue-template bonds, so it still carries bonds the diabat has broken and
    knows nothing of the ones it has formed. Auditing against it therefore reports the
    broken bond as a defect while missing the formed one entirely.

    Includes ``HarmonicBondForce`` entries and distance constraints -- the latter because
    ``constraints=HBonds`` moves X-H bonds out of the bond force.

    ``CustomBondForce`` is deliberately **excluded**. It carries two-body terms that are
    not bonds: the Q6 soft-pair acts on a reactive pair precisely in the state where that
    pair is *not* bonded, so counting it would resurrect a broken bond and invent 1-3
    neighbours from it. Call this on a per-state ``System``, not on an EVB scaffold.

    .. note::
       Under ``constraints=HAngles`` some constraints span 1-3 pairs rather than bonds and
       would be misread here as bonds. Audit such systems against an authoritative
       Topology instead.
    """
    if _looks_like_evb_scaffold(system):
        warnings.warn(_SCAFFOLD_MISUSE, stacklevel=2)
    pairs: set[BondKey] = set()
    for index in range(system.getNumConstraints()):
        i, j, _ = system.getConstraintParameters(index)
        pairs.add((min(i, j), max(i, j)))
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            for index in range(force.getNumBonds()):
                i, j, *_ = force.getBondParameters(index)
                pairs.add((min(i, j), max(i, j)))
    return pairs


def _covered_angles(system: openmm.System) -> set[AngleKey]:
    """Angles carrying a three-body term.

    1-3 constraints (HAngles / rigidWater) are handled by the caller, which has the
    constrained-pair set to hand.
    """
    covered: set[AngleKey] = set()
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicAngleForce):
            for index in range(force.getNumAngles()):
                i, j, k, *_ = force.getAngleParameters(index)
                covered.add((min(i, k), j, max(i, k)))
        elif isinstance(force, openmm.CustomAngleForce):
            for index in range(force.getNumAngles()):
                i, j, k, *_ = force.getAngleParameters(index)
                covered.add((min(i, k), j, max(i, k)))
    return covered


def audit_dropped_bonded_terms(
    system: openmm.System,
    topology,
    *,
    check_angles: bool = True,
    check_bonds: bool = True,
    bonds: Iterable[BondKey] | None = None,
    ignore_residue_names: Iterable[str] | None = None,
    raise_on_error: bool = False,
    warn: bool = True,
) -> DroppedTermReport:
    """Report bonded terms in ``topology`` that have no counterpart in ``system``.

    Parameters
    ----------
    system :
        A built ``openmm.System``.
    topology :
        The ``openmm.app.Topology`` it was built from. Its bond graph is the reference:
        every bond, and every angle implied by that graph, must be covered.
    check_angles :
        Audit angles as well as bonds. Angles are derived from the bond graph, so this
        is the check that catches a retyped centre losing its angle parameters.
    check_bonds :
        Audit bonds. Set ``False`` when ``bonds`` was derived from this same System
        (e.g. via :func:`bond_graph_from_system`): reference and coverage are then the
        same set and the check cannot fail.
    bonds :
        Reference connectivity. Defaults to the topology's bond graph. Pass a
        System-derived graph when the Topology is not authoritative for this state --
        a PDB-derived Topology carries standard template bonds, not a diabat's.
    ignore_residue_names :
        Residue names to skip entirely (e.g. ``{"HOH", "WAT"}`` for a solvated system,
        where rigid water is expected to carry constraints rather than terms). Any term
        with **all** atoms inside an ignored residue is skipped.
    raise_on_error :
        Raise ``RuntimeError`` instead of warning when terms are missing.
    warn :
        Emit a ``UserWarning`` when terms are missing. Ignored if ``raise_on_error``.

    Returns
    -------
    DroppedTermReport

    Notes
    -----
    Constraints count as coverage: ``constraints=HBonds`` moves X-H bonds out of the
    ``HarmonicBondForce``, and ``HAngles`` / ``rigidWater`` turn angles into 1-3
    constraints between the outer atoms. Counting only force entries would flag every
    such term as dropped.
    """
    if warn and _looks_like_evb_scaffold(system):
        warnings.warn(_SCAFFOLD_MISUSE, stacklevel=2)

    ignored = {name.upper() for name in (ignore_residue_names or ())}
    skip_atoms: set[int] = set()
    if ignored:
        for residue in topology.residues():
            if residue.name.upper() in ignored:
                skip_atoms.update(atom.index for atom in residue.atoms())

    reference_bonds = (
        [(min(i, j), max(i, j)) for i, j in bonds] if bonds is not None
        else [(min(a.index, b.index), max(a.index, b.index)) for a, b in topology.bonds()]
    )
    adjacency: dict[int, set[int]] = {}
    audited_bonds: list[BondKey] = []
    for i, j in reference_bonds:
        adjacency.setdefault(i, set()).add(j)
        adjacency.setdefault(j, set()).add(i)
        if i in skip_atoms and j in skip_atoms:
            continue
        audited_bonds.append((i, j))

    if check_bonds:
        covered_bonds = _covered_bonds(system)
        dropped_bonds = tuple(sorted(set(audited_bonds) - covered_bonds))
    else:
        # Tautological when `bonds` was derived from this same System: the reference and
        # the coverage would be the same set, so report nothing rather than a count that
        # implies bonds were verified.
        dropped_bonds = ()
        audited_bonds = []

    dropped_angles: tuple[AngleKey, ...] = ()
    checked_angles = 0
    if check_angles:
        constrained_pairs: set[BondKey] = set()
        for index in range(system.getNumConstraints()):
            i, j, _ = system.getConstraintParameters(index)
            constrained_pairs.add((min(i, j), max(i, j)))

        expected: set[AngleKey] = set()
        for centre, neighbours in adjacency.items():
            ordered = sorted(neighbours)
            for a in range(len(ordered)):
                for b in range(a + 1, len(ordered)):
                    i, k = ordered[a], ordered[b]
                    if {i, centre, k} <= skip_atoms:
                        continue
                    expected.add((i, centre, k))
        checked_angles = len(expected)
        covered_angles = _covered_angles(system)
        missing = {
            angle for angle in expected - covered_angles
            # HAngles / rigidWater replace an angle with a 1-3 constraint.
            if (angle[0], angle[2]) not in constrained_pairs
        }
        dropped_angles = tuple(sorted(missing))

    report = DroppedTermReport(
        dropped_bonds=dropped_bonds,
        dropped_angles=dropped_angles,
        checked_bonds=len(audited_bonds),
        checked_angles=checked_angles,
    )

    if not report.ok:
        message = (
            "Force field silently dropped bonded terms: OpenMM's bonded generators skip "
            "atom-class combinations with no matching parameter instead of raising, so "
            "createSystem succeeded but these terms do not exist in the System. This is "
            "common after retyping atoms for a reactive state (e.g. sp2 -> sp3 at a "
            "forming bond), which orphans parameters that existed in the other state.\n"
            + report.describe(topology)
        )
        if raise_on_error:
            raise RuntimeError(message)
        if warn:
            warnings.warn(message, stacklevel=2)
    return report


@dataclass(frozen=True)
class MissingExclusionReport:
    """1-2 / 1-3 pairs implied by the bond graph that carry no nonbonded exclusion."""

    missing_one_two: tuple[BondKey, ...]
    missing_one_three: tuple[BondKey, ...]
    checked_one_two: int
    checked_one_three: int

    @property
    def ok(self) -> bool:
        return not self.missing_one_two and not self.missing_one_three

    def describe(self, topology=None, *, max_items: int = 15) -> str:
        if self.ok:
            return (f"No missing exclusions ({self.checked_one_two} 1-2, "
                    f"{self.checked_one_three} 1-3 pairs checked).")
        names: Mapping[int, str] = {}
        if topology is not None:
            names = {a.index: f"{a.residue.name}{a.residue.id}:{a.name}"
                     for a in topology.atoms()}

        def label(i: int) -> str:
            return names.get(i, str(i))

        lines: list[str] = []
        for title, pairs, total in (
            ("1-2", self.missing_one_two, self.checked_one_two),
            ("1-3", self.missing_one_three, self.checked_one_three),
        ):
            if not pairs:
                continue
            lines.append(f"{len(pairs)} {title} pair(s) with no exclusion (of {total}):")
            for i, j in pairs[:max_items]:
                lines.append(f"    {label(i)} ... {label(j)}")
            if len(pairs) > max_items:
                lines.append(f"    ... and {len(pairs) - max_items} more")
        return "\n".join(lines)


def _excluded_pairs(system: openmm.System) -> set[BondKey]:
    """Pairs whose nonbonded interaction is fully switched off in this System.

    A NonbondedForce exception counts only when BOTH chargeProd and epsilon are zero:
    scaled 1-4 exceptions are not exclusions. CustomNonbondedForce exclusions count
    directly.
    """
    excluded: set[BondKey] = set()
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for index in range(force.getNumExceptions()):
                i, j, charge_prod, _sigma, epsilon = force.getExceptionParameters(index)
                q = abs(charge_prod.value_in_unit(unit.elementary_charge ** 2))
                e = abs(epsilon.value_in_unit(unit.kilojoule_per_mole))
                if q < 1e-12 and e < 1e-12:
                    excluded.add((min(i, j), max(i, j)))
        elif isinstance(force, openmm.CustomNonbondedForce):
            for index in range(force.getNumExclusions()):
                i, j = force.getExclusionParticles(index)
                excluded.add((min(i, j), max(i, j)))
    return excluded


def audit_missing_exclusions(
    system: openmm.System,
    topology,
    *,
    bonds: Iterable[BondKey] | None = None,
    ignore_residue_names: Iterable[str] | None = None,
    raise_on_error: bool = False,
    warn: bool = True,
) -> MissingExclusionReport:
    """Report 1-2 / 1-3 pairs from the bond graph that are not excluded from nonbonded.

    A correct force field excludes every 1-2 and 1-3 pair from nonbonded interactions
    (1-4 pairs are *scaled*, not excluded, and are therefore not audited here). When a
    System is assembled by hand -- for example by editing a prmtop's bond list to create
    a diabatic state -- the exclusion list is easy to leave stale: the new bond exists,
    but the 1-3 neighbours it creates are still counted at full Coulomb.

    The consequence is severe and silent. In a product diabat carrying large formal
    charges, an un-excluded 1-3 pair between two formally negative atoms is a full
    Coulomb repulsion at ~3 A, which can drive the reacting fragments apart so the
    reaction is never sampled -- while the energy gap, and hence the free energy, still
    looks perfectly finite.

    Parameters
    ----------
    system, topology :
        The built System and the Topology whose bond graph defines the expected
        exclusions.
    ignore_residue_names :
        Residue names to skip (pairs with all atoms inside an ignored residue).
    raise_on_error :
        Raise ``RuntimeError`` instead of warning.
    warn :
        Emit a ``UserWarning`` when exclusions are missing. Ignored if ``raise_on_error``.

    Returns
    -------
    MissingExclusionReport
    """
    if warn and _looks_like_evb_scaffold(system):
        warnings.warn(_SCAFFOLD_MISUSE, stacklevel=2)

    ignored = {name.upper() for name in (ignore_residue_names or ())}
    skip_atoms: set[int] = set()
    if ignored:
        for residue in topology.residues():
            if residue.name.upper() in ignored:
                skip_atoms.update(atom.index for atom in residue.atoms())

    reference_bonds = (
        [(i, j) for i, j in bonds] if bonds is not None
        else [(a.index, b.index) for a, b in topology.bonds()]
    )
    adjacency: dict[int, set[int]] = {}
    for i, j in reference_bonds:
        adjacency.setdefault(i, set()).add(j)
        adjacency.setdefault(j, set()).add(i)

    one_two: set[BondKey] = set()
    for i, neighbours in adjacency.items():
        for j in neighbours:
            if i in skip_atoms and j in skip_atoms:
                continue
            one_two.add((min(i, j), max(i, j)))

    one_three: set[BondKey] = set()
    for centre, neighbours in adjacency.items():
        ordered = sorted(neighbours)
        for a in range(len(ordered)):
            for b in range(a + 1, len(ordered)):
                i, k = ordered[a], ordered[b]
                if k in adjacency.get(i, ()):  # also 1-2: already covered above
                    continue
                if i in skip_atoms and k in skip_atoms and centre in skip_atoms:
                    continue
                one_three.add((i, k))

    excluded = _excluded_pairs(system)
    report = MissingExclusionReport(
        missing_one_two=tuple(sorted(one_two - excluded)),
        missing_one_three=tuple(sorted(one_three - excluded)),
        checked_one_two=len(one_two),
        checked_one_three=len(one_three),
    )

    if not report.ok:
        message = (
            "Nonbonded exclusions are inconsistent with the bond graph: pairs that are "
            "1-2 or 1-3 in this topology are not excluded from nonbonded interactions. "
            "This is the signature of a hand-assembled diabatic state whose bond list "
            "was edited without regenerating the exclusion list, and it counts "
            "forming-bond neighbours at full Coulomb.\n" + report.describe(topology)
        )
        if raise_on_error:
            raise RuntimeError(message)
        if warn:
            warnings.warn(message, stacklevel=2)
    return report
