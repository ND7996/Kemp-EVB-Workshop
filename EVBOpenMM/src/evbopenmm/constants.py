from __future__ import annotations

from typing import Tuple

import openmm

STATE_ENERGY_TERMS: Tuple[str, ...] = ("nonbonded", "bond", "angle", "torsion")
ONE_4PI_EPS0 = 138.935456  # kJ/mol * nm / e^2
KJ_PER_KCAL = 4.184
REPORTING_ENERGY_UNIT = "kcal/mol"
SOLUTE_RESTRAINT_GLOBAL_PARAMETER = "k_solute_posres"
REFERENCE_FORCE_CATEGORIES: Tuple[Tuple[type, str], ...] = (
    (openmm.HarmonicBondForce, "bond"),
    (openmm.CustomBondForce, "bond"),
    (openmm.HarmonicAngleForce, "angle"),
    (openmm.CustomAngleForce, "angle"),
    (openmm.PeriodicTorsionForce, "torsion"),
    (openmm.CustomTorsionForce, "torsion"),
    (openmm.CMAPTorsionForce, "torsion"),
    (openmm.NonbondedForce, "nonbonded"),
)

__all__ = [
    "STATE_ENERGY_TERMS",
    "ONE_4PI_EPS0",
    "KJ_PER_KCAL",
    "REPORTING_ENERGY_UNIT",
    "SOLUTE_RESTRAINT_GLOBAL_PARAMETER",
    "REFERENCE_FORCE_CATEGORIES",
]
