"""OpenMM-style reporters that write EVB-specific traces during MD.

Each reporter follows the OpenMM Reporter protocol
(``describeNextReport``/``report`` on a Simulation), so they attach via
``simulation.reporters.append(reporter)``.

Modules
-------
adiabatic   – trace of the direct adiabatic-EVB integrator's mixing scalars
              (E_U, E_S1, E_S2, w1, w2, H12, umbrella term).
energy_gap  – ⟨ΔU⟩ sampling on a *mapping-λ* Simulation: at each
              interval switches through two probe λ combinations, reads
              the potential (optionally per force group), restores the
              mapping λ, and writes one CSV row.

SCAAS-specific reporters (SCAASPolarizationReporter etc.) intentionally
stay in ``evbopenmm.scaas`` — SCAAS is its own module orthogonal to EVB.
"""
from __future__ import annotations

from .adiabatic import AdiabaticReporter
from .energy_gap import (
    EnergyGapReporter,
    block_averaged_sem,
    load_gap_samples,
    running_mean,
)

__all__ = [
    "AdiabaticReporter",
    "EnergyGapReporter",
    "load_gap_samples",
    "running_mean",
    "block_averaged_sem",
]
