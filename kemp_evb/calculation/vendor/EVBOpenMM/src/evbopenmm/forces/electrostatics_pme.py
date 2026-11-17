"""Periodic (PME) electrostatics builders for the EVB scaffold.

EVB electrostatics is a lambda-weighted sum over diabatic states,
``E_elec = sum_s lambda_s * E_elec(q_s)``, where each state ``s`` has its own
atomic charge vector ``q_s``. A single PME pass cannot reproduce this because the
Ewald energy is *quadratic* in charge (``E(sum lambda q) != sum lambda E(q)``).

Stage 1 (exact, this module) builds **one PME-only ``openmm.NonbondedForce`` per
diabatic state**, each carrying that state's complete charge vector. The whole
Ewald energy of state ``s`` is scaled by ``lambda_s`` through a global parameter
``g_s`` driven to ``sqrt(lambda_s)``: charges are offset to ``g_s * q_s`` so every
bilinear Ewald term (real, reciprocal, self, exclusion correction) scales by
``g_s^2 = lambda_s`` exactly. Van der Waals stays in the custom CutoffPeriodic
forces; these forces set ``epsilon = 0`` so they contribute Coulomb only.
"""

from __future__ import annotations

from collections.abc import Sequence

import openmm
from openmm import unit


def build_state_pme_force(
    *,
    state_charges: Sequence[float],
    g_param: str,
    lambda_param: str,
    cutoff_distance: unit.Quantity,
    ewald_error_tolerance: float,
) -> openmm.NonbondedForce:
    """Build a Coulomb-only PME ``NonbondedForce`` for one diabatic state.

    Parameters
    ----------
    state_charges :
        Complete per-particle charge (elementary units) for this diabatic state,
        indexed by scaffold particle index.
    g_param :
        Name of the global parameter (driven to ``sqrt(lambda_state)``) that scales
        every particle charge via :meth:`addParticleParameterOffset`.
    lambda_param :
        Name of this state's lambda parameter. Registered as an (otherwise unused)
        global so it always exists on the context — even when no other force is
        state-specific — letting ``set_state_lambdas`` set it and the ``g_param``
        sync derive ``sqrt(lambda_state)`` from it.
    cutoff_distance :
        Real-space cutoff (the reciprocal sum handles the long range).
    ewald_error_tolerance :
        PME accuracy target.

    Notes
    -----
    The caller is responsible for adding exclusions (as zero-``chargeProd``
    exceptions) for every pair handled elsewhere — intramolecular 1-2/1-3/1-4
    pairs and reactive soft-pair pairs — so the reciprocal exclusion correction
    (which scales as ``g^2`` from the offset charges) lines up with the
    lambda-scaled direct-Coulomb restoration done by the exception bond forces.
    """

    force = openmm.NonbondedForce()
    force.setNonbondedMethod(openmm.NonbondedForce.PME)
    force.setCutoffDistance(cutoff_distance)
    force.setEwaldErrorTolerance(float(ewald_error_tolerance))
    force.setUseDispersionCorrection(False)
    force.addGlobalParameter(g_param, 0.0)
    if lambda_param != g_param:
        # Declared so the context always exposes lambda_state; not referenced by the
        # force itself (the g_param offset carries the actual charge scaling).
        force.addGlobalParameter(lambda_param, 0.0)
    for charge in state_charges:
        # Base charge 0; the live charge is g_param * charge, so the whole Ewald
        # energy scales as g_param^2 = lambda_state.
        index = force.addParticle(0.0, 1.0, 0.0)
        if abs(float(charge)) > 0.0:
            force.addParticleParameterOffset(g_param, index, float(charge), 0.0, 0.0)
    return force


__all__ = ["build_state_pme_force"]
