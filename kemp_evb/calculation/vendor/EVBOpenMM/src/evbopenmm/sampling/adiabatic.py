"""Direct adiabatic ground-state EVB dynamics.

This module provides the runtime building blocks for propagating molecular
dynamics directly on the adiabatic EVB ground-state surface
:math:`E_-(x)`, the lower eigenvalue of the 2×2 EVB Hamiltonian.

The implementation rests on a Hellmann–Feynman identity: for both constant
and Gaussian-in-gap couplings, the gradient of :math:`E_-` is a linear
combination of the diabatic gradients with weights :math:`w_1, w_2` that
sum to one. The mixing scalars are computed each step from the per-state
diabatic energies, and OpenMM's :class:`openmm.CustomIntegrator` evaluates
the per-force-group forces in a single pass.

OpenMM's CustomIntegrator forbids any single computation step from
referencing more than one force group. So we accumulate per-group energies
and per-group forces into intermediate variables (``E_U``, ``E_S1``,
``E_S2``, ``F_adia``) one group at a time, then perform the mixing using
only those scalar/vector intermediates.

See :doc:`/docs/evb_adiabatic_dynamics_plan` for the math and the wider
plan.
"""

from __future__ import annotations

import math
from typing import Dict, Sequence, Tuple

import openmm
from openmm import unit


def compute_adiabatic_mixing(
    *,
    state1_energy_kj_mol: float,
    state2_energy_kj_mol: float,
    coupling_amplitude_kj_mol: float,
    coupling_width_per_kj2_mol2: float = 0.0,
    diabatic_offset_kj_mol: float = 0.0,
    umbrella_force_constant_per_kj_mol: float = 0.0,
    umbrella_dv_center_kj_mol: float = 0.0,
) -> Dict[str, float]:
    """Pure-Python reference for the per-step adiabatic mixing scalars.

    The keys (``dV``, ``H12``, ``denom``, ``c1sq``, ``c2sq``, ``beta``,
    ``w1``, ``w2``, ``w1_eff``, ``w2_eff``, ``U_umbrella``, ``E_minus``)
    mirror the global variables held inside the CustomIntegrator one-to-one,
    so this function is the canonical reference for the symbolic expressions
    baked into the integrator.

    Conventions follow :func:`evbopenmm.analysis.evb_ground_state_energy_kj_mol`:
    ``diabatic_offset_kj_mol`` (``α``) is folded onto state 2 only, i.e.
    :math:`H_{11}=V_1`, :math:`H_{22}=V_2 + α`. Setting
    ``coupling_width_per_kj2_mol2 = 0`` recovers constant coupling
    :math:`H_{12} = A`.

    The optional umbrella restraint
    ``U_umbrella = 0.5 * K * (ΔV − ΔV_center)²`` shifts the per-state
    weights to ``w1_eff = w1 − K(ΔV − ΔV_center)`` and
    ``w2_eff = w2 + K(ΔV − ΔV_center)``. With ``K = 0`` (the default)
    ``w*_eff == w*``.
    """

    V1 = float(state1_energy_kj_mol)
    V2 = float(state2_energy_kj_mol)
    A = float(coupling_amplitude_kj_mol)
    mu = float(coupling_width_per_kj2_mol2)
    alpha = float(diabatic_offset_kj_mol)
    K = float(umbrella_force_constant_per_kj_mol)
    dV_center = float(umbrella_dv_center_kj_mol)

    dV = V2 - V1 + alpha
    H12 = A * math.exp(-mu * dV * dV)
    denom_sq = dV * dV + 4.0 * H12 * H12
    denom = math.sqrt(denom_sq)
    if denom == 0.0:
        # Fully degenerate (dV = 0 and A = 0) — both eigenvalues equal V1.
        # Any unit-norm superposition is admissible; pick the symmetric mix.
        c1sq = 0.5
        c2sq = 0.5
        beta = 0.0
    else:
        c1sq = 0.5 * (1.0 + dV / denom)
        c2sq = 0.5 * (1.0 - dV / denom)
        beta = 4.0 * mu * dV * H12 * H12 / denom
    w1 = c1sq - beta
    w2 = c2sq + beta
    umbrella_shift = K * (dV - dV_center)
    w1_eff = w1 - umbrella_shift
    w2_eff = w2 + umbrella_shift
    u_umbrella = 0.5 * K * (dV - dV_center) * (dV - dV_center)
    e_minus = 0.5 * (V1 + V2 + alpha - denom)
    return {
        "dV": dV,
        "H12": H12,
        "denom": denom,
        "c1sq": c1sq,
        "c2sq": c2sq,
        "beta": beta,
        "w1": w1,
        "w2": w2,
        "w1_eff": w1_eff,
        "w2_eff": w2_eff,
        "U_umbrella": u_umbrella,
        "E_minus": e_minus,
    }


def _mask_to_term_indices(mask: int) -> Tuple[int, ...]:
    if int(mask) < 0:
        raise ValueError(f"Force-group mask must be non-negative; got {mask!r}.")
    return tuple(idx for idx in range(32) if mask & (1 << idx))


_INTEGRATOR_GLOBALS: Tuple[str, ...] = (
    "dV",
    "H12",
    "denom",
    "c1sq",
    "c2sq",
    "beta",
    "w1",
    "w2",
    "w1_eff",
    "w2_eff",
    "E_U",
    "E_S1",
    "E_S2",
    "U_umbrella",
)


def _add_mixing_globals(integrator: openmm.CustomIntegrator) -> None:
    for name in _INTEGRATOR_GLOBALS:
        integrator.addGlobalVariable(name, 0.0)


def _accumulate_group_energy(
    integrator: openmm.CustomIntegrator,
    *,
    target: str,
    indices: Sequence[int],
) -> None:
    """Accumulate ``target = sum(energy{i} for i in indices)``.

    Each ``addComputeGlobal`` step references at most one force group, so the
    per-group sum is built incrementally with one step per group.
    """

    if not indices:
        integrator.addComputeGlobal(target, "0")
        return
    integrator.addComputeGlobal(target, f"energy{indices[0]}")
    for idx in indices[1:]:
        integrator.addComputeGlobal(target, f"{target} + energy{idx}")


def _accumulate_group_force(
    integrator: openmm.CustomIntegrator,
    *,
    target: str,
    contributions: Sequence[Tuple[int, str | None]],
    reset_first: bool,
) -> None:
    """Accumulate ``target = sum(weight*f{i} for (i, weight) in contributions)``.

    ``weight=None`` means an implicit weight of 1. Each step references at
    most one force group.
    """

    if reset_first:
        integrator.addComputePerDof(target, "0")
    for index, weight_expr in contributions:
        if weight_expr is None:
            integrator.addComputePerDof(target, f"{target} + f{index}")
        else:
            integrator.addComputePerDof(target, f"{target} + ({weight_expr})*f{index}")


def _add_diabatic_energies_step(
    integrator: openmm.CustomIntegrator,
    *,
    u_mask: int,
    s1_mask: int,
    s2_mask: int,
) -> None:
    """Compute E_U, E_S1, E_S2 in three independent accumulator phases."""

    _accumulate_group_energy(
        integrator, target="E_U", indices=_mask_to_term_indices(u_mask)
    )
    _accumulate_group_energy(
        integrator, target="E_S1", indices=_mask_to_term_indices(s1_mask)
    )
    _accumulate_group_energy(
        integrator, target="E_S2", indices=_mask_to_term_indices(s2_mask)
    )


def _add_mixing_step(integrator: openmm.CustomIntegrator) -> None:
    """Compute mixing scalars from already-accumulated E_S1 / E_S2.

    Includes the optional umbrella-on-ΔV restraint
    ``U_umbrella = 0.5 * K_umbrella * (dV - dV_center)^2``. The restraint
    contributes a force-shift term ``-K_umbrella * (dV - dV_center)`` to
    ``w_1`` and ``+K_umbrella * (dV - dV_center)`` to ``w_2``, so the
    effective mixing weights are ``w1_eff`` and ``w2_eff``. With
    ``K_umbrella = 0`` (the default) ``w*_eff == w*`` and the integrator
    is identical to the unrestrained adiabatic dynamics from A2/A3.
    """

    integrator.addComputeGlobal("dV", "E_S2 - E_S1 + alpha_kj")
    integrator.addComputeGlobal("H12", "A_kj * exp(-mu_per_kj2 * dV*dV)")
    integrator.addComputeGlobal("denom", "sqrt(dV*dV + 4*H12*H12)")
    integrator.addComputeGlobal("c1sq", "0.5 * (1 + dV/denom)")
    integrator.addComputeGlobal("c2sq", "0.5 * (1 - dV/denom)")
    integrator.addComputeGlobal("beta", "4 * mu_per_kj2 * dV * H12 * H12 / denom")
    integrator.addComputeGlobal("w1", "c1sq - beta")
    integrator.addComputeGlobal("w2", "c2sq + beta")
    # Umbrella-on-ΔV (Hellmann–Feynman with parametric restraint):
    # U_umb = 0.5 * K * (dV - dV_0)^2
    # ∇ U_umb = K * (dV - dV_0) * ∇ dV = K * (dV - dV_0) * (∇V_2 - ∇V_1)
    # Add the corresponding shifts to the per-state weights.
    integrator.addComputeGlobal(
        "U_umbrella", "0.5 * K_umbrella * (dV - dV_center) * (dV - dV_center)"
    )
    integrator.addComputeGlobal("w1_eff", "w1 - K_umbrella * (dV - dV_center)")
    integrator.addComputeGlobal("w2_eff", "w2 + K_umbrella * (dV - dV_center)")


def _add_adiabatic_force_step(
    integrator: openmm.CustomIntegrator,
    *,
    u_mask: int,
    s1_mask: int,
    s2_mask: int,
) -> None:
    """Compute ``F_adia = sum(fU) + w1_eff*sum(fS1) + w2_eff*sum(fS2)``.

    The effective weights ``w1_eff = w1 - K*(dV - dV_center)`` and
    ``w2_eff = w2 + K*(dV - dV_center)`` fold in the optional umbrella
    restraint on ΔV. With ``K_umbrella = 0`` the effective weights are
    identical to the bare ``w_i`` from the closed-form 2x2 mixing.
    """

    contributions: list[Tuple[int, str | None]] = []
    for idx in _mask_to_term_indices(u_mask):
        contributions.append((idx, None))
    for idx in _mask_to_term_indices(s1_mask):
        contributions.append((idx, "w1_eff"))
    for idx in _mask_to_term_indices(s2_mask):
        contributions.append((idx, "w2_eff"))
    _accumulate_group_force(
        integrator,
        target="F_adia",
        contributions=contributions,
        reset_first=True,
    )


def _validate_force_group_masks(force_group_masks: Tuple[int, int, int]) -> Tuple[int, int, int]:
    if len(force_group_masks) != 3:
        raise ValueError(
            "force_group_masks must be a 3-tuple (universal, state1, state2)."
        )
    u_mask, s1_mask, s2_mask = (int(value) for value in force_group_masks)
    if u_mask & s1_mask or u_mask & s2_mask or s1_mask & s2_mask:
        raise ValueError(
            "force_group_masks overlap: universal/state-1/state-2 must be pairwise disjoint."
        )
    return u_mask, s1_mask, s2_mask


def _resolve_timestep(timestep: unit.Quantity | float) -> unit.Quantity:
    if isinstance(timestep, unit.Quantity):
        return timestep
    return float(timestep) * unit.femtosecond


def _add_coupling_globals(
    integrator: openmm.CustomIntegrator,
    *,
    coupling_amplitude_kj_mol: float,
    coupling_width_per_kj2_mol2: float,
    diabatic_offset_kj_mol: float,
    umbrella_force_constant_per_kj_mol: float = 0.0,
    umbrella_dv_center_kj_mol: float = 0.0,
) -> None:
    integrator.addGlobalVariable("A_kj", float(coupling_amplitude_kj_mol))
    integrator.addGlobalVariable("mu_per_kj2", float(coupling_width_per_kj2_mol2))
    integrator.addGlobalVariable("alpha_kj", float(diabatic_offset_kj_mol))
    integrator.addGlobalVariable("K_umbrella", float(umbrella_force_constant_per_kj_mol))
    integrator.addGlobalVariable("dV_center", float(umbrella_dv_center_kj_mol))


def build_nve_adiabatic_integrator(
    *,
    timestep: unit.Quantity | float,
    force_group_masks: Tuple[int, int, int],
    coupling_amplitude_kj_mol: float,
    coupling_width_per_kj2_mol2: float = 0.0,
    diabatic_offset_kj_mol: float = 0.0,
    umbrella_force_constant_per_kj_mol: float = 0.0,
    umbrella_dv_center_kj_mol: float = 0.0,
) -> openmm.CustomIntegrator:
    """Build a velocity-Verlet CustomIntegrator on the EVB adiabatic ground state.

    Parameters
    ----------
    timestep
        Integration timestep; either an OpenMM ``Quantity`` or a value in
        femtoseconds. Default for adiabatic mode is **0.5 fs** (caller's
        responsibility); larger steps validated empirically in A4.
    force_group_masks
        ``(universal, state1, state2)`` bitmasks from
        :meth:`evbopenmm.EVBSystem._adiabatic_force_group_masks`.
    coupling_amplitude_kj_mol
        ``A`` in :math:`H_{12}(\\Delta V) = A\\,e^{-\\mu \\Delta V^2}`, in
        kJ/mol. For constant coupling pass ``coupling_width_per_kj2_mol2=0``
        and ``A`` equal to ``H_{12}``.
    coupling_width_per_kj2_mol2
        ``μ`` in the Gaussian-in-gap form, in (kJ/mol)⁻². Defaults to 0,
        which recovers constant coupling.
    diabatic_offset_kj_mol
        Diabatic offset ``α = Δ_2 − Δ_1`` folded onto state 2 only. Convention
        matches :func:`evbopenmm.analysis.evb_ground_state_energy_kj_mol`.
    umbrella_force_constant_per_kj_mol
        Force constant of the optional harmonic umbrella restraint on
        ``ΔV``: ``U_umbrella = 0.5 * K * (ΔV − ΔV_center)²``. Units are
        (kJ/mol)⁻¹ — the restraint energy is in kJ/mol when ΔV is in
        kJ/mol and ``K`` is in (kJ/mol)⁻¹. Defaults to 0 (no restraint).
    umbrella_dv_center_kj_mol
        Center of the umbrella restraint, in kJ/mol. Ignored when
        ``umbrella_force_constant_per_kj_mol = 0``.
    """

    u_mask, s1_mask, s2_mask = _validate_force_group_masks(force_group_masks)
    integrator = openmm.CustomIntegrator(_resolve_timestep(timestep))
    _add_coupling_globals(
        integrator,
        coupling_amplitude_kj_mol=coupling_amplitude_kj_mol,
        coupling_width_per_kj2_mol2=coupling_width_per_kj2_mol2,
        diabatic_offset_kj_mol=diabatic_offset_kj_mol,
        umbrella_force_constant_per_kj_mol=umbrella_force_constant_per_kj_mol,
        umbrella_dv_center_kj_mol=umbrella_dv_center_kj_mol,
    )
    _add_mixing_globals(integrator)
    integrator.addPerDofVariable("F_adia", 0)
    integrator.addPerDofVariable("x1", 0)

    integrator.addUpdateContextState()
    # First half-kick: compute diabatic energies, mix, accumulate adiabatic force.
    _add_diabatic_energies_step(integrator, u_mask=u_mask, s1_mask=s1_mask, s2_mask=s2_mask)
    _add_mixing_step(integrator)
    _add_adiabatic_force_step(integrator, u_mask=u_mask, s1_mask=s1_mask, s2_mask=s2_mask)
    integrator.addComputePerDof("v", "v + 0.5*dt*F_adia/m")
    # Drift.
    integrator.addComputePerDof("x", "x + dt*v")
    integrator.addComputePerDof("x1", "x")
    integrator.addConstrainPositions()
    # Recompute everything at new positions for the second half-kick.
    _add_diabatic_energies_step(integrator, u_mask=u_mask, s1_mask=s1_mask, s2_mask=s2_mask)
    _add_mixing_step(integrator)
    _add_adiabatic_force_step(integrator, u_mask=u_mask, s1_mask=s1_mask, s2_mask=s2_mask)
    # Constraint correction folds in via (x - x1)/dt.
    integrator.addComputePerDof("v", "v + (x - x1)/dt + 0.5*dt*F_adia/m")
    integrator.addConstrainVelocities()
    return integrator


def build_langevin_middle_adiabatic_integrator(
    *,
    timestep: unit.Quantity | float,
    temperature: unit.Quantity | float,
    friction: unit.Quantity | float,
    force_group_masks: Tuple[int, int, int],
    coupling_amplitude_kj_mol: float,
    coupling_width_per_kj2_mol2: float = 0.0,
    diabatic_offset_kj_mol: float = 0.0,
    umbrella_force_constant_per_kj_mol: float = 0.0,
    umbrella_dv_center_kj_mol: float = 0.0,
) -> openmm.CustomIntegrator:
    """Langevin-Middle (single force eval/step) on the adiabatic surface.

    The drift–O–drift sandwich follows OpenMM's stock LangevinMiddleIntegrator
    pattern, with the deterministic kick driven by the adiabatically mixed
    force ``F_adia = fU + w1_eff·fS1 + w2_eff·fS2`` where ``w*_eff``
    optionally fold in a harmonic umbrella restraint on ``ΔV``.

    See :func:`build_nve_adiabatic_integrator` for the meaning of the
    ``umbrella_*`` parameters; ``K = 0`` (the default) is a no-op.
    """

    u_mask, s1_mask, s2_mask = _validate_force_group_masks(force_group_masks)
    if isinstance(temperature, unit.Quantity):
        temperature_q = temperature
    else:
        temperature_q = float(temperature) * unit.kelvin
    if isinstance(friction, unit.Quantity):
        friction_q = friction
    else:
        friction_q = float(friction) / unit.picosecond

    temperature_kelvin = float(temperature_q.value_in_unit(unit.kelvin))
    friction_per_ps = float(friction_q.value_in_unit(unit.picosecond ** -1))
    boltzmann_kj_per_mol_per_k = float(
        unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoules_per_mole / unit.kelvin)
    )
    kT_kj_per_mol = boltzmann_kj_per_mol_per_k * temperature_kelvin

    integrator = openmm.CustomIntegrator(_resolve_timestep(timestep))
    _add_coupling_globals(
        integrator,
        coupling_amplitude_kj_mol=coupling_amplitude_kj_mol,
        coupling_width_per_kj2_mol2=coupling_width_per_kj2_mol2,
        diabatic_offset_kj_mol=diabatic_offset_kj_mol,
        umbrella_force_constant_per_kj_mol=umbrella_force_constant_per_kj_mol,
        umbrella_dv_center_kj_mol=umbrella_dv_center_kj_mol,
    )
    integrator.addGlobalVariable("kT", kT_kj_per_mol)
    integrator.addGlobalVariable("friction", friction_per_ps)
    _add_mixing_globals(integrator)
    integrator.addPerDofVariable("F_adia", 0)
    integrator.addPerDofVariable("x1", 0)

    integrator.addUpdateContextState()
    _add_diabatic_energies_step(integrator, u_mask=u_mask, s1_mask=s1_mask, s2_mask=s2_mask)
    _add_mixing_step(integrator)
    _add_adiabatic_force_step(integrator, u_mask=u_mask, s1_mask=s1_mask, s2_mask=s2_mask)
    # Full kick using mixed adiabatic force.
    integrator.addComputePerDof("v", "v + dt*F_adia/m")
    # First half-drift.
    integrator.addComputePerDof("x", "x + 0.5*dt*v")
    # O step (Ornstein–Uhlenbeck) — friction in the middle of the step.
    integrator.addComputePerDof(
        "v",
        "exp(-friction*dt)*v + sqrt(kT*(1 - exp(-2*friction*dt))/m)*gaussian",
    )
    # Second half-drift.
    integrator.addComputePerDof("x", "x + 0.5*dt*v")
    integrator.addComputePerDof("x1", "x")
    integrator.addConstrainPositions()
    integrator.addComputePerDof("v", "v + (x - x1)/dt")
    integrator.addConstrainVelocities()
    return integrator


# AdiabaticReporter was moved to evbopenmm.reporters.adiabatic on 2026-07-15
# when the reporters sub-package was introduced. We re-export it here so
# existing imports (`from evbopenmm.sampling.adiabatic import AdiabaticReporter`
# and `from evbopenmm.sampling import AdiabaticReporter`) keep working.
from ..reporters.adiabatic import AdiabaticReporter, _TRACE_FIELDS, _format_value  # noqa: F401


__all__ = [
    "compute_adiabatic_mixing",
    "build_nve_adiabatic_integrator",
    "build_langevin_middle_adiabatic_integrator",
    "AdiabaticReporter",
]
