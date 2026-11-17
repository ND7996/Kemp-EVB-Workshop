from __future__ import annotations

import math
import warnings
from collections.abc import Callable, Iterable, Mapping, Sequence

import openmm
from openmm import unit

from ..types import CanonicalPairKey


def build_reactive_softpair_force(
    *,
    lambda_param: str,
    state_parameters: Mapping[CanonicalPairKey, tuple[float, float, float]],
    initial_state: str,
    particle_indices_from_canonical_pair: Callable[[str, CanonicalPairKey], tuple[int, int]],
    coulomb_constant_kj_mol_nm: float = 138.935456,
) -> openmm.CustomBondForce:
    """Reactive soft-pair: Buckingham repulsion + Coulomb on forming/breaking pairs.

    Validated against Q6's binary semantics (``src/nonbonded.f90`` ``nbe_qq``
    and ``src/simprep.f90`` ``precompute_set_values_qq``):

    ``V(r) = lambda_state * [ A * exp(-B * r) + k_coul * q_i * q_j / r ]``

    where:

    * ``A`` [kJ/mol] is the cross-atom soft-pair amplitude with *product*
      combining ``A_cross = C_i * C_j`` (per-atom ``C_i`` from
      ``soft_pair_parameters_by_atom_type``; no sqrt, no mean — matches the
      Q6 ``vdw=2`` softpair mode).
    * ``B`` [1/nm] is the cross-atom exponent with the same product combining
      ``B_cross = a_i * a_j``.
    * ``k_coul`` is the Coulomb constant in OpenMM units (kJ/mol·nm/e²).

    This form is bounded at short r (``exp(-B·r) → 0`` as ``r → 0``) and
    replaces the LJ + Coulomb for the pair entirely. The caller MUST also
    exclude the pair from the base NonbondedForce so there is no double-
    counted LJ/Coulomb; :meth:`EVBSystem.add_reactive_softpair_vdw_forces`
    handles that AND filters out reactive pairs that are 1-2 bonded (Morse) in
    the current state (the Morse term carries the bonded-state energy, so
    softpair Coulomb on the bonded pair would double-count).

    ``state_parameters`` maps ``CanonicalPairKey -> (A_kj_mol, B_per_nm, q_i*q_j_e2)``.
    The charge product uses the state-specific per-atom charges and may differ
    between states; this is intentional (each diabatic state contributes its
    own soft-pair Coulomb, matching the parity-validated form).
    """
    force = openmm.CustomBondForce(
        f"{lambda_param} * (A * exp(-B * r) + k_coul * qiqj / r)"
    )
    force.addGlobalParameter(lambda_param, 0.0)
    force.addGlobalParameter("k_coul", float(coulomb_constant_kj_mol_nm))
    force.addPerBondParameter("A")
    force.addPerBondParameter("B")
    force.addPerBondParameter("qiqj")
    for pair, (A_kj, B_per_nm, qiqj) in sorted(state_parameters.items()):
        idx_a, idx_b = particle_indices_from_canonical_pair(initial_state, pair)
        force.addBond(idx_a, idx_b, (float(A_kj), float(B_per_nm), float(qiqj)))
    return force


def build_reactive_morse_bond_force(
    *,
    lambda_param: str,
    state_parameters: Mapping[CanonicalPairKey, tuple[float, float, float]],
    initial_state: str,
    particle_indices_from_canonical_pair: Callable[[str, CanonicalPairKey], tuple[int, int]],
) -> openmm.CustomBondForce:
    """Lambda-scaled Morse bond force for forming/breaking reactive pairs.

    Energy form per bond::

        lambda_state * de * (1 - exp(-alpha * (r - r0)))^2

    Parameter tuple is ``(r0_nm, de_kj_mol, alpha_per_nm)``. The Morse
    parameters describe the bond's dissociation curve and **must come from
    a QM (or Q6 ``.fep``) source** — they are not derivable from a
    force-field parameterisation step, which would only produce a harmonic
    bond.

    Replaces the FF-derived harmonic entry on the reactive pair so the
    bonded-state energy reaches a finite dissociation asymptote rather
    than rising quadratically forever.
    """
    force = openmm.CustomBondForce(
        f"{lambda_param} * de * (1 - exp(-alpha * (r - r0)))^2"
    )
    force.addGlobalParameter(lambda_param, 0.0)
    force.addPerBondParameter("r0")
    force.addPerBondParameter("de")
    force.addPerBondParameter("alpha")
    for pair, (r0_nm, de_kj_mol, alpha_per_nm) in sorted(state_parameters.items()):
        idx_a, idx_b = particle_indices_from_canonical_pair(initial_state, pair)
        force.addBond(idx_a, idx_b, (float(r0_nm), float(de_kj_mol), float(alpha_per_nm)))
    return force


def build_reactive_harmonic_bond_force(
    *,
    lambda_param: str,
    regular_ff_param: str,
    state_parameters: Mapping[CanonicalPairKey, tuple[float, float]],
    initial_state: str,
    particle_indices_from_canonical_pair: Callable[[str, CanonicalPairKey], tuple[int, int]],
) -> openmm.CustomBondForce:
    """Harmonic bond matched to the Morse minimum, for *endpoint integration only*.

    Energy form per bond::

        regular_ff * lambda_state * 0.5 * k * (r - r0)^2

    Parameter tuple is ``(r0_nm, k_kj_mol_nm2)`` with ``k = 2 * de * alpha^2`` so
    the curvature matches the Morse at its minimum (``r0``). This is *not* part of
    the EVB diabatic energy: the global parameter ``regular_ff_param`` defaults to
    0, so the force contributes nothing to ``minimizeEnergy``/``getState`` and to
    the masked diabatic-energy sampling. It is switched on (``regular_ff_param=1``)
    ONLY when :meth:`EVBSystem.restrict_integration_to_endpoint` is called with
    ``regular_ff=True`` and is selected for integration in place of the Morse,
    which has a flattening (dissociative) restoring force that goes unstable at a
    2 fs timestep. Near ``r0`` the harmonic and Morse agree, so the endpoint
    ensemble stays consistent with the EVB diabat while integrating stably.
    """
    force = openmm.CustomBondForce(
        f"{regular_ff_param} * {lambda_param} * 0.5 * k * (r - r0)^2"
    )
    force.addGlobalParameter(regular_ff_param, 0.0)
    force.addGlobalParameter(lambda_param, 0.0)
    force.addPerBondParameter("r0")
    force.addPerBondParameter("k")
    for pair, (r0_nm, k_kj_mol_nm2) in sorted(state_parameters.items()):
        idx_a, idx_b = particle_indices_from_canonical_pair(initial_state, pair)
        force.addBond(idx_a, idx_b, (float(r0_nm), float(k_kj_mol_nm2)))
    return force


def promote_reactive_constraints_to_bonds(
    system: openmm.System,
    reactive_atom_pairs: Iterable[tuple[int, int]],
    *,
    warn: bool = True,
    default_k_kj_mol_nm2: float = 462750.4,
) -> int:
    """Convert HBonds *constraints* on reactive atom pairs into HarmonicBondForce terms.

    When a per-diabat ``System`` is built with ``constraints=HBonds``, any reactive
    bond that involves hydrogen (e.g. a transferring proton's donor O–H / acceptor
    N–H) is turned into an OpenMM **constraint** and removed from the
    ``HarmonicBondForce``. That breaks reactive EVB in two ways:

    1. Reactive-bond detection (``EVBSystem.get_forming_breaking_bond_pairs`` /
       ``_forming_breaking_bond_pairs_by_index``) compares *bond-force* terms, so it
       never sees the constrained reactive bond and the pair is silently dropped from
       the Morse setup.
    2. When the same atom (a shared transferring H) is constrained to different
       partners in different diabats, ``add_shared_constraints_to_evb_system`` keeps
       only constraints shared+equal across all states, so the conflicting reactive
       constraints are discarded — leaving the atom with *no* bonded restraint, so it
       dissociates during dynamics while the diabatic energy gap stays finite (the
       failure hides).

    This helper restores each such pair as a flexible harmonic bond so detection and
    Morse replacement operate normally. Call it on **each per-state System before
    constructing the EVBSystem**.

    Parameters
    ----------
    system : openmm.System
        Per-diabat System, mutated in place.
    reactive_atom_pairs : iterable of (int, int)
        Atom-index pairs of the reactive bonds (donor/acceptor of a transferring
        proton, forming/breaking heavy-atom bonds, ...). Order within a pair is
        irrelevant; pairs not currently constrained are ignored.
    warn : bool
        If True (default), emit a :class:`UserWarning` per constraint removed so the
        caller knows the bond is now FLEXIBLE (integrated explicitly, needs <= 1 fs)
        and will not be held rigid after the EVB is built.
    default_k_kj_mol_nm2 : float
        Force constant for the re-added harmonic bond. The reactive-bond machinery
        zeroes this ``k`` and installs the Morse, so it is only a placeholder; the
        default is a typical X–H value so the bond stays sane if Morse replacement is
        skipped.

    Returns
    -------
    int
        Number of constraints promoted to harmonic bonds.
    """
    targets = {frozenset((int(a), int(b))) for a, b in reactive_atom_pairs}
    if not targets:
        return 0
    hbond_force: openmm.HarmonicBondForce | None = None
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            hbond_force = force
            break
    if hbond_force is None:
        hbond_force = openmm.HarmonicBondForce()
        system.addForce(hbond_force)
    promoted = 0
    # Reverse order: removeConstraint renumbers all higher-index constraints.
    for i in range(system.getNumConstraints() - 1, -1, -1):
        atom1, atom2, distance = system.getConstraintParameters(i)
        if frozenset((atom1, atom2)) in targets:
            r0_nm = distance.value_in_unit(unit.nanometer)
            hbond_force.addBond(int(atom1), int(atom2), float(r0_nm),
                                float(default_k_kj_mol_nm2))
            system.removeConstraint(i)
            promoted += 1
            if warn:
                warnings.warn(
                    "promote_reactive_constraints_to_bonds: removed the HBonds "
                    f"constraint on reactive pair ({atom1}, {atom2}) "
                    f"(r0={r0_nm:.4f} nm) and re-added it as a flexible harmonic bond. "
                    "This reactive bond is now integrated explicitly (use <= 1 fs) and "
                    "will NOT be held rigid once the EVB Context is built.",
                    stacklevel=2,
                )
    return promoted


__all__ = [
    "build_reactive_harmonic_bond_force",
    "build_reactive_morse_bond_force",
    "build_reactive_softpair_force",
    "promote_reactive_constraints_to_bonds",
]


# ---------------------------------------------------------------------------
# Q6 soft-pair parameter sanity
# ---------------------------------------------------------------------------
#
# PROVENANCE of the conventional per-element (C, a) values, since they are easy to
# copy without knowing which are measured and which are convention:
#
#   C = 91.0, a = 2.5   carbon    -- Q6 ``ccl_cl.fep``; the only value with a
#                                    documented upstream source.
#   C = 85.0, a = 2.5   oxygen    -- convention, carried between systems.
#   C = 40.0, a = 2.0   hydrogen  -- an estimate for a transferring proton.
#   default 90.0 / 2.5            -- Q6's generic per-atom values.
#
# They are NOT fitted to the Lennard-Jones they replace, and must not be. The soft-pair
# exists to be much *softer* than LJ at bonding distance so the bond can form; fitting
# it to the r**-12 wall reproduces the wall. Two independent reasons this cannot work:
# a polar transferring hydrogen typically has sigma = epsilon = 0 in a protein force
# field, so there is no LJ to fit at all; and for pairs that do have LJ, any fit over a
# non-bonded range must be extrapolated inward across the steepest part of r**-12, so
# the fit window sets the answer rather than the physics.
#
# The cheap check that does work is below: evaluate the soft-pair at the forming bond's
# Morse r0. It should be a few kcal/mol at most.

_DEFAULT_SP_C_KCAL_MOL = 90.0
_DEFAULT_SP_A_PER_ANGSTROM = 2.5


def soft_pair_energy_kcal_mol(
    element_i: str,
    element_j: str,
    r_nm: float,
    *,
    parameters_by_atom_type: Mapping[str, tuple[float, float]] | None = None,
    default_c_kcal_mol: float = _DEFAULT_SP_C_KCAL_MOL,
    default_a_per_angstrom: float = _DEFAULT_SP_A_PER_ANGSTROM,
) -> float:
    """Q6 soft-pair repulsion ``A*exp(-B*r)`` for a pair, in kcal/mol.

    Per-atom values combine as ``A = C_i*C_j`` and ``B = a_i*a_j``, matching
    :meth:`EVBSystem.add_reactive_softpair_vdw_forces`. Elements absent from
    ``parameters_by_atom_type`` fall back to the Q6 defaults.
    """
    table = parameters_by_atom_type or {}
    c_i, a_i = table.get(element_i, (default_c_kcal_mol, default_a_per_angstrom))
    c_j, a_j = table.get(element_j, (default_c_kcal_mol, default_a_per_angstrom))
    a_cross = float(c_i) * float(c_j)
    b_cross = float(a_i) * float(a_j)
    return a_cross * math.exp(-b_cross * (float(r_nm) * 10.0))


def check_soft_pair_at_bond_distances(
    pairs: Iterable[tuple[str, str, float]],
    *,
    parameters_by_atom_type: Mapping[str, tuple[float, float]] | None = None,
    max_kcal_mol: float = 50.0,
    labels: Sequence[str] | None = None,
    raise_on_error: bool = False,
    warn: bool = True,
) -> dict[str, float]:
    """Evaluate the soft-pair at each forming bond's equilibrium distance.

    A soft-pair that is stiff at the distance where the bond is supposed to form
    rebuilds the very wall it exists to remove. Because it is a pure exponential, this
    is checkable in closed form before running anything.

    The check is sharpest for a **forming** bond: the soft-pair acts in the state where
    that pair is not yet bonded, so r0 is exactly the distance the atoms must reach.
    For a **breaking** bond it is a worst case rather than an operating point -- the
    soft-pair acts in the state where the pair has already separated, so the atoms
    typically sit well beyond r0 and never pay the value reported here.

    Scale expectations accordingly. With the conventional per-element values a forming
    heavy-atom bond lands around 1 kcal/mol, while a pair involving the transferring
    hydrogen evaluates to a few tens at its bonded r0 and is still healthy. Hundreds
    means the soft-pair has rebuilt the wall.

    Parameters
    ----------
    pairs :
        ``(element_i, element_j, r0_nm)`` per reactive bond -- ``r0_nm`` being the Morse
        equilibrium distance, i.e. where the bond forms.
    parameters_by_atom_type :
        Per-element ``(C, a)``; missing elements use the Q6 defaults.
    max_kcal_mol :
        Threshold above which the soft-pair is reported as too stiff. The default is
        deliberately loose enough to pass the conventional per-element values (whose
        X-H pairs reach ~30 kcal/mol at their bonded r0) while still catching a
        soft-pair fitted to Lennard-Jones, which lands in the hundreds.
    labels :
        Optional names for the pairs, used in the message.

    Returns
    -------
    dict[str, float]
        ``{label: energy_kcal_mol}`` for every pair checked.
    """
    entries = list(pairs)
    if labels is None:
        names = [f"{i}-{j} @ {r:.4f} nm" for i, j, r in entries]
    else:
        names = list(labels)
        if len(names) != len(entries):
            raise ValueError(
                f"labels has {len(names)} entries but {len(entries)} pairs were given; "
                "zipping them would silently drop pairs from the check."
            )
    energies: dict[str, float] = {}
    offenders: list[str] = []
    for name, (element_i, element_j, r_nm) in zip(names, entries):
        value = soft_pair_energy_kcal_mol(
            element_i, element_j, r_nm,
            parameters_by_atom_type=parameters_by_atom_type,
        )
        energies[name] = value
        if value > max_kcal_mol:
            offenders.append(f"    {name}: {value:.1f} kcal/mol")

    if offenders:
        message = (
            "Reactive soft-pair is stiff at the distance where the bond forms "
            f"(> {max_kcal_mol} kcal/mol). The soft-pair replaces the pair's "
            "Lennard-Jones precisely so the atoms can reach bonding distance; if it is "
            "this repulsive at the Morse r0 it rebuilds the wall it was meant to "
            "remove, and the reaction may not be sampled at all. Use softer per-element "
            "(C, a) values.\n" + "\n".join(offenders)
        )
        if raise_on_error:
            raise RuntimeError(message)
        if warn:
            warnings.warn(message, stacklevel=2)
    return energies
