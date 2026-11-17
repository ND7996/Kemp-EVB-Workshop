"""OPLS-AA/M ligand parameterization for EVB per-state systems.

Ingest a Schrodinger ``ffld_server`` OPLS-AA/M assignment (its text output) and
patch an existing OpenMM :class:`~openmm.System` so a chosen ligand carries full
OPLS-AA/M parameters: Lennard-Jones (sigma/epsilon), the complete bonded set
(bond/angle/proper/improper), and OPLS 1-4 exception scaling with geometric
mixing. Ligand<->environment cross terms are left at Lorentz-Berthelot -- the
accepted "Tier A" approximation in which the ligand-*internal* nonbonded is full
OPLS while cross-LJ keeps the host System's combining rule.

``ffld_server`` (Academic Maestro) is a *local preparation* step: this module
only parses its text output, so the library never depends on Maestro at runtime.
The parsed force field is a plain dataclass and is independently testable.

Two stages, deliberately separated so each is testable on its own:

* :func:`parse_ffld` -- text -> :class:`FfldFF` (no unit conversion).
* :func:`apply_opls_ligand` -- :class:`FfldFF` + a live System -> patched System.

Atom order in the ffld VDW block follows the input SDF/MOL atom order, so the
ligand is assumed to occupy a contiguous index block ``[ligand_start,
ligand_start + len(ffld.atoms))`` in the System, in that same order. Use
:func:`ligand_start_from_topology` to locate it from an OpenMM topology.

Unit conventions in the ffld file (OPLS/Desmond) and the OpenMM conversions:

===============  ==========================  ================================
quantity         ffld unit                   OpenMM conversion
===============  ==========================  ================================
sigma            Angstrom                    ``* 0.1`` -> nm
epsilon          kcal/mol                    ``* 4.184`` -> kJ/mol
bond ``k``       kcal/mol/A^2 (no 1/2)       ``* 836.8`` -> kJ/mol/nm^2 (1/2)
bond ``r0``      Angstrom                    ``* 0.1`` -> nm
angle ``k``      kcal/mol/rad^2 (no 1/2)     ``* 8.368`` -> kJ/mol/rad^2 (1/2)
angle ``theta0`` degrees                     ``* pi/180`` -> radian
torsion ``Vn``   kcal/mol (OPLS Fourier)     ``Vn/2 * 4.184`` -> per-term k
===============  ==========================  ================================

The OPLS proper Fourier form ``E = sum_n Vn/2 [1 +- cos(n phi)]`` maps to
OpenMM ``PeriodicTorsionForce`` terms with phase ``0`` for odd ``n`` and ``pi``
for even ``n``. An ffld improper is a single ``n=2`` term (its central atom is
the *third* index). OPLS 1-4 scaling uses ``fudgeQQ = fudgeLJ = 0.5`` with
geometric LJ mixing: ``chargeProd = 0.5 qi qj``, ``sigma14 = sqrt(si sj)``,
``eps14 = 0.5 sqrt(ei ej)``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import openmm
from openmm import app, unit

__all__ = [
    "FfldAtom",
    "FfldFF",
    "OplsPatchReport",
    "parse_ffld",
    "apply_opls_ligand",
    "ligand_start_from_topology",
    "ligand_system_energy",
]

# --- unit conversions (ffld OPLS/Desmond -> OpenMM native) ---
_KCAL_TO_KJ = 4.184
_DEG_TO_RAD = math.pi / 180.0
# OPLS harmonic potentials carry no 1/2; OpenMM Harmonic{Bond,Angle} do -> x2.
_BOND_K_FACTOR = 2.0 * _KCAL_TO_KJ / (0.1 ** 2)   # 836.8  kcal/mol/A^2 -> kJ/mol/nm^2
_ANGLE_K_FACTOR = 2.0 * _KCAL_TO_KJ               # 8.368  kcal/mol/rad^2 -> kJ/mol/rad^2


# =====================================================================
# Stage 1 -- parse ffld_server text
# =====================================================================
@dataclass
class FfldAtom:
    """One atom's OPLS-AA/M assignment from an ffld_server VDW block."""

    name: str
    opls_type: int
    vdw_symbol: str
    charge: float          # e
    sigma: float           # Angstrom
    epsilon: float         # kcal/mol
    comment: str = ""


@dataclass
class FfldFF:
    """Parsed OPLS-AA/M force field for a single molecule.

    ``atoms`` is in ffld/SDF order. Bonded lists reference atoms by their ffld
    ``name`` (resolve to indices with :meth:`index_of`).
    """

    molecule: str
    total_charge: float
    atoms: List[FfldAtom]
    stretch: List[Tuple[str, str, float, float]] = field(default_factory=list)          # a,b,k,r0
    bending: List[Tuple[str, str, str, float, float]] = field(default_factory=list)     # a,b,c,k,th0
    proper: List[Tuple[str, str, str, str, float, float, float, float]] = field(default_factory=list)
    improper: List[Tuple[str, str, str, str, float]] = field(default_factory=list)      # a,b,c,d,V2

    def index_of(self) -> Dict[str, int]:
        """Map ffld atom name -> 0-based position in :attr:`atoms`."""
        return {a.name: i for i, a in enumerate(self.atoms)}


def _sections(lines: List[str]) -> Dict[str, Tuple[int, int]]:
    """Return ``{section: (header_line, end)}`` for the recognised blocks."""
    marks = []
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("OPLSAA FORCE FIELD TYPE ASSIGNED"):
            marks.append(("vdw", i))
        elif s.startswith("Stretch"):
            marks.append(("stretch", i))
        elif s.startswith("Bending"):
            marks.append(("bending", i))
        elif s.startswith("proper Torsion"):
            marks.append(("proper", i))
        elif s.startswith("improper Torsion"):
            marks.append(("improper", i))
    bounds = {}
    for j, (name, i) in enumerate(marks):
        end = marks[j + 1][1] if j + 1 < len(marks) else len(lines)
        bounds[name] = (i, end)
    return bounds


def _is_dashes(s: str) -> bool:
    t = s.strip()
    return bool(t) and set(t) <= set("-=")


def parse_ffld(path: str | Path) -> FfldFF:
    """Parse an ``ffld_server`` OPLS-AA/M text output into an :class:`FfldFF`.

    Raises ``ValueError`` if the mandatory VDW / per-atom block is absent.
    """
    lines = Path(path).read_text().splitlines()
    # molecule name: first non-blank line after the "for molecule:" banner
    mol = "unknown"
    for i, ln in enumerate(lines):
        if "for molecule" in ln:
            for k in range(i + 1, min(i + 4, len(lines))):
                if lines[k].strip():
                    mol = lines[k].strip()
                    break
            break
    total_q = 0.0
    for ln in lines:
        if "Total charge of the system" in ln:
            total_q = float(ln.split("=")[-1])
            break

    b = _sections(lines)
    if "vdw" not in b:
        raise ValueError(f"{path}: no 'OPLSAA FORCE FIELD TYPE ASSIGNED' block found")

    # ---- VDW / per-atom block ----
    atoms: List[FfldAtom] = []
    s, e = b["vdw"]
    for ln in lines[s:e]:
        t = ln.rstrip()
        if not t.strip() or _is_dashes(t):
            continue
        p = t.split()
        if p[0].lower() == "atom":            # header row
            continue
        if len(p) < 7:                        # name type vdw symbol charge sigma epsilon [...]
            continue
        try:
            atoms.append(FfldAtom(
                name=p[0], opls_type=int(p[1]), vdw_symbol=p[2], charge=float(p[4]),
                sigma=float(p[5]), epsilon=float(p[6]),
                comment=" ".join(p[8:]) if len(p) > 8 else ""))
        except ValueError:
            continue

    ff = FfldFF(molecule=mol, total_charge=total_q, atoms=atoms)

    def rows(section: str):
        if section not in b:
            return
        s, e = b[section]
        for ln in lines[s:e]:
            t = ln.rstrip()
            if not t.strip() or _is_dashes(t):
                continue
            p = t.split()
            if p[0] in (section.capitalize(), "Stretch", "Bending", "proper", "improper"):
                continue
            yield p

    for p in rows("stretch"):
        if len(p) >= 4:
            try:
                ff.stretch.append((p[0], p[1], float(p[2]), float(p[3])))
            except ValueError:
                pass
    for p in rows("bending"):
        if len(p) >= 5:
            try:
                ff.bending.append((p[0], p[1], p[2], float(p[3]), float(p[4])))
            except ValueError:
                pass
    for p in rows("proper"):
        if len(p) >= 8:
            try:
                ff.proper.append((p[0], p[1], p[2], p[3],
                                  float(p[4]), float(p[5]), float(p[6]), float(p[7])))
            except ValueError:
                pass
    for p in rows("improper"):
        if len(p) >= 5:
            try:
                ff.improper.append((p[0], p[1], p[2], p[3], float(p[4])))
            except ValueError:
                pass
    return ff


# =====================================================================
# Stage 2 -- patch an OpenMM System's ligand to full OPLS-AA/M
# =====================================================================
@dataclass
class OplsPatchReport:
    """Summary of what :func:`apply_opls_ligand` changed."""

    ligand_start: int
    n_atoms: int
    n_14: int
    n_exclusions: int
    bonds_replaced: int
    bonds_missing: Tuple[Tuple[int, int], ...]
    angles_replaced: int
    angles_missing: Tuple[Tuple[int, int, int], ...]
    propers: int
    impropers: int


def _single_force(system: openmm.System, cls) -> object:
    forces = [f for f in system.getForces() if isinstance(f, cls)]
    if len(forces) != 1:
        raise ValueError(
            f"apply_opls_ligand expects exactly one {cls.__name__}; found {len(forces)}")
    return forces[0]


def _patch_nonbonded(system, ff: FfldFF, base: int, n: int, set_charges: bool):
    nb = _single_force(system, openmm.NonbondedForce)
    sig: Dict[int, float] = {}
    eps: Dict[int, float] = {}
    q: Dict[int, float] = {}
    for k in range(n):
        a = ff.atoms[k]
        gi = base + k
        oldq, _, _ = nb.getParticleParameters(gi)
        newq = a.charge if set_charges else oldq.value_in_unit(unit.elementary_charge)
        s = a.sigma * 0.1
        e = a.epsilon * _KCAL_TO_KJ
        nb.setParticleParameters(gi, newq * unit.elementary_charge,
                                 s * unit.nanometer, e * unit.kilojoule_per_mole)
        sig[gi] = s
        eps[gi] = e
        q[gi] = newq
    # rebuild ligand 1-4 exceptions with OPLS scaling + geometric mixing;
    # leave 1-2 / 1-3 full exclusions untouched
    n14 = nexcl = 0
    for ei in range(nb.getNumExceptions()):
        i, j, cprod, s14, e14 = nb.getExceptionParameters(ei)
        if not (base <= i < base + n and base <= j < base + n):
            continue
        is_excl = (abs(cprod.value_in_unit(unit.elementary_charge ** 2)) < 1e-12 and
                   e14.value_in_unit(unit.kilojoule_per_mole) < 1e-12)
        if is_excl:
            nexcl += 1
            continue
        cp = 0.5 * q[i] * q[j]
        sg = math.sqrt(sig[i] * sig[j])
        ep = 0.5 * math.sqrt(eps[i] * eps[j])
        nb.setExceptionParameters(ei, i, j, cp * unit.elementary_charge ** 2,
                                  sg * unit.nanometer, ep * unit.kilojoule_per_mole)
        n14 += 1
    return n14, nexcl


def _rebuild_bonds(system, ff: FfldFF, base: int, n: int):
    old = _single_force(system, openmm.HarmonicBondForce)
    new = openmm.HarmonicBondForce()
    new.setUsesPeriodicBoundaryConditions(old.usesPeriodicBoundaryConditions())
    io = ff.index_of()
    ff_bond = {tuple(sorted((base + io[a], base + io[b]))): (r0 * 0.1, k * _BOND_K_FACTOR)
               for a, b, k, r0 in ff.stretch}
    repl = 0
    missing: List[Tuple[int, int]] = []
    for bi in range(old.getNumBonds()):
        i, j, r0, k = old.getBondParameters(bi)
        if base <= i < base + n and base <= j < base + n:
            key = tuple(sorted((i, j)))
            if key in ff_bond:
                r0n, kn = ff_bond[key]
                new.addBond(i, j, r0n * unit.nanometer,
                            kn * unit.kilojoule_per_mole / unit.nanometer ** 2)
                repl += 1
            else:
                missing.append((i - base, j - base))
                new.addBond(i, j, r0, k)
        else:
            new.addBond(i, j, r0, k)
    return new, repl, missing


def _rebuild_angles(system, ff: FfldFF, base: int, n: int):
    old = _single_force(system, openmm.HarmonicAngleForce)
    new = openmm.HarmonicAngleForce()
    new.setUsesPeriodicBoundaryConditions(old.usesPeriodicBoundaryConditions())
    io = ff.index_of()
    ff_ang: Dict[Tuple[int, int, int], Tuple[float, float]] = {}
    for a, b, c, k, th in ff.bending:
        key = (base + io[a], base + io[b], base + io[c])
        ff_ang[key] = (th * _DEG_TO_RAD, k * _ANGLE_K_FACTOR)
        ff_ang[key[::-1]] = (th * _DEG_TO_RAD, k * _ANGLE_K_FACTOR)
    repl = 0
    missing: List[Tuple[int, int, int]] = []
    for ai in range(old.getNumAngles()):
        i, j, k_, th, kk = old.getAngleParameters(ai)
        if all(base <= x < base + n for x in (i, j, k_)):
            key = (i, j, k_)
            if key in ff_ang:
                thn, kn = ff_ang[key]
                new.addAngle(i, j, k_, thn * unit.radian,
                             kn * unit.kilojoule_per_mole / unit.radian ** 2)
                repl += 1
            else:
                missing.append((i - base, j - base, k_ - base))
                new.addAngle(i, j, k_, th, kk)
        else:
            new.addAngle(i, j, k_, th, kk)
    return new, repl, missing


def _rebuild_torsions(system, ff: FfldFF, base: int, n: int):
    old = _single_force(system, openmm.PeriodicTorsionForce)
    new = openmm.PeriodicTorsionForce()
    new.setUsesPeriodicBoundaryConditions(old.usesPeriodicBoundaryConditions())
    io = ff.index_of()
    # keep non-ligand torsions; drop ligand ones (replaced by OPLS below)
    for ti in range(old.getNumTorsions()):
        i, j, k, l, per, ph, kk = old.getTorsionParameters(ti)
        if not all(base <= x < base + n for x in (i, j, k, l)):
            new.addTorsion(i, j, k, l, per, ph, kk)
    nprop = nimp = 0
    for a, b, c, d, V1, V2, V3, V4 in ff.proper:
        i, j, k, l = (base + io[a], base + io[b], base + io[c], base + io[d])
        for per, V, phase in ((1, V1, 0.0), (2, V2, math.pi), (3, V3, 0.0), (4, V4, math.pi)):
            if abs(V) > 1e-9:
                new.addTorsion(i, j, k, l, per, phase * unit.radian,
                               (V / 2.0 * _KCAL_TO_KJ) * unit.kilojoule_per_mole)
                nprop += 1
    for a, b, c, d, V2 in ff.improper:
        i, j, k, l = (base + io[a], base + io[b], base + io[c], base + io[d])
        if abs(V2) > 1e-9:
            new.addTorsion(i, j, k, l, 2, math.pi * unit.radian,
                           (V2 / 2.0 * _KCAL_TO_KJ) * unit.kilojoule_per_mole)
            nimp += 1
    return new, nprop, nimp


def _swap_single_force(system: openmm.System, cls, new) -> None:
    """Remove the System's single force of type ``cls`` and append ``new``.

    Identity (``is``) comparison is unreliable for OpenMM force wrappers (each
    ``getForce`` returns a fresh proxy), so remove by type. Safe because this is
    called before ``new`` is added, when exactly one ``cls`` force exists.
    """
    for idx in range(system.getNumForces()):
        if isinstance(system.getForce(idx), cls):
            system.removeForce(idx)
            break
    system.addForce(new)


def apply_opls_ligand(
    system: openmm.System,
    ffld: FfldFF,
    ligand_start: int,
    *,
    set_charges: bool = False,
) -> OplsPatchReport:
    """Patch ``system`` in place so the ligand carries full OPLS-AA/M parameters.

    The ligand is the contiguous atom block
    ``[ligand_start, ligand_start + len(ffld.atoms))`` in ffld/SDF order. Its
    per-atom LJ, bonded terms, and 1-4 exceptions are overwritten with the
    OPLS-AA/M values; 1-2/1-3 exclusions, all non-ligand terms, and (by default)
    the existing per-atom charges are left untouched.

    Parameters
    ----------
    system
        An OpenMM System with exactly one each of ``NonbondedForce``,
        ``HarmonicBondForce``, ``HarmonicAngleForce`` and
        ``PeriodicTorsionForce`` (raises ``ValueError`` otherwise).
    ffld
        Parsed OPLS-AA/M force field (see :func:`parse_ffld`).
    ligand_start
        Index of the first ligand atom in the System.
    set_charges
        If True, also overwrite ligand charges with the ffld values; if False
        (default) keep the System's charges (assumed already OPLS).

    Returns
    -------
    OplsPatchReport
        Counts of replaced / missing terms; a non-empty ``bonds_missing`` or
        ``angles_missing`` means the ffld lacked a term for a ligand
        bond/angle present in the System (that term keeps its original value).
    """
    n = len(ffld.atoms)
    if n == 0:
        raise ValueError("ffld force field has no atoms")
    if ligand_start < 0:
        raise ValueError(f"ligand_start must be >= 0, got {ligand_start}")
    if ligand_start + n > system.getNumParticles():
        raise ValueError(
            f"ligand block [{ligand_start}, {ligand_start + n}) exceeds System size "
            f"{system.getNumParticles()}")

    n14, nexcl = _patch_nonbonded(system, ffld, ligand_start, n, set_charges)
    # build replacement forces (reading the originals) before swapping any in
    nb_bonds, brepl, bmiss = _rebuild_bonds(system, ffld, ligand_start, n)
    nb_ang, arepl, amiss = _rebuild_angles(system, ffld, ligand_start, n)
    nb_tor, nprop, nimp = _rebuild_torsions(system, ffld, ligand_start, n)
    _swap_single_force(system, openmm.HarmonicBondForce, nb_bonds)
    _swap_single_force(system, openmm.HarmonicAngleForce, nb_ang)
    _swap_single_force(system, openmm.PeriodicTorsionForce, nb_tor)
    return OplsPatchReport(
        ligand_start=ligand_start, n_atoms=n, n_14=n14, n_exclusions=nexcl,
        bonds_replaced=brepl, bonds_missing=tuple(bmiss),
        angles_replaced=arepl, angles_missing=tuple(amiss),
        propers=nprop, impropers=nimp)


def ligand_start_from_topology(topology: app.Topology, residue_name: str) -> int:
    """Return the first atom index of ``residue_name`` in ``topology``.

    Validates that all atoms of that residue form a single contiguous block
    (as OPLS/ffld ordering requires). Raises ``ValueError`` if the residue is
    absent or its atoms are not contiguous.
    """
    idx = [a.index for a in topology.atoms() if a.residue.name == residue_name]
    if not idx:
        raise ValueError(f"residue '{residue_name}' not found in topology")
    if idx != list(range(min(idx), max(idx) + 1)):
        raise ValueError(f"residue '{residue_name}' atoms are not contiguous: {idx}")
    return min(idx)


def ligand_system_energy(
    system: openmm.System,
    positions,
    ligand_start: int,
    n_atoms: int,
) -> Tuple[float, float]:
    """Return ``(potential_energy_kj_mol, ligand_net_charge_e)`` for a quick check.

    Builds a throwaway CPU/Reference context, so it is a validation helper, not
    a hot-path routine.
    """
    integ = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    try:
        plat = openmm.Platform.getPlatformByName("CPU")
    except Exception:
        plat = openmm.Platform.getPlatformByName("Reference")
    ctx = openmm.Context(system, integ, plat)
    try:
        ctx.setPositions(positions)
        e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        nb = _single_force(system, openmm.NonbondedForce)
        qlig = sum(nb.getParticleParameters(ligand_start + k)[0].value_in_unit(unit.elementary_charge)
                   for k in range(n_atoms))
    finally:
        del ctx, integ
    return e, qlig
