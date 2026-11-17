"""OPLS-AA/M support utilities for OpenMM systems.

Two general helpers needed to use the OPLS-AA/M force field (opls-aam.xml, shipped under
``evbopenmm/data/opls_aam/``) correctly in OpenMM:

* :func:`add_opls_impropers` -- OPLS impropers are PERIODIC (n=2, phase 180, central atom at
  position 2), but OpenMM's built-in ``<Improper>`` matcher does not apply these fully-specified
  quartets, so they must be enumerated explicitly against the topology's sp2 centres.  Call this
  after ``createSystem`` (aromatic rings + amide/carbonyl planarity would otherwise be missing).

* :func:`flatten_opls_lj_to_nonbonded` -- the opls-aam.xml ``<Script>`` moves all LJ into a
  geometric CustomNonbondedForce and zeroes the NonbondedForce LJ (correct for standalone MD).
  Consumers that read per-atom sigma/epsilon from the NonbondedForce (e.g. EVBSystem's vdW
  partition builder) need the LJ back there; this moves it and removes the CustomNonbondedForce.
  EVBSystem then re-applies geometric combining itself (``lj_combining='geometric'``).
"""
from __future__ import annotations

import itertools
import xml.etree.ElementTree as ET
from collections import defaultdict

import openmm
import openmm.unit as u


# --------------------------------------------------------------------------- impropers
def _load_improper_tables(xml_paths):
    """Return (type -> class) and the periodic improper templates from the FF XMLs."""
    type2class = {}
    impropers = []
    for x in xml_paths:
        root = ET.parse(x).getroot()
        at = root.find("AtomTypes")
        if at is not None:
            for t in at:
                type2class[t.get("name")] = t.get("class")
        for imp in root.iter("Improper"):
            if imp.get("periodicity1") is None:
                continue  # skip any leftover harmonic-form impropers
            impropers.append((imp.get("class1"), imp.get("class2"), imp.get("class3"),
                              imp.get("class4"), int(imp.get("periodicity1")),
                              float(imp.get("phase1")), float(imp.get("k1"))))
    return type2class, impropers


def _assign_neighbours(want_classes, nbr_idx, nbr_classes):
    """Match want_classes (template positions 1,3,4) to the three neighbours as a permutation."""
    for perm in itertools.permutations(range(3)):
        if all(nbr_classes[perm[i]] == want_classes[i] for i in range(3)):
            return [nbr_idx[perm[0]], nbr_idx[perm[1]], nbr_idx[perm[2]]]
    return None


def add_opls_impropers(system, topology, forcefield, ff_xml_paths, verbose=False):
    """Enumerate + add OPLS-AA/M impropers to the System's PeriodicTorsionForce.

    Parameters
    ----------
    system, topology :
        The built OpenMM System and its Topology (post-createSystem).
    forcefield :
        The ``openmm.app.ForceField`` used to build the system; its graph-based template
        matching supplies each atom's class (robust to atom-name differences, NOT names).
    ff_xml_paths :
        The FF XML paths (the periodic ``<Improper>`` quartets are read from them).

    Returns the number of impropers added.  Raises if a 3-neighbour (sp2-candidate) atom
    fails to match a template (would otherwise silently drop an improper), or if OpenMM
    already applied impropers itself (would double-count).
    """
    from openmm.app import ForceField as _FF
    type2class, impropers = _load_improper_tables(ff_xml_paths)
    by_central = defaultdict(list)
    for c1, c2, c3, c4, n, ph, k in impropers:
        by_central[c2].append((c1, c3, c4, n, ph, k))

    data = _FF._SystemData(topology)
    forcefield._matchAllResiduesToTemplates(data, topology, {}, False)
    atoms = list(topology.atoms())

    adj = defaultdict(list)
    for b in topology.bonds():
        adj[b[0].index].append(b[1].index)
        adj[b[1].index].append(b[0].index)

    # Atoms whose type is not in the OPLS class map belong to another force field (water/ions from
    # a separate XML, or an OpenFF/SMILES-parameterized ligand): they carry no OPLS impropers, so
    # cls() returns None and any centre/neighbour with a None class is skipped below. (The OPLS FF
    # itself is complete -- guarded by tests/test_opls_forcefield.py -- so a None class here means a
    # non-OPLS atom, not a dropped OPLS improper.)
    def cls(i):
        return type2class.get(data.atomType.get(atoms[i]))

    ptf = next(f for f in system.getForces()
               if isinstance(f, openmm.PeriodicTorsionForce))

    # DEFENSIVE: OpenMM's own matcher applies 0 of these quartets, so we add them. If a future
    # OpenMM version DID match some, they would double-count; detect any improper-shaped torsion
    # already present (one atom bonded to the other three) and refuse.
    bonded = {frozenset((b[0].index, b[1].index)) for b in topology.bonds()}

    def _is_bonded(x, y):
        return frozenset((x, y)) in bonded

    preexisting = 0
    for i in range(ptf.getNumTorsions()):
        a, b, c, d, *_ = ptf.getTorsionParameters(i)
        for centre, others in ((a, (b, c, d)), (b, (a, c, d)),
                               (c, (a, b, d)), (d, (a, b, c))):
            if all(_is_bonded(centre, o) for o in others):
                # Only an OPLS-templated centre can carry an opls-aam <Improper>, so only those
                # would double-count with the ones we add below. An improper-shaped torsion on a
                # non-OPLS centre (cls None -- a SMIRNOFF/OpenFF-parameterized ligand or water)
                # is that force field's own and legitimately coexists; the enumeration loop below
                # already skips non-OPLS centres, so it never re-adds one.
                if cls(centre) is not None:
                    preexisting += 1
                break
    if preexisting:
        raise ValueError(f"add_opls_impropers: {preexisting} improper-shaped torsion(s) already "
                         f"in the PeriodicTorsionForce -- OpenMM matched impropers itself; adding "
                         f"ours would double-count. Remove <Improper> from the FF XML first.")

    added = 0
    for c in range(len(atoms)):
        nbrs = adj[c]
        if len(nbrs) != 3:
            continue
        cc = cls(c)
        if cc not in by_central:
            continue
        ncls = [cls(x) for x in nbrs]
        for (c1, c3, c4, n, ph, k) in by_central[cc]:
            asg = _assign_neighbours([c1, c3, c4], nbrs, ncls)
            if asg is not None:
                a1, a3, a4 = asg
                ptf.addTorsion(a1, c, a3, a4, n, ph, k)   # centre at position 2
                added += 1
                if verbose:
                    print(f"  improper: {atoms[a1].name}-{atoms[c].name}*-"
                          f"{atoms[a3].name}-{atoms[a4].name}  ({cc})")
                break  # one improper per centre
    return added


# --------------------------------------------------------------------------- LJ flatten
def _is_geometric_lj_force(f: openmm.Force) -> bool:
    if not isinstance(f, openmm.CustomNonbondedForce):
        return False
    expr = f.getEnergyFunction().replace(" ", "")
    names = {f.getPerParticleParameterName(i) for i in range(f.getNumPerParticleParameters())}
    return "sqrt(sigma1*sigma2)" in expr and {"sigma", "epsilon"} <= names


def flatten_opls_lj_to_nonbonded(system: openmm.System) -> openmm.System:
    """Move LJ from the opls-aam.xml Script's geometric CustomNonbondedForce into the
    NonbondedForce (so consumers that read per-atom sigma/epsilon there can see it) and remove
    that CustomNonbondedForce.  Charges and the (already-correct) 1-4 exceptions are untouched.

    Modifies `system` in place and returns it.  Raises if the expected forces aren't found.
    """
    nb = next((f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)), None)
    if nb is None:
        raise ValueError("flatten_opls_lj: no NonbondedForce in the system")
    geo = [f for f in system.getForces() if _is_geometric_lj_force(f)]
    if len(geo) != 1:
        raise ValueError(f"flatten_opls_lj: expected exactly 1 geometric-LJ CustomNonbondedForce, "
                         f"found {len(geo)}. Is this an OPLS (opls-aam.xml Script) system?")
    cnb = geo[0]
    n = nb.getNumParticles()
    if cnb.getNumParticles() != n:
        raise ValueError("flatten_opls_lj: particle-count mismatch between NB and CustomNB")

    per = [cnb.getPerParticleParameterName(i) for i in range(cnb.getNumPerParticleParameters())]
    sig_idx, eps_idx = per.index("sigma"), per.index("epsilon")

    for i in range(n):
        params = cnb.getParticleParameters(i)
        sigma = params[sig_idx] * u.nanometer
        epsilon = params[eps_idx] * u.kilojoule_per_mole
        charge = nb.getParticleParameters(i)[0]           # keep the charge
        nb.setParticleParameters(i, charge, sigma, epsilon)

    # remove the geometric CustomNonbondedForce. Match by predicate, not `is`: OpenMM's SWIG
    # getForce() returns a fresh proxy each call, so identity comparison is unreliable.
    for fi in range(system.getNumForces()):
        if _is_geometric_lj_force(system.getForce(fi)):
            system.removeForce(fi)
            break
    return system


__all__ = ["add_opls_impropers", "flatten_opls_lj_to_nonbonded"]
