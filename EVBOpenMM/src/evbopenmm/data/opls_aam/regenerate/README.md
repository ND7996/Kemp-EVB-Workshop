# Regenerating `opls-aam.xml`

`../opls-aam.xml` (the OPLS-AA/M force field behind `force_field="opls-aam"`) is generated from the
Jorgensen-group CHARMM-format OPLS-AA/M parameter files with a **fixed** converter.  The upstream
converter (leelasd/OPLS-AAM `processCharmmForceField.py`, itself derived from OpenMM/devtools)
produced a file that passed `createSystem` but was silently wrong in six ways; `build_opls_aam.py`
here carries the fixes so the XML is reproducible instead of hand-edited.

## Regenerate

```bash
python3 build_opls_aam.py top_opls_aam.inp par_opls_aam.inp > ../opls-aam.xml.base
```

`top_opls_aam.inp` / `par_opls_aam.inp` are the Yale files (http://zarbi.chem.yale.edu/oplsaam.html)
with CRLF stripped.

## The six fixes baked into `build_opls_aam.py`

1. **Python 3** — `len(fields)//2`, `size//2` (was Py2 `/`).
2. **Singular section headers** — the par file uses `BOND`/`ANGLE`/`DIHEDRAL`; the parser now
   normalizes them to the plural names it expects (else all bonds/angles/dihedrals silently drop).
3. **Case-insensitive `END`** — the top file ends with uppercase `END`.
4. **Periodic impropers** — OPLS impropers are periodic (amplitude, n=2, phase 180), emitted as
   `<Improper periodicity1 phase1 k1>` inside the `PeriodicTorsionForce` (not harmonic in a
   `CustomTorsionForce`, which is never applied). OpenMM's matcher still won't apply them, so
   `evbopenmm.forces.add_opls_impropers` reads them and applies them explicitly (central atom = pos 2).
5. **Real LJ, not placeholders** — the `<Script>` now assigns each atom its self sigma/epsilon from
   the class LJ-matrix diagonal, and reassigns the 1-4 LJ exceptions from the sigma14/epsilon14
   matrices (both were left at the sigma=1/epsilon=1 placeholder).
6. **PME method map** — the LJ `CustomNonbondedForce` maps a periodic Coulomb method (Ewald/PME/
   LJPME) to `CutoffPeriodic` (copying `PME` verbatim raised "Illegal value for nonbonded method").

Verified: the regenerated base passes the six-point checklist in `tests/test_opls_forcefield.py`.

## Ca²⁺ add-on (not in the CHARMM inputs)

The shipped `../opls-aam.xml` = the regenerated base **plus** a Ca²⁺ ion (OPLS-AA/M has none):
add `<Type name="900" class="CAL" element="Ca" mass="40.08"/>` to `<AtomTypes>`, a
`<Atom type="900" charge="2.0" sigma="1.0" epsilon="1.0"/>` to `<NonbondedForce>`, and set the
Script's `EXTRA_LJ = {'900': (0.243657, 0.502080)}` (CHARMM CAL: Rmin/2 1.367 Å, eps 0.120 kcal/mol).
The `CA` residue template lives in `../opls-aam-patches.xml`.
