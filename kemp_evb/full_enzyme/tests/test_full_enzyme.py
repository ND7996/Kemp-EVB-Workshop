"""Import contracts for the explicitly different full-enzyme adaptation."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from kemp_inputs import CASES, case_paths, read, sections
from q_topology import parse_topology
from full_enzyme import state_system, patched_terms


class FullInputTests(unittest.TestCase):
    def test_complete_native_payloads(self):
        for case in CASES:
            d=parse_topology(read(case_paths(case)[2]))
            self.assertEqual(len(d['coordinates_angstrom']),d['atoms'])
            self.assertEqual(len(d['charges_e']),d['atoms'])
            self.assertAlmostEqual(sum(d['charges_e']),round(sum(d['charges_e'])),places=3)
            self.assertEqual(d['vdw_combination_rule'],1)
            self.assertEqual(d['coulomb_14_scale'],.5)
            excluded={tuple(p) for key in ('short_pairs','long_pairs') for p in d['pairs']['excluded'][key]}
            for o in range(d['solute_atoms']+1,d['atoms']+1,3):
                self.assertTrue({(o,o+1),(o,o+2),(o+1,o+2)} <= excluded)

    def test_imported_lj_matches_native_pair_coefficients(self):
        # Native tables contain sqrt(Aii), sqrt(Bii), not Aii and Bii.
        # Test energy at a finite separation, independently of sigma/epsilon conversion.
        import openmm
        from openmm import unit
        d=parse_topology(read(case_paths(CASES[0])[2]));f=sections(read(case_paths(CASES[0])[0]))
        system,*_=state_system(d,f,0)
        nb=next(force for force in system.getForces() if isinstance(force,openmm.NonbondedForce))
        oxygen=d['solute_atoms'];code=d['atom_type_codes'][oxygen]-1
        a=d['native_LJ_tables']['sqrt (Aii) normal:'][code]
        b=d['native_LJ_tables']['sqrt (Bii) normal:'][code]
        _,sigma,epsilon=nb.getParticleParameters(oxygen)
        sigma=sigma.value_in_unit(unit.nanometer);epsilon=epsilon.value_in_unit(unit.kilojoule_per_mole)
        for distance_angstrom in (2.8,3.2,4.0):
            native=(a*a/distance_angstrom**12-b*b/distance_angstrom**6)*4.184
            ratio=sigma/(distance_angstrom*.1)
            converted=4*epsilon*(ratio**12-ratio**6)
            self.assertAlmostEqual(native,converted,places=9)
        self.assertGreater(epsilon,.1)

    def test_state_charge_and_reactive_morse_convention(self):
        d=parse_topology(read(case_paths(CASES[0])[2]));f=sections(read(case_paths(CASES[0])[0]))
        states=[state_system(d,f,i) for i in (0,1)]
        self.assertEqual([len(s[1]) for s in states],[2,1])
        self.assertAlmostEqual(sum(v[0] for v in states[0][1].values())-sum(v[0] for v in states[1][1].values()),44.5)
        self.assertEqual([len(s[2]) for s in states],[1,0])
        self.assertEqual([s[0].getNumConstraints() for s in states],[2313,2313])
        import openmm
        from openmm import unit
        totals=[]
        for s,*_ in states:
            nb=next(force for force in s.getForces() if isinstance(force,openmm.NonbondedForce))
            totals.append(sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(s.getNumParticles())))
        self.assertAlmostEqual(totals[0],totals[1],places=5)
        self.assertAlmostEqual(totals[0],-7.,places=3)

if __name__=='__main__':unittest.main()
