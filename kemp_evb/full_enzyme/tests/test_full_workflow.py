"""Numerical minimization and transactional restart tests, without long MD."""
import json
import sqlite3
import sys
import tempfile
from pathlib import Path
import unittest
import numpy as np
import openmm as mm
from openmm import app, unit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from full_workflow import DEFAULTS, Store, bounded_minimize, diagnose, finite_state, load_config


def toy():
    topology=app.Topology();chain=topology.addChain('A');res=topology.addResidue('TST',chain)
    a=topology.addAtom('C',app.element.carbon,res);b=topology.addAtom('O',app.element.oxygen,res);topology.addBond(a,b)
    system=mm.System();system.addParticle(0);system.addParticle(16)
    force=mm.HarmonicBondForce();force.addBond(0,1,.1,1000.);system.addForce(force)
    sim=app.Simulation(topology,system,mm.VerletIntegrator(.0001),mm.Platform.getPlatformByName('Reference'))
    sim.context.setPositions([[0,0,0],[.2,0,0]])
    return sim


class WorkflowTests(unittest.TestCase):
    def test_rejects_bad_sampling_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'config.json'
            for bad in ({'sample_interval':True},{'lambda_schedule':[1.,.5,.5,0.]},{'timestep_fs':float('nan')},{'sample_interval':2000},{'platform_typo':'CPU'}):
                path.write_text(json.dumps(bad))
                with self.assertRaises(ValueError):load_config(path)
            path.write_text('{}');self.assertEqual(load_config(path)['seed'],101)

    def test_bounded_minimization_reduces_energy_keeps_frozen_atom(self):
        sim=toy()
        with tempfile.TemporaryDirectory() as directory:
            result=bounded_minimize(sim,{**DEFAULTS,'relaxation_steps':30,'relaxation_lbfgs_iterations':0,'relaxation_force_tolerance_kj_mol_nm':1.},Path(directory))
            self.assertLess(result['final_energy_kj_mol'],result['initial_energy_kj_mol'])
            self.assertFalse(result['converged_force_tolerance'])
            self.assertEqual(result['stop_reason'],'iteration_limit')
            np.testing.assert_array_equal(finite_state(sim)[2][0],[0.,0.,0.])
            rows=(Path(directory)/'relaxation.csv').read_text().splitlines()[1:]
            energies=[float(row.split(',')[1]) for row in rows]
            self.assertTrue(all(a>b for a,b in zip(energies,energies[1:])))

    def test_presets_are_valid_and_long_template_has_documented_duration(self):
        config_root=Path(__file__).resolve().parents[1]/'configs'
        for name in ('full_enzyme_short.json','full_enzyme_smoke.json','full_enzyme_production.json'):
            cfg=load_config(config_root/name)
            self.assertEqual(cfg['lambda_schedule'][0],1.)
            self.assertEqual(cfg['lambda_schedule'][-1],0.)
        self.assertAlmostEqual(len(cfg['lambda_schedule'])*cfg['production_steps_per_window']*cfg['timestep_fs']/1e6,10.2)

    def test_fep_diagnostics_constant_gap_does_not_certify_convergence(self):
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)
            sample={'E_R_kj_mol':0.,'E_P_kj_mol':4.}
            raw={'config':{'temperature_kelvin':300.},'cursor':{'phase':'complete'},'windows':[
                {'index':0,'lambda_R':1.,'samples':[sample,sample]},
                {'index':1,'lambda_R':0.,'samples':[sample,sample]}]}
            (out/'replica_101_raw.json').write_text(json.dumps(raw))
            result=diagnose(out,101);edge=result['edges'][0]
            self.assertAlmostEqual(edge['forward_deltaG_kj_mol'],4.)
            self.assertAlmostEqual(edge['reverse_deltaG_kj_mol'],4.)
            self.assertFalse(result['convergence_established'])
            self.assertFalse(result['calibrated'])

    def test_samples_and_restart_cursor_are_atomic(self):
        sim=toy()
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'run.sqlite';store=Store(path,'test',False)
            sim.step(3);expected=finite_state(sim)[2].copy()
            cursor={'phase':'production','window':0,'production_done':3}
            row=(0,3,3,1.,4.,5.,4.)
            store.save(sim,cursor,[row]);checkpoint=store.get('checkpoint')
            sim.step(2)
            with self.assertRaises(sqlite3.IntegrityError):store.save(sim,{'phase':'complete'},[row])
            self.assertEqual(store.cursor(),cursor);self.assertEqual(store.get('checkpoint'),checkpoint)
            store.db.close();resumed=Store(path,'test',True);resumed.restore(sim)
            self.assertEqual(sim.currentStep,3);np.testing.assert_array_equal(finite_state(sim)[2],expected)
            self.assertEqual(resumed.db.execute('SELECT COUNT(*) FROM samples').fetchone()[0],1)
            resumed.db.close()
            with self.assertRaises(ValueError):Store(path,'different model',True)
            with self.assertRaises(ValueError):Store(path,'test',False)

if __name__=='__main__':unittest.main()
