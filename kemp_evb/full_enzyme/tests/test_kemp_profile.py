"""Offline reconstruction of recorded full-system mapping energies."""
import json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from kemp_profile import reconstruct,write_outputs,fit_target


class ProfileTests(unittest.TestCase):
    def test_recorded_mapping_profile_keeps_calibration_and_coverage_limits(self):
        raw=json.loads((Path(__file__).resolve().parents[1]/'demo_results/full_enzyme_engineering/replica_101_raw.json').read_text())
        analysis,summary=reconstruct(raw,10.,0.,6,1)
        self.assertFalse(summary['physical_barrier_available'])
        self.assertFalse(summary['calibrated'])
        self.assertAlmostEqual(summary['total_product_offset_kcal_mol'],44.5)
        with tempfile.TemporaryDirectory() as directory:
            write_outputs(analysis,summary,directory)
            svg=(Path(directory)/'profile.svg').read_text()
            self.assertIn('UNCALIBRATED TRIAL PROFILE',svg)
            self.assertNotIn('nan',svg.lower())
            self.assertTrue((Path(directory)/'profile.csv').exists())
            self.assertEqual(json.loads((Path(directory)/'summary.json').read_text())['internal_missing_bins'],summary['internal_missing_bins'])

    def test_target_fit_matches_declared_target_without_certifying_model(self):
        import numpy as np
        rng=np.random.default_rng(72);windows=[]
        for index,lam in enumerate(np.linspace(1,0,21)):
            samples=[]
            for step,q in enumerate(rng.normal(1-2*lam,(.008314462618*300/100)**.5,300)):
                samples.append({'step':step,'E_R_kj_mol':50*(q+1)**2,'E_P_kj_mol':50*(q-1)**2})
            windows.append({'index':index,'lambda_R':float(lam),'samples':samples})
        raw={'schema':'kemp_full_mapping_records_v1','metadata':{'morse_P_minus_R_correction_kcal_mol':0.},
             'config':{'temperature_kelvin':300.},'cursor':{'phase':'complete'},'windows':windows}
        analysis,summary=fit_target(raw,2.,5.,18,10)
        self.assertAlmostEqual(summary['conditional_barrier_estimate_kcal_mol'],5.,places=4)
        self.assertFalse(summary['physical_barrier_available'])
        self.assertFalse(summary['calibration_validated'])
        self.assertEqual(summary['reference'],'reactant')
        with tempfile.TemporaryDirectory() as directory:
            write_outputs(analysis,summary,directory)
            self.assertIn('FITTED TO CHOSEN TARGET',(Path(directory)/'profile.svg').read_text())

    def test_diagnostics_and_nonfinite_parameters_are_rejected(self):
        with self.assertRaises(ValueError):reconstruct({'production_complete':True},10.,0.,40,20)
        with self.assertRaises(ValueError):reconstruct({'schema':'kemp_full_mapping_records_v1'},float('nan'),0.,40,20)

if __name__=='__main__':unittest.main()
