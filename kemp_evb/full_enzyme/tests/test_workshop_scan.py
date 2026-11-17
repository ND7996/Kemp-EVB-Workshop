"""Workshop target matching must remain labelled and visibly non-unique."""
import hashlib,json,sys,tempfile,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from workshop_scan import scan,write_scan


def synthetic_records():
    # Independent samples from harmonic mapping ensembles; not enzyme data.
    rng=np.random.default_rng(44);windows=[]
    for i,lam in enumerate(np.linspace(1,0,21)):
        samples=[]
        for j,g in enumerate(rng.normal(400*(1-2*lam),(.596*800)**.5,80)):
            samples.append({'step':j,'E_R_kj_mol':float((g+400)**2/1600*4.184),'E_P_kj_mol':float((g-400)**2/1600*4.184)})
        windows.append({'index':i,'lambda_R':float(lam),'samples':samples})
    return {'schema':'kemp_full_mapping_records_v1','metadata':{'morse_P_minus_R_correction_kcal_mol':44.5},
        'config':{'temperature_kelvin':300.},'cursor':{'phase':'complete'},'windows':windows}


class WorkshopScanTests(unittest.TestCase):
    def test_multiple_target_matches_are_not_unique_calibration(self):
        raw=synthetic_records();before=json.dumps(raw,sort_keys=True)
        result,analysis,summary=scan(raw,[40.,60.,80.,100.,120.,140.],[-120.,-110.],13.2,18,8)
        self.assertEqual(len(result['grid_candidates']),12)
        self.assertEqual(len(result['target_matches']),2)
        for match in result['target_matches']:self.assertAlmostEqual(match['barrier_kcal_mol'],13.2,places=4)
        self.assertNotAlmostEqual(result['target_matches'][0]['reaction_dg_kcal_mol'],result['target_matches'][1]['reaction_dg_kcal_mol'])
        self.assertEqual(result['best']['native_alpha_kcal_mol'],-110.)
        for row in result['target_matches']+[result['baseline']]:
            peak=row['profile_peak'];i=peak['bin_index']
            self.assertAlmostEqual(peak['raw_gap_kcal_mol'],row['x'][i])
            self.assertAlmostEqual(peak['offset_adjusted_gap_kcal_mol'],row['x'][i]-row['total_product_offset_kcal_mol'])
            self.assertAlmostEqual(peak['free_energy_relative_reactant_kcal_mol'],row['barrier_kcal_mol'])
            self.assertEqual(peak['raw_peak_bin_contains_zero'],peak['raw_bin_lower_kcal_mol']<=0<=peak['raw_bin_upper_kcal_mol'])
        self.assertFalse(result['calibrated']);self.assertTrue(summary['workshop_only'])
        self.assertFalse(summary['physical_barrier_available'])
        self.assertEqual(summary['native_qfep']['coupling_kcal_mol'],100.4)
        self.assertEqual(json.dumps(raw,sort_keys=True),before)
        with tempfile.TemporaryDirectory() as directory:
            write_scan(result,analysis,summary,directory)
            page=(Path(directory)/'scan.html').read_text()
            self.assertIn('NOT A PREDICTION',page)
            self.assertNotIn('src="https:',page)
            self.assertIn('legend-archived',page)
            self.assertIn('Centre selected profile maximum',page)
            self.assertIn('display shift only',page)
            self.assertIn('zero-line',page)
            self.assertTrue((Path(directory)/'combinations.csv').exists())
            self.assertIn('WORKSHOP TARGET MATCH',(Path(directory)/'selected_example/profile.svg').read_text())

    def test_invalid_and_unresolved_inputs_do_not_create_a_barrier(self):
        raw=synthetic_records()
        for couplings,offsets in (([],[-110.]),([0.],[-110.]),([10.],[float('nan')])):
            with self.assertRaises(ValueError):scan(raw,couplings,offsets)
        with self.assertRaises(ValueError):scan(raw,[80.,100.],[-110.],13.2,18,10000)

if __name__=='__main__':unittest.main()
