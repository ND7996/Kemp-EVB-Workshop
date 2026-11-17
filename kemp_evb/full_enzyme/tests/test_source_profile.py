"""Archived-input interpretation and frame-wise exponential-average checks."""
import json,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from source_profile import read_qfep,reconstruct_source,frame_bin_free_energy


class SourceProfileTests(unittest.TestCase):
    def test_author_input_and_explicit_morse_conversion(self):
        p=read_qfep()
        self.assertEqual(p['coupling_kcal_mol'],100.4)
        self.assertEqual(p['native_alpha_kcal_mol'],-110.)
        self.assertEqual(p['coupling_form'],'constant')
        self.assertEqual(p['native_windows'],51)
        self.assertEqual(p['native_skip_samples_per_window'],4000)
        self.assertEqual(p['native_min_points_per_window_bin'],1000)
        self.assertAlmostEqual(p['native_alpha_kcal_mol']+44.5,-65.5)
        self.assertFalse(p['transfer_validated_for_current_model'])

    def test_frame_bin_exponential_average_and_gauge_invariance(self):
        kt=.6;values=[0.,kt*np.log(2.)]
        expected=3.-kt*np.log(.75)
        self.assertAlmostEqual(frame_bin_free_energy(values,kt,3.),expected)
        self.assertAlmostEqual(frame_bin_free_energy(np.array(values)+7.,kt,3.-7.),expected)
        with self.assertRaises(ValueError):frame_bin_free_energy([],kt,0.)

    def test_common_per_frame_energy_changes_leave_profile_unchanged(self):
        path=Path(__file__).resolve().parents[1]/'demo_results/full_enzyme_engineering/replica_101_raw.json'
        raw=json.loads(path.read_text());shifted=json.loads(path.read_text())
        for i,window in enumerate(shifted['windows']):
            for j,sample in enumerate(window['samples']):
                sample['E_R_kj_mol']+=100.*(i+j)
                sample['E_P_kj_mol']+=100.*(i+j)
        a,s=reconstruct_source(raw,6,1);b,t=reconstruct_source(shifted,6,1)
        va=[np.nan if v is None else v for v in a['profile']['state_free_energies_kcal_mol']['evb_ground']]
        vb=[np.nan if v is None else v for v in b['profile']['state_free_energies_kcal_mol']['evb_ground']]
        np.testing.assert_allclose(va,vb,atol=1e-8,equal_nan=True)
        self.assertFalse(s['calibrated']);self.assertTrue(s['author_parameters_recovered'])
        self.assertEqual(s['total_product_offset_kcal_mol'],-65.5)
        if s['internal_missing_bins']:self.assertIsNone(s['conditional_barrier_estimate_kcal_mol'])
        _,empty=reconstruct_source(raw,6,1000)
        self.assertEqual(empty['resolved_bins'],0)
        self.assertIsNone(empty['conditional_barrier_estimate_kcal_mol'])
        self.assertTrue(any('Fewer than three' in warning for warning in empty['warnings']))
        with self.assertRaises(ValueError):reconstruct_source(raw,6,1,4000)

    def test_frame_mixing_and_fep_normalization_match_independent_eigenvalues(self):
        # Symmetric energy increments2,4,6 give averaged forward/reverse deltaF=4kJ.
        # Each bin has one frame from each window; compute both estimates directly.
        source=read_qfep();source['coupling_kcal_mol']=1.;source['native_alpha_kcal_mol']=0.
        samples=[{'step':i,'E_R_kj_mol':0.,'E_P_kj_mol':p} for i,p in enumerate((2.,4.,6.))]
        raw={'schema':'kemp_full_mapping_records_v1','metadata':{'morse_P_minus_R_correction_kcal_mol':0.},
            'config':{'temperature_kelvin':300.},'cursor':{'phase':'complete'},'windows':[
                {'index':0,'lambda_R':1.,'samples':samples},{'index':1,'lambda_R':0.,'samples':samples}]}
        analysis,summary=reconstruct_source(raw,3,1,parameters=source)
        expected=[]
        for p in (6.,4.,2.):
            eigenvalue=np.linalg.eigvalsh(np.array([[0.,4.184],[4.184,p]]))[0]
            expected.append((eigenvalue+(eigenvalue-p+4.))/2./4.184)
        np.testing.assert_allclose(analysis['profile']['state_free_energies_kcal_mol']['evb_ground'],expected,atol=1e-10)
        self.assertEqual(summary['accepted_bin_counts'],[2,2,2])

    def test_semantic_fep_match_survives_different_pdb_indices(self):
        from kemp_inputs import ROOT,CASES,case_paths,read,sections
        reference=sections((ROOT/'inputs/protocol_reference/GNCA4-WT.fep').read_text())
        def normalize(data):
            inverse={int(row[1]):int(row[0]) for row in data['atoms']};normalized={}
            for key,rows in data.items():
                if key=='atoms':continue
                width={'change_bonds':2,'change_angles':3,'change_torsions':4,'change_impropers':4}.get(key)
                normalized[key]=sorted(tuple([str(inverse[int(x)]) for x in row[:width]]+row[width:]) if width else tuple(row) for row in rows)
            return normalized
        for case in CASES:self.assertEqual(normalize(reference),normalize(sections(read(case_paths(case)[0]))))

if __name__=='__main__':unittest.main()
