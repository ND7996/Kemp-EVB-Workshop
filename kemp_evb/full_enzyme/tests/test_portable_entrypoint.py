"""Exercise the user-facing full-enzyme CLI and real saved energy record."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

WORKSHOP=Path(__file__).resolve().parents[2]
ENTRY=WORKSHOP/'calculation/scripts/run_kemp_evb.py'
FULL=WORKSHOP/'full_enzyme'

class PortableEntryPoint(unittest.TestCase):
    def test_saved_full_system_record_reconstructs_original_comparison(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'analysis'
            result=subprocess.run([sys.executable,str(ENTRY),'--mode=profile',f'--config={FULL/"configs/full_enzyme_short.json"}',
                '--records',str(FULL/'demo_results/source_parameter_comparison/replica_101_raw.json'),
                '--output-dir',str(out),'--bins','18'],capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            path=out/'profile_source_qfep_bins18_min20_skip0'
            summary=json.loads((path/'summary.json').read_text())
            self.assertEqual(summary['sampled_windows'],21)
            self.assertEqual(summary['sample_count'],1050)
            self.assertAlmostEqual(summary['conditional_barrier_estimate_kcal_mol'],9.572629012215984,places=8)
            self.assertFalse(summary['convergence_established'])
            self.assertIn('separate GNCA4-WT',summary['parameter_provenance_note'])
            self.assertIn('Grey: published Q6 barrier 13.2',(path/'profile.svg').read_text())

    def test_explicit_parameters_use_the_same_framewise_reconstruction(self):
        with tempfile.TemporaryDirectory() as td:
            result=subprocess.run([sys.executable,str(ENTRY),'--mode','profile','--coupling','100.4','--offset','-110',
                '--records',str(FULL/'demo_results/source_parameter_comparison/replica_101_raw.json'),
                '--output-dir',td,'--bins','18'],capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            summary=json.loads((Path(td)/'profile_trial_H100.4_offset-110_bins18_min20/summary.json').read_text())
            self.assertAlmostEqual(summary['conditional_barrier_estimate_kcal_mol'],9.572629012215984,places=8)
            self.assertFalse(summary['author_parameters_recovered'])
            self.assertFalse(summary['calibrated'])

    def test_legacy_command_now_advertises_the_full_workflow(self):
        path=WORKSHOP/'calculation/scripts/run_short_evb_demo.py'
        result=subprocess.run([sys.executable,str(path),'--help'],capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('production',result.stdout)
        self.assertIn('--resume',result.stdout)
        self.assertIn('--records',result.stdout)

if __name__=='__main__':unittest.main()
