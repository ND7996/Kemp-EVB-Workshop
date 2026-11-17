# Run the full-enzyme Kemp calculation

The primary calculation now uses the complete 6,330-atom Des27.7-IN system. The five-atom fixed-path illustration is separate and is not a substitute for enzyme sampling.

Use the [public installation/run instructions](https://nd7996.github.io/Kemp-EVB-Workshop/#linux) or the [full-enzyme guide](../full_enzyme/README.md).

From the repository root:

```bash
conda env update -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python -m pip install --no-deps --no-build-isolation ./EVBOpenMM
python kemp_evb/calculation/scripts/run_kemp_evb.py
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode diagnose
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode profile --bins 18 --open
```

Use `conda env create` for a new environment. Existing outputs are preserved; use `--resume` with unchanged settings to continue, or choose a new `--output-dir`.

The full workflow prepares the matching states, minimizes/heats/equilibrates, collects mapping-window energies, writes trajectories and reconstructs the profile. Short settings do not establish convergence; boundary electrostatics and transferred parameters differ from native Q6. The result is not a validated reproduction of the published activation barrier.

The older `run_short_evb_demo.py` command now forwards to this full-system workflow. For the explicitly separate five-atom illustration only, use `run_five_atom_example.py`.
