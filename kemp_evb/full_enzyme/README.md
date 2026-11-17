# Full-enzyme Kemp EVBOpenMM workflow

This restores the full-system OpenMM adaptation from EVB_MHETase commit `e933495` into the standalone Kemp repository. It uses the bundled top-level EVBOpenMM source. The default Des27.7-IN model contains **6,330 atoms, 771 water molecules, 22 mapped reactive atoms and three forming/breaking Morse pairs**. It is not the five-atom illustration.

The imported Q/OPLS coefficients, charges and bonded state changes come from the deposited Q6 inputs. This adaptation uses a **changed finite-droplet Hamiltonian**: direct nonperiodic electrostatics, a water-oxygen wall and frozen outer protein. It does not reproduce Q6's LRF/SCAAS boundary treatment or all native restraints. Source coupling/offset values for profile analysis are taken from the separate GNCA4-WT QFEP reference archive; their transfer to Des27.7-IN is unvalidated. No exact Q6 or published-barrier reproduction is claimed.

## Install

From the repository root:

```bash
conda env update -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python -m pip install --no-deps --no-build-isolation ./EVBOpenMM
```

Use `conda env create` instead of `update` for a new environment. All paths are relative; no private-repository access is required.

## Run the complete short full-enzyme workflow

```bash
python kemp_evb/calculation/scripts/run_kemp_evb.py
```

This builds matching states, performs bounded descent/L-BFGS relaxation, heats and equilibrates the reactant endpoint, then samples **21 mapping windows**. It writes `topology.pdb`, segmented DCD/energy traces, `replica_101_raw.json`, diagnostics metadata and transactional restarts under `kemp_evb/full_enzyme/runs/full_enzyme_short/`.

The short preset is a workflow test, not converged sampling: timestep 0.25 fs, 200 heating steps per temperature, 1,000 endpoint-equilibration steps, 200 equilibration + 1,000 production steps/window. Twenty-one windows × 1,000 × 0.25 fs is only **5.25 ps total production**. Benchmark on your hardware; no overnight completion time is guaranteed.

## Diagnose and view the new profile

```bash
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode diagnose
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode profile --bins 18 --open
```

The profile command reads your full-system diabatic-energy samples and reconstructs an EVB free-energy estimate. It does not draw the prescribed five-atom path. Blue denotes the computed OpenMM curve; the grey line/band is the published **Q6 computed barrier mean ± replica SD** for the selected case, not a complete reference curve. Missing bins remain gaps; a disconnected profile does not get an invented barrier. Coarser binning is not additional sampling.

## Reconstruct the previously saved full-enzyme result immediately

```bash
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode profile --records kemp_evb/full_enzyme/demo_results/source_parameter_comparison/replica_101_raw.json --output-dir kemp_evb/full_enzyme/runs/recorded_comparison --bins 18 --open
```

The original raw SHA256 matches the archived plot manifest. This calculation reproduces the saved **9.572629 kcal/mol parameter-conditional estimate**, versus the paper's 13.2 ± 0.6 kcal/mol Q6 result. It is not a new MD run, a fit to that barrier, or evidence of convergence. Different boundary, parameter-transfer and sampling choices all limit the comparison; longer sampling alone does not establish equivalence.

## Separate stages and resuming

The one-command workflow is recommended. For separate stages use a fresh output directory and the same config throughout:

```bash
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode prepare
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode relax --resume
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode equilibrate --resume
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode production --resume
```

Do not repeat `prepare` over an existing run. A completed short run is resumed with:

```bash
python kemp_evb/calculation/scripts/run_kemp_evb.py --resume --max-runtime-seconds 10800
```

Omit `--resume` for a fresh timed run. The time budget is checked between chunks, not a hard scheduler deadline. `run.sqlite` commits checkpoint state and sample cursors together. Preserve it and keep code, config, inputs, OpenMM and platform consistent. Do not concatenate trajectory segments blindly after a crash.

## Longer runs and independent replicas

Edit a copied `configs/full_enzyme_short.json` and select it with `--config <copy.json>`. Increase `equilibration_steps`, `window_equilibration_steps` and `production_steps_per_window`; refine `lambda_schedule` according to measured overlap. Time in ps = steps × timestep_fs / 1000. `full_enzyme_production.json` is a longer **unvalidated template**, not a guarantee of reproduction or convergence. Give each independent seed a separate output directory. Never change the config to resume an existing run.

Validate the boundary/force-field transfer, stability, target calibration, effective sampling and replica uncertainty before interpreting a physical barrier. The default QFEP-style postprocessor applies saved H12=100.4 and alpha=−110 kcal/mol plus the documented +44.5 kcal/mol Morse energy-zero correction. These are not established Des27.7-IN calibration parameters. Do not tune them simply to manufacture the paper's barrier.

## Tests

```bash
python -m unittest discover -s kemp_evb/full_enzyme/tests -v
```

Raw-input and numerical tests cover coefficient conversion, state charges/constraints, minimization, transactional restart and profile reconstruction. They do not prove native Q6 Hamiltonian parity. `RESTORATION_PROVENANCE.json` records the original source files and hashes; `docs/VALIDATION.md` records the newly executed checks.
