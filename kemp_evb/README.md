# Kemp EVB workshop

Open the [live workshop website](https://nd7996.github.io/Kemp-EVB-Workshop/?v=portable-setup-9060d6f), or open [`index.html`](index.html) locally. It explains EVB and spherical boundaries, shows the Kemp reaction scheme, then shows a clearly labeled Q6 vs OpenMM comparison and its input/output files before the α/H₁₂ examples.

## Q6 inputs and outputs

- Inputs: starting coordinates (`.pdb.gz`), topology and force-field terms (`.top.gz`), and the two-state reaction definition and atom map (`.fep`). A full dynamics/analysis run also needs Qdyn stage/window control files (`.inp`) and a QFEP analysis control file.
- Outputs: Qdyn produces trajectory/restart data and energies for each mapping window (for example, `fep_100.en`). QFEP reads those window energies and its analysis settings, then writes `qfep.out` and the free-energy profile.
- The workshop’s bundled Q6 inputs are the coordinates, topology and reaction definition; they are not a complete production Q6 run setup.

## Short OpenMM example inputs and outputs

- Inputs: the example checks the bundled Q6 coordinate, topology and FEP files for a consistent atom map. It uses five mapped atom positions and selected reactive terms, plus the `--coupling` and `--offset` settings. The scripts are `run_short_evb_demo.py`, `kemp_inputs.py` and `evb_teaching.py`; `environment.yml` lists the dependencies.
- Outputs: a new run folder with `teaching_evb.html` (interactive potential-energy illustration) and `teaching_summary.json`. This fixed-path teaching example does not produce Qdyn trajectories, QFEP `.en` files or a sampled free-energy profile.
- The α/H₁₂ parameter scan on the page reuses saved energies and does not start another run.

## Installation

Use the [portable setup instructions](calculation/README.md): create the Conda environment, activate `kemp-short`, then install `calculation/requirements-evb.txt` with that environment’s Python. No personal source path is needed.
