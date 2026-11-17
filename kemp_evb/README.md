# Kemp EVB workshop

Open [`index.html`](index.html) to start. It explains EVB and spherical boundaries, shows the Kemp reaction scheme, then compares the paper’s Q6 result with the short OpenMM example before the α/H₁₂ parameter examples.

## Q6 inputs and outputs

- Inputs: starting coordinates (`.pdb.gz`), topology and force-field terms (`.top.gz`), and the two-state reaction definition and atom map (`.fep`). A full dynamics/analysis run also needs Qdyn stage/window control files (`.inp`) and a QFEP analysis control file.
- Outputs: Qdyn produces trajectory/restart data and energies for each mapping window (for example, `fep_100.en`). QFEP reads those window energies and its analysis settings, then writes `qfep.out` and the free-energy profile.
- The workshop’s bundled Q6 inputs are the coordinates, topology and reaction definition; they are not a complete production Q6 run setup.

## Short OpenMM example inputs and outputs

- Inputs: the example checks the bundled Q6 coordinate, topology and FEP files for a consistent atom map. It uses five mapped atom positions and selected reactive terms, plus the `--coupling` and `--offset` settings. The scripts are `run_short_evb_demo.py`, `kemp_inputs.py` and `evb_teaching.py`; `environment.yml` lists the dependencies.
- Outputs: a new run folder with `teaching_evb.html` (interactive potential-energy illustration) and `teaching_summary.json`. This fixed-path teaching example does not produce Qdyn trajectories, QFEP `.en` files or a sampled free-energy profile.
- The α/H₁₂ parameter scan on the page reuses saved energies and does not start another run.
