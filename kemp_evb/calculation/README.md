# Optional short EVB calculation

This Linux-friendly example uses OpenMM to evaluate an EVB model for five mapped atoms along a prescribed geometry path. It checks the lower EVB energy against EVBOpenMM. It is a software exercise, not molecular dynamics or a free-energy barrier calculation; the coupling and offset are teaching values.

## Inputs and outputs

The bundled Q6 reference inputs are a coordinate file (`.pdb.gz`), topology (`.top.gz`) and two-state reaction/atom-map file (`.fep`). The check command reads these to confirm that the topology, structure and reactive atom map agree. The reduced OpenMM calculation takes five mapped atom positions and selected reactive terms from these inputs, plus the `--coupling` and `--offset` values. It does not build a complete enzyme system from the Q6 topology.

The run command creates a new folder under `kemp_evb/calculation/runs/` with `teaching_evb.html` and `teaching_summary.json`. The HTML shows the evaluated potential-energy curves along the fixed teaching path; the JSON records the settings and run summary. The run does not create Qdyn window energies (`fep_*.en`), QFEP output (`qfep.out`) or a sampled free-energy profile.

## Included files

- scripts/run_short_evb_demo.py: input check and run command.
- scripts/evb_teaching.py and scripts/kemp_inputs.py: reduced model and source reader.
- the mapped Q6 coordinates, topology, and reaction file.
- environment.yml: Python, OpenMM, NumPy and installation tools from Conda.
- requirements-evb.txt: bundled Martin-Floor EVBOpenMM source snapshot, installed locally with pip.

## Install and run on Linux

Install Git and Conda or Miniforge. In a terminal:

~~~bash
git clone https://github.com/ND7996/Kemp-EVB-Workshop.git
cd Kemp-EVB-Workshop
conda env create -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python -m pip install --no-deps --no-build-isolation ./EVBOpenMM
python kemp_evb/calculation/scripts/run_short_evb_demo.py --check
python kemp_evb/calculation/scripts/run_short_evb_demo.py
~~~

No personal checkout path or `EVBOPENMM_SRC` setting is needed. Run these commands from the repository root. The pip step uses the active environment’s Python, installs the bundled EVBOpenMM source, and keeps Conda’s OpenMM/NumPy packages (`--no-deps`). No GitHub credentials are needed for the local package installation.

If the Conda step succeeds but the pip step fails, the environment remains usable: check that the top-level EVBOpenMM folder is present, return to the repository root and retry only the pip command. If `kemp-short` already exists, update it with `conda env update -f kemp_evb/calculation/environment.yml`, then activate it and run the pip step.

The final command saves teaching_evb.html and teaching_summary.json in a new folder under kemp_evb/calculation/runs/. Open the HTML file in a browser. A CPU is enough; no GPU is needed. To explore a different teaching setting, add --coupling 12 --offset 2 to the final command. These values are not fitted parameters.

The historical full-enzyme backend is available in a [fixed repository snapshot](https://github.com/ND7996/EVB_MHETase/tree/e933495/Kemp_EVB). The live workshop example here remains a reduced teaching calculation.
