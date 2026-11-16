# Short Kemp EVB calculation

This folder contains a runnable, reduced five-atom Kemp reaction calculation for the hands-on. It evaluates an 81-point, user-defined geometric path with OpenMM and checks the lower two-state energy against the pinned EVBOpenMM implementation. The input geometry and Q6 force-field terms come from the published Des27.7-IN files cited below.

This is a software and EVB demonstration. It is **not** molecular dynamics, a full-enzyme calculation, an optimized reaction path, or a free-energy/barrier estimate. Protein, solvent, most bonded/nonbonded terms, entropy, and path optimization are omitted. The H12=10 kcal/mol and extra product offset=0 kcal/mol defaults are teaching values, not fitted parameters.

## Files included here

- `scripts/run_short_evb_demo.py`: command-line entry point and input check.
- `scripts/evb_teaching.py`, `scripts/kemp_inputs.py`: reduced OpenMM model and source-input reader.
- `inputs/q6/Des27_7_IN/Des27_7_IN.fep`, `.pdb.gz`, `.top.gz`: the matched reaction mapping, coordinates and topology needed by these scripts.
- `environment.yml`: Python, OpenMM, NumPy and the pinned EVBOpenMM dependency.
- `cluster/short_evb_demo.slurm`: ready-to-adapt Slurm job for the short calculation.
- `cluster/full_enzyme_short.slurm`, `full_enzyme_production.slurm`: templates for the distinct full-system backend in `EVB_MHETase/Kemp_EVB`.

All calculation-specific inputs and scripts for the short demo are in this GitHub repository. The EVBOpenMM library is installed as a pinned public dependency; it is not duplicated in this repository.

## Run locally on Linux

From the root of the cloned `Kemp-EVB-Workshop` repository:

```bash
conda env create -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python kemp_evb/calculation/scripts/run_short_evb_demo.py --check
python kemp_evb/calculation/scripts/run_short_evb_demo.py
```

The last command creates a fresh folder under `kemp_evb/calculation/runs/` containing `teaching_evb.csv`, `teaching_summary.json` and `teaching_evb.html`. Open the HTML file in a browser. To choose another teaching setting, pass `--coupling 12 --offset 2`; use a new output path if you supply `--output`.

## Run the short demo on a Slurm cluster

Use the login node only to clone/pull the workshop, create the environment, and install dependencies. Run the calculation on a compute node through Slurm. If your cluster requires a module, load the site’s Conda/Miniforge module before creating the environment, or identify its base directory for `CONDA_BASE`.

From the workshop repository root on the cluster:

```bash
conda env create -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python kemp_evb/calculation/scripts/run_short_evb_demo.py --check
CONDA_BASE="$(conda info --base)" sbatch kemp_evb/calculation/cluster/short_evb_demo.slurm
```

`CONDA_BASE` is passed into the job so the non-interactive shell can activate the environment. If your cluster does not expose `conda` in the login shell, set it to the installation path shown by your cluster documentation. The short demo uses OpenMM’s CPU `Reference` platform; no GPU is required. Check the job with `squeue -u "$USER"`; output and errors go to `kemp-evb-JOBID.out` and `.err`, and calculation results go under `kemp_evb/calculation/runs/`.

## Full-enzyme OpenMM/EVBOpenMM route

The full-system workflow is maintained in the separate [EVB_MHETase/Kemp_EVB repository](https://github.com/ND7996/EVB_MHETase/tree/main/Kemp_EVB). It uses a changed finite-droplet model and is not an exact Q6 reproduction. Clone both source repositories on the login node and install the compatible EVBOpenMM revision into the backend environment. Adjust the working-directory path if your cluster stores projects somewhere other than `$HOME/projects`:

```bash
mkdir -p "$HOME/projects"
cd "$HOME/projects"
git clone https://github.com/ND7996/EVB_MHETase.git
git clone https://github.com/Martin-Floor/EVBOpenMM.git
git -C EVBOpenMM checkout dedf301eb9b0021a4be200327cea53f5b7f034c6
cd EVB_MHETase/Kemp_EVB
conda env create -f environment.yml
conda activate kemp-openmm
python -m pip install -e "$HOME/projects/EVBOpenMM"
```

Set the paths to your actual Conda and backend checkout. From this workshop repository root, submit the short full-system mapping workflow:

```bash
CONDA_BASE="/path/to/miniforge3" \
KEMP_EVB_DIR="/path/to/EVB_MHETase/Kemp_EVB" \
sbatch kemp_evb/calculation/cluster/full_enzyme_short.slurm
```

That preset has the full 6,330-atom system, 21 mapping windows and one short replica. It exercises preparation, relaxation, heating, equilibration and mapping-data collection; it is not a convergence study or a validated physical barrier.

For longer production, create a persistent run directory on cluster storage and submit:

```bash
mkdir -p /path/to/persistent/kemp-full-production
CONDA_BASE="/path/to/miniforge3" \
KEMP_EVB_DIR="/path/to/EVB_MHETase/Kemp_EVB" \
KEMP_RUN_DIR="/path/to/persistent/kemp-full-production" \
sbatch kemp_evb/calculation/cluster/full_enzyme_production.slurm
```

The production config samples 51 windows × 800,000 steps at 0.25 fs (10.2 ns total mapped production, plus preparation/equilibration). The example job runs up to 11 hours, saves progress, and should be resubmitted to resume if it reaches the scheduler limit:

```bash
CONDA_BASE="/path/to/miniforge3" \
KEMP_EVB_DIR="/path/to/EVB_MHETase/Kemp_EVB" \
KEMP_RUN_DIR="/path/to/persistent/kemp-full-production" \
RESUME=1 sbatch kemp_evb/calculation/cluster/full_enzyme_production.slurm
```

Check the cluster’s wall-time, memory, CPU, module, and scratch-storage rules and edit the `#SBATCH` lines in the templates when needed. The production length does not by itself establish convergence or calibration. Keep run data out of GitHub; upload source/configuration and compact plots/tables, not trajectories or checkpoints.

## Source and attribution

The selected `Des27_7_IN` structure, topology and FEP mapping are from the Kemp EVB deposited input dataset, [Zenodo DOI 10.5281/zenodo.14563437](https://doi.org/10.5281/zenodo.14563437), associated with Listov et al., “Complete computational design of high-efficiency Kemp elimination enzymes,” *Nature* (2025), [10.1038/s41586-025-09136-2](https://doi.org/10.1038/s41586-025-09136-2). The compact calculation preserves the source terms used by the teaching script; it does not reproduce the paper’s full Hamiltonian or reported barrier.
