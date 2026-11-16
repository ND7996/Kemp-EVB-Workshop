# Cluster quick start

The workshop includes an actual short calculation, not only saved plots. The first route below runs a reduced five-atom Kemp EVB scan using the inputs bundled in the workshop repository. The second route runs the separate full-enzyme OpenMM/EVBOpenMM workflow. The full workflow needs its own backend checkout.

## 1. Prepare once on the cluster login node

Use your cluster’s documented Conda/Miniforge module if needed. From the workshop repository root:

```bash
conda env create -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python kemp_evb/calculation/scripts/run_short_evb_demo.py --check
```

Environment creation installs Python 3.12, OpenMM 8.6.1, NumPy and EVBOpenMM at the workshop’s pinned revision. If compute nodes cannot access the internet, this still works because installation occurs on the login node. If login nodes also block downloads, install from the cluster’s approved package cache or stage the pinned dependency through an approved route before creating the environment.

## 2. Submit the short live EVB demo

From the repository root, pass the Conda base path and submit the included Slurm script:

```bash
CONDA_BASE="$(conda info --base)" sbatch kemp_evb/calculation/cluster/short_evb_demo.slurm
squeue -u "$USER"
```

The job uses one CPU node, two CPU threads and OpenMM’s CPU `Reference` platform. It evaluates 81 geometries and compares OpenMM’s lower two-state energy with EVBOpenMM at each point. It is a fast, reduced teaching scan; it does not calculate an enzyme activation free energy. Results are saved under `kemp_evb/calculation/runs/evb-demo-JOBID/`.

If Conda is provided as a module and the shell does not know `conda info --base`, use the module’s documented installation path instead:

```bash
CONDA_BASE="/cluster/path/to/miniforge3" sbatch kemp_evb/calculation/cluster/short_evb_demo.slurm
```

## 3. Submit a short full-enzyme workflow

On the login node, clone the backend and its EVBOpenMM dependency, then install the matching package into the backend environment. Replace `$HOME/projects` with your project directory if needed:

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

Point `KEMP_EVB_DIR` to the `Kemp_EVB` directory in that checkout. Then submit:

```bash
CONDA_BASE="/cluster/path/to/miniforge3" \
KEMP_EVB_DIR="/cluster/path/to/EVB_MHETase/Kemp_EVB" \
sbatch kemp_evb/calculation/cluster/full_enzyme_short.slurm
```

This runs `configs/full_enzyme_short.json`: 6,330 atoms, 21 windows and one short replica, including preparation, relaxation, heating, equilibration and production mapping. It tests the end-to-end workflow but does not establish profile convergence or physical calibration.

## 4. Submit longer mapping production

Use persistent cluster storage (not a temporary node-local directory) because the simulation writes restart and trajectory files. Make an output directory, then submit:

```bash
mkdir -p /cluster/scratch/yourname/kemp-production
CONDA_BASE="/cluster/path/to/miniforge3" \
KEMP_EVB_DIR="/cluster/path/to/EVB_MHETase/Kemp_EVB" \
KEMP_RUN_DIR="/cluster/scratch/yourname/kemp-production" \
sbatch kemp_evb/calculation/cluster/full_enzyme_production.slurm
```

The production config contains 51 windows × 800,000 steps at 0.25 fs, or 10.2 ns of mapped production, in addition to preparation and equilibration. The template runs for at most 11 hours per job. If it stops at the time limit, submit the same job with `RESUME=1` and the same output directory:

```bash
CONDA_BASE="/cluster/path/to/miniforge3" \
KEMP_EVB_DIR="/cluster/path/to/EVB_MHETase/Kemp_EVB" \
KEMP_RUN_DIR="/cluster/scratch/yourname/kemp-production" \
RESUME=1 sbatch kemp_evb/calculation/cluster/full_enzyme_production.slurm
```

Use `squeue -u "$USER"` to check jobs and `scancel JOBID` to stop a job. Add your cluster’s required account, partition, constraint, module-load and storage settings to the `#SBATCH` header or job body. The examples request CPU resources; the current configurations select OpenMM CPU. Ask the cluster support team about local policy before increasing wall time or memory.

## What is and is not in the public workshop repository

The short five-atom calculation is self-contained here: Python scripts, matched FEP/PDB/topology inputs, environment specification, and Slurm job. It does not require downloading the enzyme backend source. The full-enzyme scripts, full protein system, and longer configurations are maintained in `EVB_MHETase/Kemp_EVB` and are fetched separately. Neither path is a validated reproduction of the paper’s free-energy profile; see the calculation README for scope and source details.
