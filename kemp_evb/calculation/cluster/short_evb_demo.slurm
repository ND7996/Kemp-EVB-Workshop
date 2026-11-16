#!/bin/bash
#SBATCH --job-name=kemp-evb-demo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=kemp-evb-%j.out
#SBATCH --error=kemp-evb-%j.err
set -euo pipefail

# Supply CONDA_BASE when submitting; see guides/CLUSTER_QUICKSTART.md.
: "${CONDA_BASE:?Set CONDA_BASE to your Miniforge/Conda installation directory}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate kemp-short
cd "$SLURM_SUBMIT_DIR"

srun python kemp_evb/calculation/scripts/run_short_evb_demo.py \
  --case Des27_7_IN \
  --coupling 10 \
  --offset 0 \
  --output "kemp_evb/calculation/runs/evb-demo-${SLURM_JOB_ID}"
