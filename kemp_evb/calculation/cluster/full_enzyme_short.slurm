#!/bin/bash
#SBATCH --job-name=kemp-full-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=kemp-full-short-%j.out
#SBATCH --error=kemp-full-short-%j.err
set -euo pipefail

# This template runs code/configuration in the separate EVB_MHETase repository.
: "${CONDA_BASE:?Set CONDA_BASE to your Miniforge/Conda installation directory}"
: "${KEMP_EVB_DIR:?Set KEMP_EVB_DIR to the cloned EVB_MHETase/Kemp_EVB directory}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate kemp-openmm
cd "$KEMP_EVB_DIR"

srun python scripts/full_enzyme.py \
  --mode run \
  --config configs/full_enzyme_short.json \
  --output-dir "runs/full_enzyme_short-${SLURM_JOB_ID}"
