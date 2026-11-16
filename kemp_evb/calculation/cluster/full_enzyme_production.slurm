#!/bin/bash
#SBATCH --job-name=kemp-full-production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=kemp-full-production-%j.out
#SBATCH --error=kemp-full-production-%j.err
set -euo pipefail

# Submit a fresh run with RESUME=0. Resubmit with RESUME=1 to continue it.
: "${CONDA_BASE:?Set CONDA_BASE to your Miniforge/Conda installation directory}"
: "${KEMP_EVB_DIR:?Set KEMP_EVB_DIR to the cloned EVB_MHETase/Kemp_EVB directory}"
: "${KEMP_RUN_DIR:?Set KEMP_RUN_DIR to a persistent output directory}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate kemp-openmm
cd "$KEMP_EVB_DIR"

resume_args=()
if [[ "${RESUME:-0}" == "1" ]]; then resume_args+=(--resume); fi
srun python scripts/full_enzyme.py \
  --mode run \
  --config configs/full_enzyme_production.json \
  --output-dir "$KEMP_RUN_DIR" \
  --max-runtime-seconds 39600 \
  "${resume_args[@]}"
