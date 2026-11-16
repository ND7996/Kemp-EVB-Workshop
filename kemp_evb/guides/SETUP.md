# Cluster setup and calculation routes

The repository has two ways to work through the Kemp example:

1. **Run the bundled reduced demo.** This is the recommended first calculation. It includes its mapped input files, scripts, dependency list, and Slurm submission file. It checks an 81-position, five-atom EVB model with OpenMM and EVBOpenMM.
2. **Run the full-enzyme backend.** This uses a separate checkout of `EVB_MHETase/Kemp_EVB`, where the complete system and the preparation, relaxation, heating/equilibration, mapping/FEP, and analysis scripts are maintained.

For cluster commands, Slurm templates, software setup, and the short/long full-enzyme commands, go to [`CLUSTER_QUICKSTART.md`](CLUSTER_QUICKSTART.md). The short calculation is a reduced software demonstration. The full-system workflow uses a changed finite-droplet model and is not a validated reproduction of Q6.

## Dependencies for the bundled demo

The environment file installs Python 3.12, OpenMM 8.6.1, NumPy, and EVBOpenMM at a pinned public revision. Create the environment on a login node or another node with approved package access:

```bash
conda env create -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python kemp_evb/calculation/scripts/run_short_evb_demo.py --check
```

Submit the bundled calculation from the workshop repository root:

```bash
CONDA_BASE="$(conda info --base)" sbatch kemp_evb/calculation/cluster/short_evb_demo.slurm
```

## Full-enzyme workflow stages

The backend’s `full_enzyme.py --mode run` command performs starting-structure preparation, restrained relaxation/minimization, heating, equilibration, and mapping-window production. Its analysis/profile modes connect the sampled state energies and inspect the resulting profile. A shorter preset checks that these stages run; the production preset uses substantially longer sampling. Exact setup and resubmission commands are in the cluster guide and backend documentation.
