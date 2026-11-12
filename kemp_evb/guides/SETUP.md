# Installation pre-work

Finish this before the60-minute session. Participants who only use the saved HTML viewers need a browser; they can open index.html directly. Python is needed for the launcher. Numerical recomputation additionally needs the full repository, NumPy, OpenMM and EVBOpenMM.

## Existing checkout

From the EVB_MHETase repository root:

```bash
git pull --no-rebase
conda env create -f Kemp_EVB/environment.yml
conda activate kemp-openmm
python -m pip install -e ../EVBOpenMM
cd workshop/kemp_evb
python run_practical.py --check
```

Skip environment creation if `kemp-openmm` already exists. The pip command assumes EVBOpenMM is a sibling checkout; replace that path with its actual location. The host's established location is `/home/hp/nayanika/github/EVBOpenMM`. No recurring EVBOPENMM_SRC environment variable is required after installing the package.

## Fresh installation

Install Conda/Miniforge and Git first. Choose a working directory, then:

```bash
git clone https://github.com/ND7996/EVB_MHETase.git
git clone https://github.com/Martin-Floor/EVBOpenMM.git
git -C EVBOpenMM checkout dedf301eb9b0021a4be200327cea53f5b7f034c6
cd EVB_MHETase
conda env create -f Kemp_EVB/environment.yml
conda activate kemp-openmm
python -m pip install -e ../EVBOpenMM
cd workshop/kemp_evb
python run_practical.py --check
```

That EVBOpenMM revision matches the inspected/tested core API. Use these clone/checkout commands only for a fresh installation; do not reset an existing checkout containing your work.

## Dependencies

The environment specifies Python3.12, OpenMM8.6.1, NumPy and pip. EVBOpenMM supplies the reaction-force and analysis layer. CPU execution is sufficient for the short teaching-core evaluation and saved-energy analysis. A GPU is optional for later full-enzyme dynamics.

Q6, AmberTools, Rosetta, a Voth/LAMMPS build, SciPy and Matplotlib are **not required for the practical**. The workshop explains Q6; it does not execute it. AmberTools is only relevant to the separate optional Amber-input rebuilding exercise. Browser plots need no plotting-package installation.

If installation is incomplete, open index.html and complete the saved-viewer route in HANDS_ON.md. Do not spend the practical downloading packages or starting long simulations.
