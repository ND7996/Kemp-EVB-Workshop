# Kemp workshop setup and calculation stages

Everything in this repository is for the **Kemp reaction**. The hands-on uses saved examples and needs only a current web browser. Open `index.html`; no installation, account, or simulation run is required.

To check that the downloaded workshop files are intact, install Python 3 and run this from the `kemp_evb` folder:

```bash
python start_workshop.py --check
```

The interactive pages are self-contained HTML. The optional launcher also uses only Python's standard library.

## Full Kemp EVB calculation

The full calculation has three stages. The stage order matters because each stage depends on checked output from the previous one.

### 1. Prepare the reactant and product states

Begin with the Kemp enzyme and substrate structure. Build a reactant (R) and a product (P) state with the same atoms in the same order. Assign the state-specific bond and charge descriptions, force-field terms, and EVB parameters. Check the atom map, total charge, proton placement, residue identities, solvent and ions. Save the structure and parameter sources so both endpoints can be reproduced. Do not proceed until R and P load correctly and their mapping is verified.

### 2. Relax the endpoints

Minimize and equilibrate R and P separately, using the same environment and a documented restraint schedule. Check for clashes, distorted geometry, unstable energies, and sensible temperature and density. Relaxation prepares the endpoint structures; it does not calculate the reaction barrier.

### 3. Run FEP and analyze the profile

Use the checked R/P atom map to define λ windows between the state energy descriptions. Sample each window, inspect neighboring-window overlap, and repeat independent replicas. Analyze the energy records to obtain a free-energy profile, activation barrier, and reaction free energy. Check window coverage and replica agreement before interpreting the result. A single trajectory or a profile that merely matches a target is not enough to establish convergence.

## What this workshop can run

The downloaded materials contain the Kemp reaction mapping, saved short-run data, and browser exercises for EVB concepts and profile analysis. They do **not** include the complete live R/P system builder, endpoint relaxation inputs, or λ-window simulation scripts. The browser examples therefore do not run the three calculation stages above. The launcher checks workshop files; its `core`, `profile`, and `scan` options open the saved Kemp examples.

## Dependencies for later live calculations

When the Kemp run scripts and inputs are added, the calculation environment will need Python, NumPy, OpenMM, and a compatible EVBOpenMM checkout. A CPU is sufficient for small checks; longer enzyme sampling may use a supported GPU. The current repository does not provide a pinned Conda environment or working live-run commands, so do not treat the saved-data exercises as a substitute for those inputs.

Q6, AmberTools, Rosetta, LAMMPS, SciPy, and Matplotlib are not needed to open or complete this browser hands-on.

## Optional workshop-file rebuilding

Rebuilding the reaction figure and participant handout requires PyMuPDF. Rebuilding the slide PDF from its source requires ReportLab. These tools are only needed to edit or regenerate the supplied workshop files; they are not needed to view the site, slides, or PDFs.
