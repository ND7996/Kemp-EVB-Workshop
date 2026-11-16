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

### 2. Relax, heat and equilibrate the starting mapping state

Follow the documented engine protocol. Native Q6 uses restrained low-temperature dynamics, heating and equilibration before sequential mapping windows. Our EVBOpenMM driver minimizes the reactant mapping state, heats and equilibrates it, then advances through the windows. These short examples do not independently equilibrate both endpoints. Check contacts, geometry, energies and temperature. Relaxation does not calculate the barrier; changing the preparation protocol requires documenting and validating that change.

### 3. Run FEP and analyze the profile

Use the checked R/P atom map to define λ windows between the state energy descriptions. Sample each window, inspect neighboring-window overlap, and repeat independent replicas. Analyze the energy records to obtain a free-energy profile, activation barrier, and reaction free energy. Check window coverage and replica agreement before interpreting the result. A single trajectory or a profile that merely matches a target is not enough to establish convergence.

## What this workshop can run

The downloaded materials contain the Kemp reaction mapping, saved short-run data, and browser exercises for EVB concepts and profile analysis. They do **not** include the complete live R/P system builder, endpoint relaxation inputs, or λ-window simulation scripts. The browser examples therefore do not run the three calculation stages above. The launcher checks workshop files; its `core`, `profile`, and `scan` options open the saved Kemp examples.

## Dependencies for later live calculations

The linked scientific repository already supplies runnable scripts, inputs and `Kemp_EVB/environment.yml` for the EVBOpenMM adaptation. It needs Python, NumPy, OpenMM and the documented compatible EVBOpenMM checkout. Native Q6 instead needs the pinned Q6 source build, GNU Fortran, make and its deposited topology/FEP inputs. See the backend README for installation and exact run commands. The standalone browser repository intentionally does not duplicate those engines.

Q6, AmberTools, Rosetta, LAMMPS, SciPy, and Matplotlib are not needed to open or complete this browser hands-on.

## Optional workshop-file rebuilding

Rebuilding the reaction figure, participant handout, or slide PDF requires PyMuPDF. Rebuilding the slide PDF also requires ReportLab. These tools are only needed to edit or regenerate the supplied workshop files; they are not needed to view the site, slides, or PDFs.

## Scientific backend and native Q6 reconstruction

The full scripts are maintained separately in [EVB_MHETase/Kemp_EVB](https://github.com/ND7996/EVB_MHETase/tree/main/Kemp_EVB). Native execution and recovered-template provenance are documented in [NATIVE_Q6_REPRODUCTION.md](https://github.com/ND7996/EVB_MHETase/blob/main/Kemp_EVB/docs/NATIVE_Q6_REPRODUCTION.md). This standalone repository remains a saved-data/browser workshop; its Python launcher opens existing viewers and does not execute those simulations.

## Replica exercise

Fresh shortened Q6 replicas can require about an hour or more each on CPU. Run them before the workshop; use the saved native outputs and replica statistics in the practical. The two-replica runner in the scientific repository handles Qdyn, QFEP and SD/SEM automatically. SD measures between-run spread; SEM = SD/√n estimates mean precision under independence assumptions. With n=2 it is very uncertain and omits common startup/model/binning bias. The published ±0.6 is SD over30replicas, not SEM.
