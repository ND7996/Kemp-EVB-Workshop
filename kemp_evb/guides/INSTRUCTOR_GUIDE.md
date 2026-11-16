# Instructor guide

Start with `index.html`, then use `presentation.html`. The saved activities open without installing software. For a live cluster run, point participants to [`CLUSTER_QUICKSTART.md`](CLUSTER_QUICKSTART.md).

## Explain the matrix simply

A matrix is a small table. EVB compares reactant and product bond descriptions for the same atoms. Their energies go on the diagonal; α moves the product energy; H₁₂ controls mixing in the two other cells. The calculation combines those values and uses the lower mixed energy. The matrix is a way to calculate energy, not a picture of two molecules and not a reaction barrier.

## Explain the calculations

The included live run is a reduced five-atom Kemp model. It uses the published mapped inputs, selected reactive terms, and an imposed 81-position geometry path. OpenMM evaluates the energies; EVBOpenMM is used to cross-check the lower two-state energy. Protein, solvent, most force-field terms, dynamics, path optimization and free-energy sampling are omitted. The demonstration settings are illustrative, not calibrated.

The full-enzyme files are kept in the separate `EVB_MHETase/Kemp_EVB` repository. Its workflow includes preparation, relaxation, heating/equilibration, mapping-window production and profile analysis. The current model uses a changed finite-droplet boundary, so do not describe it as an exact reproduction of the Q6 paper calculation. Longer sampling alone does not establish convergence or calibration.

## If participants ask about the paper plots

The paper's Figure 3e reports Q6 activation-barrier estimates, not a full free-energy curve. The saved OpenMM profiles use a different boundary and sampling. Treat them as separate examples of what each calculation reports, not as a direct engine comparison.
