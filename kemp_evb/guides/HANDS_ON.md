# Hands-on: the Kemp reaction and EVB

Start at [`../index.html`](../index.html). Follow the atom labels through the reaction, then use the saved examples to see how an EVB profile changes when its settings change.

## The matrix in plain language

The EVB matrix is a small 2 × 2 table. One diagonal cell is the energy of the reactant bonding description; the other is the product energy plus the offset α. Both off-diagonal cells contain H₁₂, which controls how much the two descriptions mix. These are two descriptions of the same atoms, not two independent molecules. The lower mixed energy is the EVB energy at a given geometry. The matrix does not itself equal a free-energy barrier.

## Run the reduced live calculation

Follow [`../calculation/README.md`](../calculation/README.md) for dependencies and local execution, or [`CLUSTER_QUICKSTART.md`](CLUSTER_QUICKSTART.md) for login-node setup and Slurm submission. The bundled calculation uses the mapped `Des27_7_IN` source input, evaluates 81 positions along an imposed path, and checks OpenMM against EVBOpenMM.

It includes only five atoms and selected reactive terms. Protein/solvent effects, most force-field terms, molecular dynamics, path optimization, and free-energy sampling are omitted. Its H₁₂ and energy offset are illustrative values, so its output is not a physical activation barrier or a reproduction of the paper.

## Run the full system

The separate [EVB_MHETase/Kemp_EVB backend](https://github.com/ND7996/EVB_MHETase/tree/main/Kemp_EVB) contains the mapped full-enzyme workflow. It has stages for preparation, relaxation, heating/equilibration, mapping-window production, and profile analysis. The cluster guide includes submission templates for the short end-to-end run and longer sampling preset. Longer sampling alone does not prove convergence or parameter calibration.

## Source data

The mapped input files bundled for the reduced demo are from the published Kemp EVB dataset, [Zenodo DOI 10.5281/zenodo.14563437](https://doi.org/10.5281/zenodo.14563437), associated with Listov et al., *Nature* (2025), [10.1038/s41586-025-09136-2](https://doi.org/10.1038/s41586-025-09136-2).
