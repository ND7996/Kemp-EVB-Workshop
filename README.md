# Kemp EVB workshop

A beginner-friendly hands-on about the Kemp reaction and empirical valence bond (EVB). Open the [public workshop](https://nd7996.github.io/Kemp-EVB-Workshop/) or open `kemp_evb/index.html` from a clone.

## Start here

1. Read the short reaction and EVB explanations on the website.
2. Use the saved browser exercises to explore the mechanism and energy profiles.
3. To run a live calculation, follow [`kemp_evb/calculation/README.md`](kemp_evb/calculation/README.md). It includes the complete reduced-demo inputs and Slurm job file.
4. For cluster setup and the separate full-enzyme workflow, follow [`kemp_evb/guides/CLUSTER_QUICKSTART.md`](kemp_evb/guides/CLUSTER_QUICKSTART.md).

The 2 × 2 EVB matrix is a score table for two bonding descriptions of the same atoms. Its diagonal entries are the reactant and product energies; α shifts the product energy; H₁₂ controls how strongly the descriptions mix. The lower mixed energy is the EVB energy at that geometry. The website keeps the equation and calculator in an optional section.

The short live demo runs an 81-position, reduced five-atom path with OpenMM and checks energies against a pinned EVBOpenMM library. It is designed to show the software working; it is not molecular dynamics, an optimized pathway, a free-energy barrier, or a reproduction of the full enzyme calculation. The separate backend repository contains the full-enzyme workflow.
