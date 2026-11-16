# Kemp EVB workshop

This folder contains the public Kemp reaction workshop: the website, presentation, saved data activities, and a small runnable OpenMM/EVB calculation.

## Begin with the website

Open [`index.html`](index.html), or visit the [public page](https://nd7996.github.io/Kemp-EVB-Workshop/). The matrix explanation uses plain language first; the equation and interactive calculator are optional.

## Run a live calculation

The reduced five-atom calculation and its reaction inputs are included in [`calculation/`](calculation/README.md). To run it on a Slurm cluster, use [`guides/CLUSTER_QUICKSTART.md`](guides/CLUSTER_QUICKSTART.md). This run checks OpenMM against EVBOpenMM along an imposed set of geometries. It is not a full-enzyme reaction free-energy calculation.

The full-enzyme workflow is maintained in the separate [`EVB_MHETase/Kemp_EVB`](https://github.com/ND7996/EVB_MHETase/tree/main/Kemp_EVB) repository. The cluster guide gives the commands for that short workflow and its longer production configuration.

The saved profile and paper comparison are teaching material; read the labels and limitations alongside each plot. The workshop does not claim that its reduced scan reproduces the paper’s reported barrier.
