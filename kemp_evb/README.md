# Kemp reaction workshop

A short, beginner-friendly hands-on for the **Kemp reaction**. Open the [public workshop website](https://nd7996.github.io/Kemp-EVB-Workshop/) or open [`index.html`](index.html) in a browser. The saved examples work offline and need no software installation.

The website explains the full calculation sequence—**prepare the reactant and product states, relax/heat/equilibrate the starting mapping state, then sample and analyze the FEP profile**—and shows which of those steps this teaching package can run. The interactive activities use saved Kemp examples.

Open [`presentation.html`](presentation.html) for the slides, [`guides/HANDS_ON.md`](guides/HANDS_ON.md) for the participant worksheet, and [`guides/SETUP.md`](guides/SETUP.md) for the available files, dependencies, and simulation-stage overview.

The matrix calculator, native-Q6 reconstruction status and cited parameter history are now included in the website. The original 10 slides remain, with six optional detail slides. The actual OpenMM driver relaxes the reactant mapping state before window sampling; general advice about independently relaxing endpoints is distinguished from what that driver executed.

<!-- BEGIN NATIVE REPLICA LINKS -->
Actual shortened native Q6 replicas: **9.82, 8.10 kcal/mol**; mean **8.96**, sample SD **1.216**, SEM **0.860** (n=2). [Open the replica exercise](viewers/native_q6_replicas.html). The paper reports13.2±0.6 (SD over30replicas); our two-replica SEM does not establish convergence or accuracy.
