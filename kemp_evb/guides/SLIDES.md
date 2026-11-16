# Slide notes for the Kemp hands-on

Open `index.html` and `presentation.html`. Saved examples run in a browser. The separate cluster guide explains how to submit the bundled live calculation and the longer full-enzyme route.

## 1. The reaction

Asp162 accepts a proton from the substrate. The C–H bond breaks as the ring opens. C3 is the substrate carbon that loses Ht; Ht is the transferred hydrogen; O1 and N2 are the ring oxygen and nitrogen; Asp162 OD1 is the accepting enzyme oxygen.

## 2. Reactant and product

The two states are different bond and charge descriptions of the same atoms. They are not two independent molecules or separate simulation paths.

## 3. The EVB table

Call the matrix a small 2 × 2 score table. The diagonal cells contain reactant and product energies. α shifts the product energy; H₁₂ in the other two cells controls how strongly the descriptions mix. The program combines the values and takes the lower mixed energy for that geometry. That energy is not itself a free-energy barrier. Only show the equation if the group wants more detail.

## 4. Q6 and OpenMM

Paper Figure 3e reports Q6 barrier estimates; the saved OpenMM profile is a different calculation with different model boundaries, electrostatics, Morse energy reference and sampling. Do not compare them as a direct engine benchmark.

## 5. Try the examples

Use the saved browser activities to explore the mechanism, profiles and parameter sensitivity. A change that makes a curve resemble a target does not validate the parameters.

## 6. Run on the cluster

The workshop repository bundles a live five-atom calculation, its mapped inputs and Slurm file. Point out that it checks software and EVB energy mixing along an imposed geometry path; it is not MD or a physical enzyme barrier. The cluster guide also gives commands for the separate full-enzyme short workflow and longer production configuration.

The full-enzyme path has stages for preparing matched states, relaxing, heating/equilibration, mapping-window production and profile analysis. Its current finite-droplet model is not an exact Q6 reproduction.
