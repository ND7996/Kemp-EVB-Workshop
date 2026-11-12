# Instructor guide

## Before the session

Open `index.html` and `presentation.html`. Check that the three links to the saved exercises work. The browser-only workshop needs no installation and does not run a simulation.

## Suggested session flow

1. Welcome participants and explain that the saved examples need no prior EVB experience.
2. Use the reaction figure to show the proton transfer and ring opening.
3. Explain how EVB connects two bonding descriptions of the same atoms. Define the offset and coupling in plain language.
4. Open the EVB viewer and invite participants to change one control at a time.
5. Read the saved profile together. Identify the axes, barrier estimate, and product-versus-reactant energy.
6. Compare parameter pairs with similar peak heights.
7. Recap what the examples show and what they cannot establish.

## Keep these points clear

- A state is a force-field description of the reactant or product, not a different molecule.
- The energy offset shifts one state relative to the other; coupling controls how strongly they mix.
- Neither parameter is the reaction barrier.
- The plotted horizontal axis is an energy-gap coordinate, not an atom distance or elapsed time.
- The short OpenMM example is for learning the workflow and reading a profile. It does not reproduce the paper's longer Q6 calculation or validate a physical barrier.

## If someone wants to run calculations

The session uses saved browser examples. The added workflow slide explains the three stages of a full Kemp calculation: prepare matching R/P states, relax both endpoints, then sample and analyze the FEP profile. Those stages are not run during the hands-on, and the live R/P and FEP scripts are not included in this workshop package. See [`SETUP.md`](SETUP.md) for the current scope and dependencies.
