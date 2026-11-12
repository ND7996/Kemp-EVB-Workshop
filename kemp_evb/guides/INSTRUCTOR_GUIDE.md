# Instructor guide

## Before the session

Open `index.html` and `presentation.html`. Check that the three links to the saved exercises work. The browser-only workshop needs no installation and does not run a simulation.

## 30-minute session

| Time | Activity |
| --- | --- |
| 0–3 min | Welcome. Tell participants they can follow the saved examples without prior EVB experience. |
| 3–8 min | Use the reaction figure to show what changes from reactant to product. |
| 8–12 min | Explain that EVB connects two bonding descriptions of the same atoms. Define the offset and coupling in plain language. |
| 12–17 min | Open the small EVB viewer. Ask participants to change one control at a time. |
| 17–22 min | Open the saved profile. Identify the axes, barrier and product-versus-reactant energy. |
| 22–28 min | Open the settings comparison. Compare two parameter pairs with similar peak heights. |
| 28–30 min | Recap what EVB shows and what the short example cannot establish. |

## Keep these points clear

- A state is a force-field description of the reactant or product, not a different molecule.
- The energy offset shifts one state relative to the other; coupling controls how strongly they mix.
- Neither parameter is the reaction barrier.
- The plotted horizontal axis is an energy-gap coordinate, not an atom distance or elapsed time.
- The short OpenMM example is for learning the workflow and reading a profile. It does not reproduce the paper's longer Q6 calculation or validate a physical barrier.

## If someone wants to run calculations

The session uses saved browser examples. Optional live calculations need the separate setup in [`SETUP.md`](SETUP.md); they are not part of the timed hands-on.
