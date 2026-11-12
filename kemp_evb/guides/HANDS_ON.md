# Participant worksheet — 30 minutes

Open `index.html` in a browser and follow the links below. The saved examples work offline; this 30-minute activity does not need a terminal or package installation.

## 0–5 min: define the reaction and states

Open `data/des27_reaction.fep` and the state diagram on index.html.

Asp162 abstracts the substrate proton. The CH–HT and OR–NR bonds break, Asp OD1–HT forms, and other bonded/nonbonded terms also change. R/P map the same atoms; the22-atom reactive region has total charge−1 in both states. Q atom13 maps to global serial4014 (HT); Q atom21 maps to2526 (Asp162 OD1). Some source comments use stale residue numbers.

**Write down:** Which bond forms? What is a diabatic state? Why must total charge be consistent between states?

## 5–10 min: explore the EVB teaching core

Open **EVB fundamentals** on index.html. Move the coupling and extra-offset sliders. Observe the avoided crossing and state mixing. The sliders update a saved teaching example; they do not launch molecular dynamics.

**Write down:** How do coupling and offset differ? Is this geometric five-atom potential curve an enzyme activation free energy?

## 10–15 min: source parameters and the saved Kemp profile

Read `data/qfep_reference.inp`, then:

Open **Archived-parameter comparison** on index.html.

The archived input gives H12=100.4 and native alpha −110 kcal/mol. Q6 uses negative-minimum Morse wells; the zero-minimum EVBOpenMM convention requires+44.5 on P, giving total P offset−65.5. This is bookkeeping, not the activation barrier. The short dataset yields about 9.57 kcal/mol at 18 coarse bins. The paper's Des27.7-IN mean is13.2±0.6 kcal/mol (replica SD), from different boundary/sampling conditions.

**Write down:** Why is H12 not the barrier? Why is9.57 not established reproduction of13.2? Which engine produced this sample dataset?

## 15–25 min: compare parameter combinations

Open **Workshop combinations** on index.html.

Use the dropdown/table to compare25 grid pairs and five refined examples matching the chosen13.2 target. Find two different H12/alpha pairs with the same forward barrier. Record their reaction free energies. One selected example holds alpha −110 and uses H12 ≈88.72; it is a target match, not an independent prediction.

**Legend:** grey dashed = archived baseline on the same saved OpenMM energies; blue = selected candidate; red dashed = chosen target. The grey curve is not an experimental trace. Leave it visible while comparing.

**Write down:** Two matching pairs and their reaction ΔG. Does one barrier determine two parameters uniquely? Does a target match establish convergence or physical calibration?

## 25–30 min: explain the maximum and coordinate zero

Select grid H12=80, alpha −90. Find the blue diamond (binned maximum) and the open grey baseline marker. The displayed raw peak is about−74.36 kcal/mol, with a coarse bin width about 52.86; the maximum is not required to be at raw gap zero.

Toggle **Centre selected profile maximum at x = 0 (display only)**. Both curves receive the same horizontal shift. The selected marker goes to0; the grey marker goes to about−52.86. The original coordinates, energies and barrier heights do not change.

**Write down:** What does raw x=0 mean? What does centring prove? Is this binned maximum a validated molecular transition state? Name one shared EVB concept and one Q6/OpenMM implementation difference.

## What to submit or discuss

A screenshot of two matching parameter examples, their reaction free energies, and one sentence explaining why they are workshop target matches rather than validated predictions. Add a sentence identifying the grey curve and explaining the display-only centring operation.

This worksheet analyzes the bundled saved Kemp data. It does not prepare R/P systems, relax endpoints, or run λ-window FEP. Those three stages are described in [`SETUP.md`](SETUP.md) and the website's calculation workflow section.
