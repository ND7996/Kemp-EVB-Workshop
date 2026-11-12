# Slide notes for the hands-on

Keep the workshop page open so participants can open each saved Kemp viewer. No package installation or simulation run is needed.

## 1. Welcome

Say: “We will use one reaction to see what EVB does. You do not need to know the equations or the software beforehand.” Point out that the examples use saved data.

## 2. The reaction

Asp162 accepts a proton from the substrate. The C–H bond breaks as the ring opens. Use the figure to distinguish the reactant and product descriptions. Avoid starting with atom serial numbers.

## 3. The two states

Explain that the reactant and product states describe the same atoms with different bonding and charge patterns. They are not two separate molecules moving through the simulation.

## 4. EVB settings

The energy offset α shifts the product description relative to the reactant. H₁₂ controls how strongly the two descriptions mix. Neither setting is the reaction barrier. Keep the explanation at this level unless the group asks for the equation.

## 5. Q6 and OpenMM

Paper Figure 3e shows activation-barrier estimates, not a free-energy profile. The workshop's OpenMM viewer shows a short profile from 21 windows and 50 records per window. The paper uses Q6 with OPLS-AA, TIP3P, SCAAS solvent, and a local reaction field. Boundary treatment, electrostatics, Morse energy reference, and sampling differ, so these are not direct profile-to-profile results. Link to the complete paper figure and OpenMM profile viewer from the slide.

## 6. The full calculation sequence

Explain that a live Kemp calculation has three jobs: prepare matching R and P states, relax each endpoint, then sample λ windows and analyze the FEP profile. The workshop describes these stages but does not include the live run scripts.

## 7. Try the core example

Ask participants to change one control at a time. Have them say what changed in the curve. The small example introduces the controls; it is not the full enzyme free-energy calculation.

## 8. Read a profile

The horizontal energy-gap axis is not a distance or time. The vertical axis is free energy relative to the reactant. The peak estimates a barrier; the product level gives the reaction free energy. A binned maximum is not itself a transition-state structure.

## 9. Compare settings

Invite participants to find two parameter pairs with similar peak heights, then compare their product energies. The point is that one barrier does not determine all model parameters.

## 10. Close

Repeat the take-home: EVB connects reactant and product bonding descriptions. Real predictions also need suitable reference data and adequate sampling. The short OpenMM example is a teaching result, not a reproduction of the paper’s Q6 calculation.
