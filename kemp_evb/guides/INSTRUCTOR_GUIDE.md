# Instructor notes

Start with the methods: classical MD moves atoms but usually keeps bonds fixed; QM/MM treats reacting atoms quantum mechanically and the surroundings classically; EVB couples reactant and product force-field states. The paper uses Q6 for EVB. The short workshop example uses OpenMM with EVBOpenMM.

Then follow the Kemp reaction and explain the 2 × 2 table once: the diagonal entries are state energies, α shifts the product energy, and H₁₂ controls mixing. Neither setting is the barrier.

The paper reports Q6 activation-barrier estimates. The workshop profile is a short OpenMM example using a different setup and sampling, so the results are not a direct software comparison.

The optional five-atom Linux example checks software behavior along a prescribed path. It is not a molecular-dynamics or free-energy calculation.
