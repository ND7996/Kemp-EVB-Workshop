# Kemp EVB workshop

Open [`index.html`](index.html) to start. The page is arranged in this order:

1. What EVB describes and what a spherical boundary means.
2. The Kemp reaction scheme and the atoms that change.
3. The paper’s Q6 result and the short OpenMM example, including the input files they use.
4. The α and H₁₂ parameter examples.

The Q6 source inputs are a coordinate file (`.pdb.gz`), topology (`.top.gz`) and two-state reaction definition (`.fep`). The short OpenMM demonstration checks these inputs and uses five mapped atoms with selected reactive terms; it is not a full enzyme simulation. The parameter scan uses saved OpenMM energies and does not start a new simulation.
