# Hands-on guide

Open ../index.html. Begin with the difference between classical MD, QM/MM MD, and EVB. Follow the atom labels through the reaction, then use the saved viewers to change the EVB offset and coupling.

The EVB matrix is a 2 × 2 score table for two bonding descriptions of the same atoms. The diagonal boxes hold state energies; α shifts the product energy; H₁₂ controls mixing. The computer combines these values into an EVB energy at that geometry. A free-energy barrier requires sampling.

For the optional Linux calculation, follow ../calculation/README.md. It is a five-atom software demonstration along a prescribed path.
