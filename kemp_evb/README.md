# Kemp EVB workshop

Use the **[public workshop page](https://nd7996.github.io/Kemp-EVB-Workshop/?v=profile-parameters)** for the introduction, chemistry, actual Q6/OpenMM comparison and full-enzyme run commands.

The main live calculation is now the [full 6,330-atom enzyme workflow](full_enzyme/README.md), using the included top-level EVBOpenMM package. It prepares matching states, relaxes/heats/equilibrates, collects mapping energies and reconstructs a profile. PDB, DCD, CSV/JSON and transactional restart outputs remain local.

The [saved full-system profile](viewers/full_enzyme_recorded_profile.html) is reconstructed from a checksum-verified historical raw record. Blue is the computed OpenMM curve; grey is the published Q6 barrier reference, not a complete paper curve. The changed boundary/parameter transfer and short sampling remain unvalidated.

The browser α/H12 scan reuses saved energies; it starts no new MD.
