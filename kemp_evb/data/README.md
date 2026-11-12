# Data used in the practical

`short_run/replica_101_raw.json` is a supplied short, completed 21-window OpenMM adaptation run, with 50 energy records per window. The workshop uses it as saved Kemp data for teaching profile analysis; it does not rerun the simulation. Its original metadata and bytes are preserved, and its SHA256 is recorded in PROVENANCE.json. It is not a native Q6 energy file or a published-paper trajectory. These short correlated samples do not establish convergence.

`des27_reaction.fep` is the original Des27.7-IN reaction mapping/parameter input from the supplied Nature-paper archive. Some comments use stale residue numbering; numerical identities are authoritative.

`qfep_reference.inp` is the original archived analysis input from the supplied GNCA Kemp-variants dataset. The archive README calls it a production input, but its format is QFEP analysis: H12=100.4, native alpha−110, constant coupling. It is used as comparison evidence, not proof of transfer to the changed OpenMM boundary.
