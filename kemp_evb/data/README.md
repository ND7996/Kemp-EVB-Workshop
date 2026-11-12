# Data used in the practical

`short_run/replica_101_raw.json` is the small user-supplied completed21-window OpenMM adaptation run, with50 energy records/window. It is included so a fresh lab checkout can recompute the source comparison and workshop scan without running long MD. Its original metadata and bytes are preserved; the SHA256 is checked against the recorded comparison and listed in PROVENANCE.json. It is not a native Q6 energy file or a published-paper trajectory. Short correlated samples do not establish convergence.

`des27_reaction.fep` is the original Des27.7-IN reaction mapping/parameter input from the supplied Nature-paper archive. Some comments use stale residue numbering; numerical identities are authoritative.

`qfep_reference.inp` is the original archived analysis input from the supplied GNCA Kemp-variants dataset. The archive README calls it a production input, but its format is QFEP analysis: H12=100.4, native alpha−110, constant coupling. It is used as comparison evidence, not proof of transfer to the changed OpenMM boundary.
