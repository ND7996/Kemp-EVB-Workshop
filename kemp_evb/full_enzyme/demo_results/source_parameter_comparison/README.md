# Archived QFEP-parameter comparison

The newly supplied `Kemp_variants.zip` includes `EVB_Simulations/Inputs_EVB/qfep.inp`. It is a QFEP analysis input, despite the archive README referring to it as production. It specifies constant H12=100.4 kcal/mol and native state-2 alpha=−110.0 kcal/mol. Applying the native-to-zero-minimum Morse bookkeeping correction+44.5 gives total P offset−65.5 kcal/mol. No fit to the13.2 enzyme barrier is performed.

The charge/type/bond/angle/torsion/improper parameters and Q-index-normalized changes in the archived GNCA FEP match all four supplied Des27 FEP files. This supports a source-parameter comparison; it does not establish host/boundary or reference-calibration equivalence. Original analysis input and provenance are stored under inputs/protocol_reference. Native binary energy outputs and an independent reference-calibration dataset are not included in the upload.

## Analysis and result

This reconstruction mixes the two energies for each frame first, then performs stable within-window/bin exponential averaging and count-weighted combination of windows, with EVBOpenMM's mapping-FEP normalization. The formula order follows the native Q6 source's QFEP energy mixing and bin averaging. This differs from the earlier illustrative lowest eigenvalue of already-averaged diabatic curves. It is a QFEP-style comparison, not a claim of exact binary QFEP parity.

For the uploaded21-window short run (50 samples/window),18 bins and minimum20 per window/bin resolve a connected curve: conditional barrier **9.572629 kcal/mol** and reaction free energy **−44.452166 kcal/mol**. These are short-run estimates in the changed finite-droplet model, not a reproduced13.2±0.6 result or a newly validated physical calibration. No uncertainty interval is inferred. Finer40-bin coverage remains fragmented; no interpolation is used.

The archived analysis uses51 energy files, kT0.596 kcal/mol, skip4000 records/window,51 bins and minimum1000 per window/bin. The short run cannot support those settings. This comparison uses its actual300K sampling temperature, zero additional skip (simulation equilibration was discarded before saving production records),18 bins and minimum20. Those differences are explicit in summary.json. More saved records alone do not establish independence.

The plot includes the full product well and a labelled detail. Its offset-adjusted energy-gap axis is E_R−E_P−(alpha+Morse correction); no transition-state position or barrier height is forced.

From the workshop repository root, reconstruct the included raw record without new dynamics:

```bash
python kemp_evb/calculation/scripts/run_kemp_evb.py --mode profile --records kemp_evb/full_enzyme/demo_results/source_parameter_comparison/replica_101_raw.json --output-dir kemp_evb/full_enzyme/runs/recorded_comparison --bins 18 --open
```

Default and explicit-parameter profile commands both use frame-wise mixing/reweighting. Explicit values are labelled user trials, not physical calibration. The grey paper barrier in the newly generated public viewer is a scalar mean/SD reference, not a full native curve. The original archived plot files in this folder remain as historical source evidence; the current public figure is generated from the checksum-matched raw record.
