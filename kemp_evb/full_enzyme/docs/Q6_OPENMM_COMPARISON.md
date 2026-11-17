# What the restored workflow preserves and changes

The common method is two-state EVB: two energy descriptions of identical nuclei, an off-diagonal coupling H12, a relative state offset alpha and free-energy sampling/reweighting. OpenMM is an engine; EVBOpenMM constructs the EVB force model. Qdyn is Q6's dynamics engine and QFEP its native free-energy postprocessor. Differences below describe this particular adaptation, not universal limits of either engine.

| Component | Native Q6 route / source evidence | Restored OpenMM adaptation |
|---|---|---|
| State setup | Deposited Q topology, FEP reaction definition and stage `.inp` files | Parses Q coefficient tables, masses, charges, bonded state changes and exceptions; constructs two prebuilt states and three reactive Morse pairs |
| Starting system | Deposited Des27.7-IN prepared coordinates and droplet topology | Same default 6,330-particle inventory, including 771 waters and 22 mapped reactive atoms |
| Boundary | Native spherical environment and Q boundary/surface treatment | Nonperiodic droplet, water-oxygen harmonic wall, native excluded outer atoms assigned zero mass |
| Electrostatics | `lrf on`; shown production template ordinary cutoffs 10 Å, Q-atom cutoff 99 Å | Direct interactions with 9.9 nm cutoff covering the finite droplet; no native LRF/SCAAS equivalence |
| Relaxation and thermostat | Native input-stage sequence, including restrained dynamics; exact author preparation is not reconstructed by the OpenMM optimizer | Displacement-limited constrained descent, L-BFGS, staged heating and LangevinMiddle dynamics; native production sequence restraints not copied exactly |
| Constraints and timestep | Shown production control enables solvent SHAKE and 1 fs steps | Three rigid-water distance constraints per water; 0.25 fs default for flexible reactive bonds |
| Production length | Shown input: 200,000 × 1 fs = 200 ps/window; 51-window reconstruction gives 10.2 ns/replica | Short preset: 21 × 1,000 × 0.25 fs = 5.25 ps total production; a longer unvalidated template is included |
| Saved samples | Native Qdyn `.en`, trajectories and restart files | JSON/CSV diabatic energies, DCD trajectories, XML state and transactional SQLite samples/checkpoint cursors |
| Analysis | Native QFEP, with case-specific control settings requiring verification | QFEP-style frame-wise EVB mixing and conditional-bin exponential reweighting; not native binary QFEP parity |

The included `inputs/native_protocol_reference/production_example.inp` is an unmodified nonempty native control template from the historical deposited-input reconstruction. It is evidence for the displayed settings, not a complete native Q6 run deck. Some deposited files were empty; an exact Q6 reproduction also needs documented handling of those missing stages and appropriate original analysis/settings.

The profile default uses H12=100.4 and alpha=−110 kcal/mol from the **separate GNCA4-WT QFEP archive** (Zenodo 12666978), not verified Des27.7-IN author calibration. A +44.5 kcal/mol product-state energy-zero conversion accounts for the different Morse convention. Semantic agreement of reactive FEP terms does not establish transferability of coupling/offset across enzymes or boundary treatments.

The saved 21-window OpenMM record reconstructs a conditional 9.572629 kcal/mol barrier at 18 coarse bins, with conditional reaction free energy about −44.45 kcal/mol. The published Des27.7-IN Q6 computed reference is 13.2 ± 0.6 kcal/mol (replica mean ± SD). The reaction energy, model transfer and sampling are not validated. These two barrier numbers must not be described as equivalent or explained solely by run length.

In particular, more sampling reduces sampling error within the selected Hamiltonian; it cannot by itself remove a changed electrostatic or restraint model. Exact reproduction would require numerical energy/force checks against Q6, consistent boundary forces and restraints, reference calibration, appropriate window lengths/overlap and independent replicas.

Sources: Listov et al., Nature (2025), DOI 10.1038/s41586-025-09136-2; deposited inputs 10.5281/zenodo.14563437; separate QFEP reference 10.5281/zenodo.12666978; the Q6 manual; and the included native-control example/source provenance. Code-level checks validate specific conversions and execution paths, not complete native Hamiltonian parity.
