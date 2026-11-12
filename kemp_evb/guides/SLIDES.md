# Kemp EVB workshop: slides and speaker notes

## 1. EVB with Kemp elimination: one method, two engines

- 30-minute talk, then30-minute practical.
- Understand states, coupling, offsets and mapping/FEP.
- Compare the paper's Q6 calculation with an OpenMM adaptation.
- Today: saved results and target-matching exercises, not a validated reproduction.

Speaker notes: 2 min. Explain what is executable now and distinguish source comparison from workshop target matching.

## 2. Define the Kemp reaction and mapped states

- Asp162 abstracts the proton from5-nitrobenzisoxazole.
- CH-HT and OR-NR break; Asp OD1-HT forms.
- Other bonded terms, charges and nonbonded parameters also change.
- R/P share the same atoms; the22-atom reactive region has charge −1 in both.

Speaker notes: 3 min. Read numerical identities from the FEP. Some source comments have stale residue numbers. The product is ring-opened nitrile/phenoxide plus protonated Asp.

## 3. What is EVB? What is a diabatic state?

- EVB describes a reaction by mixing classical bonding/charge patterns.
- A diabat is one force-field state, not a second set of atoms.
- Each state defines bonded terms, charges, LJ and exclusion rules.
- Reference calibration and ensemble sampling are needed for free energies.

Speaker notes: 2 min. EVB is a model, not a synonym for either Q6 or OpenMM. Explain why a topology can change while particle identities remain mapped.

## 4. The two-state EVB Hamiltonian

- H = [[V_R, H12], [H12, V_P + alpha]].
- Its lower eigenvalue is the coupled ground-state potential.
- Eigenvector weights describe state mixing.
- The same matrix algebra applies in Q6 and OpenMM.

Speaker notes: 3 min. Write E-minus=(VR+VP+alpha-sqrt((VR-VP-alpha)^2+4H12^2))/2. Distinguish an energy evaluation from a free-energy calculation.

## 5. Coupling, offset and barrier are different quantities

- Alpha sets the relative diabatic energy origin.
- H12 mixes the states and opens an avoided crossing.
- Neither coefficient is the activation barrier.
- Fit independent reference evidence for a physical model; target matching today is a labelled exercise.

Speaker notes: 2 min. A single barrier does not identify two parameters. Do not present a per-enzyme target match as independent validation.

## 6. What are Q6, OpenMM and EVBOpenMM?

- Q6: standalone MD/EVB toolkit, including Qdyn dynamics and QFEP analysis.
- OpenMM: programmable MD engine with Python APIs and CPU/GPU kernels.
- EVBOpenMM: the added EVB force/model/workflow layer used here.
- File formats, analysis, boundaries and precision can differ even when EVB theory is shared.

Speaker notes: 3 min. Neither engine inherently selects one force field. This paper chose OPLS for EVB and Amber for separate conventional MD. Q6 is explained, not installed for the practical.

## 7. Mapping/FEP and profile reconstruction

- Sample lambda*V_R + (1-lambda)*V_P at several mapping windows.
- Save both state energies at every sampled configuration.
- Mix energies per frame, then reweight and bin consistently.
- Lambda is not physical reaction time or the adiabatic state-mixing weight.

Speaker notes: 3 min. Nonlinear operations do not commute: the default source route mixes per-frame energies before exponential bin averaging. A geometric potential scan is not an enzyme PMF. Coverage, correlation and replicas matter.

## 8. What the Nature paper simulated

- Conventional MD: Amber ff19SB/GAFF2 with OPC.
- EVB: Q6, OPLS-AA, TIP3P, SCAAS sphere and local reaction field.
- 30 replicas/case;30 ns equilibration;51 windows x200 ps production.
- Des27.7-IN computed mean13.2 +/-0.6 kcal/mol; uncertainty is replica SD.

Speaker notes: 2 min. These are published settings/results, not the short workshop run. The paper is Listov et al., Nature2025. Experimental apparent barriers derived from kcat are distinct observables.

## 9. What this OpenMM example implements

- All6330 atoms and the supplied reactive changes are imported for Des27.7-IN.
- Workflow: minimization, heating, equilibration and restartable mapping production.
- Boundary changes: direct nonperiodic electrostatics, droplet wall and fixed outer protein.
- No Q6 LRF equivalence or converged physical barrier is claimed.

Speaker notes: 2 min. Existing dynamics/energy checks establish execution, not whole-Hamiltonian equivalence. The practical uses saved short energies rather than starting long MD.

## 10. Recovered parameters and energy conventions

- Archived qfep.inp: constant H12=100.4 and native alpha = -110 kcal/mol.
- Reactive FEP sections match the Des27 inputs after Q-index mapping.
- Q6 Morse minimum=-De; EVBOpenMM minimum=0.
- +44.5 on P preserves the relative convention: total P offset=-65.5 kcal/mol.

Speaker notes: 2 min. The archive is GNCA Kemp-variants data cited by the protocol, not replacement Des27 structures. The44.5 correction is bookkeeping, not a fitted barrier. Original independent reference targets remain unverified.

## 11. Read the source-parameter comparison

- Saved21-window short OpenMM run;50 records/window.
- At18 coarse bins the source-based comparison is about 9.57 kcal/mol.
- No fit to13.2 was used for this curve.
- Changed boundary, short correlated data and coarse binning limit interpretation.

Speaker notes: 2 min. The paper reference is not reproduced. Explain the full product-side reaction free energy, not only the forward peak. Finer bins remain fragmented.

## 12. Many parameter pairs can match one target

- Reuse the saved energies: no additional MD for the workshop.
- 25 grid combinations plus5 examples matching the chosen13.2 barrier.
- One example: H12 about 88.72, native alpha = -110 kcal/mol.
- Matching pairs have different reaction free energies; this is not unique physical calibration.

Speaker notes: 2 min. Keep the archived grey baseline visible. These are explicitly labelled workshop target matches, not predictions.

## 13. Grey baseline, profile maximum and coordinate zero

- Grey dashed: archived baseline on the same energies; blue: selected candidate.
- The diamond marks a binned maximum, not a validated molecular TS.
- Raw or corrected gap zero need not be a free-energy maximum.
- TS-centred view applies the same display shift to both curves; energies/barriers do not change.

Speaker notes: 2 min. Example grid H12=80, alpha = -90 has raw peak-74.36 and baseline-127.22. Common centring gives selected0, baseline-52.86. Then begin the participant worksheet. Talk total 30 min.

## 14. Practical: first 15 minutes

- 0-5: inspect the reaction and Q atom mapping.
- 5-10: run the core or open its saved viewer; vary coupling/offset.
- 10-15: read the archived input and recompute the source profile.
- Use run_practical.py or the offline index; the sample energies are included.

Speaker notes: Practical instruction slide, not extra talk time. Commands and browser-only alternatives are in HANDS_ON.md and PARTICIPANT_HANDOUT.pdf. Installation is pre-work.

## 15. Practical: next 15 minutes

- 15-25: compare parameter pairs matching one barrier; record reaction free energies.
- 25-30: select H12=80, alpha = -90 and toggle the common display shift.
- Identify grey baseline, blue maximum marker and coordinate zero.
- Discuss one shared EVB concept and one Q6/OpenMM implementation difference.

Speaker notes: Practical instruction slide. No long simulation is required. Use saved viewers immediately if computation fails. Instructor answers are in INSTRUCTOR_GUIDE.md.

## 16. Backup: five-atom potential exercise

- Three native Morse terms plus selected C-N harmonic terms.
- 81 imposed geometries; illustrative H12/extra offset.
- OpenMM lower eigenvalues checked against EVBOpenMM.
- Enzyme, solvent and most force-field terms are omitted; this is not a PMF.

Speaker notes: Optional. The core teaches matrix algebra and parameter roles; its controls must not be relabelled as full-enzyme parameters.

## 17. Backup: native Morse energy origins

- Q6: D*(exp(-2*a*dr)-2*exp(-a*dr)), with minimum-De.
- Zero-minimum Morse adds De for every active bond.
- For these terms R gains182.75, P gains138.25 kcal/mol.
- P-R shifts by-44.5, so add+44.5 to P for convention equivalence.

Speaker notes: Optional. This is only an energy-origin conversion. It does not calibrate the physical model.

## 18. Backup: full-enzyme execution evidence

- 6330 atoms; corrected native LJ import.
- Minimization, heating, equilibration and three production windows executed.
- 480 MD steps at 0.25 fs; checkpoint pause/resume tested.
- This saved DCD is a0.12 ps engineering test, not an observed reaction.

Speaker notes: Optional. Load the topology and both saved DCD segments in the documented order. This dataset differs from the21-window1050-energy-record sample used for profile/scan exercises.

## 19. Backup: source-profile analysis order

- Native QFEP input supplies100.4 coupling and-110 state offset.
- Constant coupling permits reanalysis from the saved state energies.
- Evaluate the lower eigenvalue per frame before exponential bin averaging.
- The sample cannot support native skip4000/minimum 1000 settings; actual settings are disclosed.

Speaker notes: Optional. No native binary energy/QFEP parity run is claimed. Keep archived parameter recovery separate from physical transfer validation.

## 20. Backup: how to run longer simulations

- Change equilibration_steps and window_equilibration_steps in the chosen JSON.
- production_steps_per_window and lambda_schedule set the production length.
- Use a new output directory for changed settings and a separate replica seed.
- The long template collects10.2 ns; CPU runtime can greatly exceed one night.

Speaker notes: Optional. This is after the workshop. Sampling alone cannot establish a missing physical reference calibration. See the full project README; do not start this during the practical.

## 21. Backup: Amber does work with OpenMM

- OpenMM reads Amber prmtop/inpcrd files.
- AmberTools is needed to rebuild the optional Amber lane, not this practical.
- That classical compatibility example is separate from OPLS EVB.
- No Q6, Rosetta or Voth/LAMMPS build is required for today.

Speaker notes: Optional. Do not mix the paper's conventional Amber MD lane with its Q6/OPLS reaction lane.

## 22. Paper, supporting information and deposited inputs

- Listov et al., Complete computational design of high-efficiency Kemp elimination enzymes, Nature (2025).
- Article: doi:10.1038/s41586-025-09136-2; associated Supplementary Information, Table 3.
- Simulation inputs: Zenodo doi:10.5281/zenodo.14563437 (EVB_Simulations.zip).
- Design code: Fleishman-Lab/denovoKemp; design code and EVB inputs have different roles.

Speaker notes: The paper reports computationally designed Kemp eliminases. EVB helps interpret catalytic energetics; the workshop focuses on Des27.7-IN. Cite the SI for protocol details and distinguish computed barriers from kcat-derived apparent experimental barriers.

## 23. Related protocol evidence: do not conflate the datasets

- Gutierrez-Rus et al., JACS (2025), doi:10.1021/jacs.4c09428, with its Supporting Information.
- GNCA Kemp inputs: Zenodo doi:10.5281/zenodo.12666978; Kemp_variants.zip.
- Recovered qfep.inp supplies H12=100.4 and alpha=-110 kcal/mol.
- GNCA is a different scaffold: these are protocol evidence, not replacement Des27 coordinates or proof of parameter transfer.

Speaker notes: The related archive resolves analysis parameters, not the validity of the changed OpenMM environment. Use SOURCES.md for attribution and the original authors’ supporting information.

## 24. Atom names: read the actual mapping

- Actual PDB identities: SUB257 and catalytic ASP162; stale FEP comments use different residue numbers.
- CH=proton-bearing substrate carbon; HT=transferred H; NR=ring nitrogen; OR=ring oxygen.
- Asp OD1 accepts HT; OD2 is the other oxygen; CG is its carboxyl carbon.
- Q indices are local reactive indices; PDB serials are global atom identifiers. Atom names are not force-field types.

Speaker notes: Use the full 22-atom table on the HTML page. Q1=CH serial4002; Q4=OR4005; Q5=NR4006; Q13=HT4014; Q21=OD12526. The PDB shows Asp162 and substrate257, not the stale comments.

## 25. Detailed Kemp mechanism and state bookkeeping

- R: CH-HT and OR-NR bonded; substrate neutral, Asp carboxylate -1.
- Coupled proton abstraction and ring opening: OD1 accepts HT, CH-NR gains nitrile character, OR gains phenoxide character.
- P: OD1-HT bonded; CH-HT and OR-NR absent; nitrile/phenoxide product -1 and neutral AspH.
- Total reactive charge stays -1. This two-state model introduces no separately sampled carbanion intermediate.

Speaker notes: These chemical events are coupled; the bullets are explanatory bookkeeping, not proof of a timed sequence. EVB mixes endpoint force fields rather than explicitly propagating electrons. Inspect geometric distances separately from the energy-gap coordinate.

## 26. Published Q6 versus our short OpenMM example

- SI Table3: Des27.7-IN Q6 barrier=13.2 +/-0.6 kcal/mol (replica SD).
- Short archived-parameter OpenMM estimate=9.57; no independent-replica uncertainty.
- Difference=-3.63 kcal/mol; different boundary and limited sampling prevent an agreement claim.
- A target-matched 13.2 workshop curve is not independent reproduction.

Speaker notes: Original article intro/mechanism/results pages and SI Table3 are available in the HTML with full attribution. Native Q6 was not rerun here. Experiment Des27.7=16.7+/-0.3 is an apparent barrier, not the same observable.
