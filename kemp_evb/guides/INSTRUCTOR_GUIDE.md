# Instructor guide — Kemp EVB,60 minutes

The workshop is ready as a teaching/demo pack. Independent physical calibration/convergence of the changed full-enzyme model remains incomplete. Do not present a target match as the paper's reproduced barrier.

## Before the audience arrives

1. Pull the complete repository and run `python start_workshop.py --check` from the workshop/kemp_evb folder.
2. Open index.html. Check all three viewers and the slides/handout PDFs. Saved browser exercises need no numerical environment or network.
3. For live numerical evaluation, complete SETUP.md and run `python run_practical.py --check`. CPU is sufficient. Run the three practical steps once with `--no-open` if testing on the presentation laptop.
4. Use slides1–13 for the30-minute talk; slides14–15 introduce the practical. Slides16 onward are optional backup. Speaker notes are in SLIDES.md.
5. Keep the browser-only route ready. Do not build environments or run full-enzyme production during the timed session. The small raw-energy sample is bundled; no participant needs an existing local simulation.

## Talk pacing (exactly30 minutes)

| Slides | Minutes | Message |
| --- | ---: | --- |
| 1–3 |7 | Outcomes, mechanism, diabatic states. |
| 4–5 |5 | Hamiltonian, coupling versus offset. |
| 6–7 |6 | Q/OpenMM roles; mapping and statistical analysis. |
| 8–9 |4 | Published setup versus our changed model. |
| 10–11 |4 | Archived values, energy-zero conversion and source comparison. |
| 12–13 |4 | Parameter non-uniqueness, grey baseline and coordinate zero. |

## Practical answer key

**0–5, reaction:** break CH–HT and OR–NR; form Asp162 OD1–HT. More bonded/nonbonded terms change than these three net bond changes. The same atoms are mapped and the22-atom reactive region has total charge−1 in both states. A diabatic state is a force-field description of a bonding/charge pattern, not another set of particle identities.

**5–10, core:** H12 mixes states/opens an avoided crossing; alpha shifts the relative state origin. Neither is the activation barrier. The five-atom imposed-path potential omits the enzyme, solvent, most force terms and ensemble entropy. It is not the enzyme PMF. The live step actually evaluates OpenMM; offline sliders reuse evaluated diabats.

**10–15, source:** original QFEP input gives H12=100.4 and alpha −110 kcal/mol. The+44.5 correction restores the relative energy convention for the supplied zero-minimum Morse wells, so total P offset is−65.5. Do not add it twice or treat44.5 as a barrier. The source-based coarse short-run barrier is about 9.57; paper Des27.7-IN computed mean is13.2±0.6 (replica SD). Boundary/sampling/reference differences and short correlated data prevent a reproduction claim.

**15–25, combinations:** examples matching the imposed13.2 forward target include approximately (H12,alpha)=(75.90,−130), (82.45,−120), (88.72,−110), (97.96,−100), (107.61,−90), all in kcal/mol. Their conditional reaction free energies differ, roughly−69.0,−58.3,−47.8,−36.6,−25.6. Matching one number cannot identify two parameters or validate the entire curve. The grey baseline uses the same saved data with archived parameters; it is not an experimental trace.

**25–30, maximum/zero:** select grid H12=80, alpha −90. Its estimated raw peak is−74.3595; the archived baseline estimate is−127.2181. A common+74.3595 horizontal shift gives selected0 and archived−52.8586. Original data/energies/barriers remain unchanged. Raw or corrected gap zero describes a diagonal-energy equality, not a requirement for a free-energy maximum.18 bins yield a coarse maximum estimate, not a validated molecular TS or a committor result.

**Q/OpenMM:** both can implement EVB states, coupling and mapping; engines/file formats/APIs/hardware differ. This paper used Q6/OPLS-AA/TIP3P with SCAAS/LRF. OpenMM is a programmable MD engine; EVBOpenMM is the added model/workflow layer. The current droplet changes the boundary/electrostatics. OpenMM supports Amber input formats, but that does not turn this paper's OPLS EVB model into an Amber model.

## Optional full-enzyme execution and trajectory

After the workshop, run the project workflow in ../../Kemp_EVB using its README. The saved full-enzyme engineering trajectory is in `../../Kemp_EVB/demo_results/full_enzyme_engineering/`. In VMD load topology.pdb, then trajectory.dcd, then trajectory_resume.dcd. It is a480-step,0.12 ps mapping test, not an observed reaction or equilibrated production. The21-window source/scan dataset is a separate completed short run,1050 energy records; do not confuse the two datasets.

For longer runs change `equilibration_steps`, `window_equilibration_steps`, `production_steps_per_window` and `lambda_schedule` in the chosen JSON before a new run. Sample/trajectory/checkpoint intervals control reporting. Use a new output directory for changed settings. The51-window long template collects10.2 ns production and can take much longer than a night on CPU. This is outside the practical; length alone does not validate the model.

## Troubleshooting and fallback

- Missing numerical packages: activate kemp-openmm and follow SETUP.md. Use index.html meanwhile.
- No browser opens: open index.html or the printed output HTML path manually.
- No local run folder: the practical runner uses the included sample; do not ask participants to create long MD first.
- Repeated core run: the runner chooses a new numbered folder. A user-specified existing core output is preserved and rejected.
- Different raw data in an output directory: select a new `--output-dir`; the runner will not overwrite it.
- Fragmented finer profile: retain gaps. Changing to18 coarse bins is a reporting choice, not new sampling or convergence.
- Changed numerical/code version: the energy data remain analysable; do not promise bitwise trajectory restart across changed platforms/code.

## Rebuild materials (instructor/developer only)

The delivered slides/viewers need no rebuild. To regenerate after editing this folder, install matplotlib, python-pptx and reportlab in a materials environment, then from the repository root:

```bash
python Kemp_EVB/workshop_tools/build_materials.py
python workshop/kemp_evb/start_workshop.py --check
```

This rebuilds slides, participant PDF, copied viewers, landing page and hash manifest. It does not run MD. The original small raw-energy sample must remain present and match its provenance. Readiness checks cover file integrity and live practical execution; they do not establish physical calibration/convergence.

For MHETase/PTE extensions, reuse the state-mapping/sampling/validation ideas, not the Kemp numerical parameters. A new reaction requires its own mechanism, reference evidence and reactive parameters.
