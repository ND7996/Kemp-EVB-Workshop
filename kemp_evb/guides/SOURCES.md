# Sources, attribution and result labels

1. **Listov et al.** Complete computational design of high-efficiency Kemp elimination enzymes. Nature (2025). https://doi.org/10.1038/s41586-025-09136-2 . This is the Kemp enzyme-design study. Supplementary Table 3 gives computed means/replica SD and experimental apparent barriers; the latter are derived from kcat and are a different observable.
2. **Nature simulation inputs:** https://doi.org/10.5281/zenodo.14563437 . The user supplied EVB_Simulations.zip. Original inputs/hashes are in https://github.com/ND7996/EVB_MHETase/tree/main/Kemp_EVB/inputs/source_manifest.json. The workshop includes an unchanged copy of the published FEP mapping; numerical atom identities override stale comments.
3. **Gutierrez-Rus et al.** Enzyme Enhancement Through Computational Stability Design Targeting NMR-Determined Catalytic Hotspots. JACS (2025). https://doi.org/10.1021/jacs.4c09428 . The Nature article cites this GNCA Kemp-eliminase protocol; the enzyme scaffold differs.
4. **GNCA protocol archive:** https://doi.org/10.5281/zenodo.12666978 . The supplied Kemp_variants.zip contains qfep.inp, with constant H12=100.4 and native state-2 alpha −110 kcal/mol. The archive MD5 matches the supplied listing. The dataset listing specifies CC BY 4.0. Exact selected source bytes, attribution and manifest are in https://github.com/ND7996/EVB_MHETase/tree/main/Kemp_EVB/inputs/protocol_reference/. The archive README calls qfep.inp production input, but it is an analysis input.
5. **Q6 source:** https://github.com/qusers/Q6 , topology-header revision f04f99ed0722be0611f22dd963263bcc5bc78290. qbond establishes the −De Morse convention; qfep.f90 establishes the analysis-input and per-frame mixing/bin-averaging interpretation. No full native binary-QFEP parity calculation is claimed.
6. **EVBOpenMM:** inspected/tested API revision dedf301eb9b0021a4be200327cea53f5b7f034c6 from https://github.com/Martin-Floor/EVBOpenMM . Its documented lack of Q6 LRF prevents an exact-boundary equivalence claim. The reduced core is checked against its two-state energy calculation; the full workflow uses its EVB force scaffold and mapping energies.
7. **OpenMM documentation:** https://docs.openmm.org/latest/userguide/ . OpenMM is the programmable MD engine, with CPU/GPU platforms, forces and reporters; it is not itself one force field.
8. **Rosetta design repository:** https://github.com/Fleishman-Lab/denovoKemp , inspected revision c2a74217e57f8e8f1ae620efa551b7e18bf7c6f0. Its design source is not renamed or republished as EVB parameters.

## Included result provenance

- The browser parameter scan compares saved OpenMM energy data under different coupling and product-offset settings; it is a sensitivity exercise and does not launch new trajectories.
- ../data/short_run/replica_101_raw.json is the user-supplied completed short OpenMM adaptation, 21 windows × 50 records, 6,330 particles. It is not an authors' native .en dataset. The source hash is verified against the recorded comparison; original metadata are preserved in the file and ../data/PROVENANCE.json.
- The source profile applies recovered archived parameters to that short dataset, with 18 coarse bins and minimum 20 per window/bin. It does not establish physical transfer or convergence.
- The parameter viewer includes 25 parameter combinations (24 grid pairs plus the archived baseline) and five examples adjusted to the 13.2 kcal/mol barrier reference. Target matches are workshop illustrations, not independent predictions or recovered author parameter values.
- The optional saved full-enzyme DCD is a separate 480-step engineering run, not a reaction movie or equilibrated production.

Workshop diagrams, comparison charts and viewers are newly authored teaching material. Complete unaltered source pages are reproduced separately in ../assets/paper/ for noncommercial teaching, with attribution, license link and page/source hashes. The Nature article has a CC BY-NC-ND notice; the workshop uses attributed numerical facts and original illustrations. Original parameter-input bytes are retained with provenance. See the project manifests for the complete input and result hashes.

## Calibration lineage checked in November2026

Hong et al., *The evolution of multiple active site configurations in a designed enzyme*, Nature Communications9,3900(2018), [doi:10.1038/s41467-018-06305-y](https://doi.org/10.1038/s41467-018-06305-y), [SI](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-06305-y/MediaObjects/41467_2018_6305_MOESM1_ESM.pdf): full reactive Tables10–21; Table10 reports H12=225, alpha=−324 kcal/mol and a21.2 kcal/mol aqueous-reference target.

Risso et al., Chemical Science11,6134–6148(2020), [doi:10.1039/D0SC01935F](https://doi.org/10.1039/D0SC01935F), [accessible full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7407621/): Methods explicitly reuse Hong2018 calibration and Q5.10. The later JACS/Q6 archive updates reactive terms and supplies100.4/−110; these model versions must not be mixed.

See the website’s parameter-examples section and `Kemp_EVB/docs/KEMP_LITERATURE_REVIEW.md` for every paper used, review limitations, reference-naming discrepancy and current versus historical calibration status.

Full review: [Kemp literature review](https://github.com/ND7996/EVB_MHETase/blob/main/Kemp_EVB/docs/KEMP_LITERATURE_REVIEW.md). Native run commands: [Q6 reconstruction](https://github.com/ND7996/EVB_MHETase/blob/main/Kemp_EVB/docs/NATIVE_Q6_REPRODUCTION.md).
