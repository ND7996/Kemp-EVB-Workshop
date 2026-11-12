# Kemp EVB workshop

**Public webpage:** [Kemp EVB workshop](https://nd7996.github.io/Kemp-EVB-Workshop/). The site is deployed from this workshop-only repository using GitHub Pages. The first deployment may take a few minutes after the workflow completes.


**Open the local `index.html` file** in Edge, Chrome or Firefox. It introduces the paper and reaction, then guides you through EVB, the OpenMM adaptation and the hands-on exercises. The session is a 30-minute introduction followed by a 30-minute practical. You can use the slides and saved data without installing software or running a simulation. Keep the folder together so its figures and interactive exercises load.


- `index.html`: paper, atom names, mechanism, preparation, relaxation, FEP and Q6/OpenMM comparison, with source links.
- `presentation.html`: browser slides; opens from the workshop page.
- `Kemp_EVB_workshop.pdf`: presentation PDF.
- `PARTICIPANT_HANDOUT.pdf`: 30-minute practical worksheet.
- `guides/`: setup, practical, instructor notes, speaker notes and attribution.
- `assets/`: figures and attributed source-page images.
- `viewers/`: three offline interactive exercises.
- `data/`: minimal inputs and saved energies for the optional practical.

For optional live recomputation, follow [setup](guides/SETUP.md), activate the environment, then run from this folder:

```bash
python run_practical.py --step core
python run_practical.py --step profile
python run_practical.py --step scan
```

`start_workshop.py` is an optional launcher, and `workshop_manifest.json` checks the workshop files. The practical uses saved data, so no long MD run is needed. Curves adjusted to a chosen target are exercises in parameter effects, not validated predictions.

The browser slides and PDF are included in this folder. Workshop build tools live separately in `Kemp_EVB/workshop_tools/`.
