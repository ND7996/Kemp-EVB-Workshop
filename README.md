# Kemp EVB workshop

Open the [public workshop website](https://nd7996.github.io/Kemp-EVB-Workshop/?v=portable-setup-9060d6f).

The website introduces EVB and spherical boundaries, explains the Kemp mechanism, compares Q6 and OpenMM, and provides an interactive α/H₁₂ scan.

## Portable OpenMM setup

No personal paths or separate EVBOpenMM checkout are required. From this repository's root:

```bash
conda env create -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python -m pip install --no-deps -r kemp_evb/calculation/requirements-evb.txt
python kemp_evb/calculation/scripts/run_short_evb_demo.py --check
python kemp_evb/calculation/scripts/run_short_evb_demo.py
```

If `kemp-short` already exists, use `conda env update -f kemp_evb/calculation/environment.yml` instead of environment creation, then activate it and run the pip command. The requirements file installs the pinned Martin-Floor EVBOpenMM package. The short example is a reduced potential-energy teaching scan, not a full-enzyme FEP calculation.

[Detailed installation and troubleshooting](kemp_evb/calculation/README.md) · [Q6 and OpenMM inputs and outputs](kemp_evb/README.md)

If an older page appears, use the website link above and refresh with Ctrl+Shift+R (Cmd+Shift+R on macOS).
