# Optional short EVB calculation

This Linux-friendly example uses OpenMM to evaluate an EVB model for five mapped atoms along a prescribed geometry path. It checks the lower EVB energy against EVBOpenMM. It is a software exercise, not molecular dynamics or a free-energy barrier calculation; the coupling and offset are teaching values.

## Included files

- scripts/run_short_evb_demo.py: input check and run command.
- scripts/evb_teaching.py and scripts/kemp_inputs.py: reduced model and source reader.
- inputs/q6/Des27_7_IN/: mapped coordinates, topology, and reaction file.
- environment.yml: Python, OpenMM, NumPy, and pinned EVBOpenMM.

## Install and run on Linux

Install Git and Conda or Miniforge. In a terminal:

~~~bash
git clone https://github.com/ND7996/Kemp-EVB-Workshop.git
cd Kemp-EVB-Workshop
conda env create -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python kemp_evb/calculation/scripts/run_short_evb_demo.py --check
python kemp_evb/calculation/scripts/run_short_evb_demo.py
~~~

The final command saves teaching_evb.html and teaching_summary.json in a new folder under kemp_evb/calculation/runs/. Open the HTML file in a browser. A CPU is enough; no GPU is needed. To explore a different teaching setting, add --coupling 12 --offset 2 to the final command. These values are not fitted parameters.

The full-enzyme preparation, relaxation, FEP sampling, and analysis workflow is maintained in the [Kemp backend](https://github.com/ND7996/EVB_MHETase/tree/main/Kemp_EVB).
