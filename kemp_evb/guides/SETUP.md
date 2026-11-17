# Linux setup for the short example

Install Git and Conda or Miniforge. From a terminal:

~~~bash
git clone https://github.com/ND7996/Kemp-EVB-Workshop.git
cd Kemp-EVB-Workshop
conda env create -f kemp_evb/calculation/environment.yml
conda activate kemp-short
python kemp_evb/calculation/scripts/run_short_evb_demo.py --check
python kemp_evb/calculation/scripts/run_short_evb_demo.py
~~~

The run saves a browser plot and JSON summary under kemp_evb/calculation/runs/. A CPU is enough.
