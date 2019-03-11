BRENDA Cheminformatics project
=======

A project for looking at substrate differences for enzyme classes.

Overview
--------

    project
    |- doc/             # documentation for the study and other explanatory material
    |  +- paper/        # manuscript(s), whether generated or not
    |
    |- data/            # raw and primary data, are not changed once created
    |  |- raw_internal/ # raw data generated in-lab or by collaborators, will not be altered
    |  |- raw_external/ # data from third-party sources, databases etc, will not be altered
    |  |- intermediate/ # intermediate data on its way from raw to final, are not changed once created
    |  +- final/        # final data for figures/visualizations, will not be altered once created
    |
    |- code/            # all programmatic code relating to the project
    |  |- template_engine_python.py           # tool for generating new Python scripts containing boilerplate code
    |  |- template_engine_jupyter-notebook.py # tool for generating new Python scripts containing boilerplate code
    |  |- template_engine_R.py                # tool for generating new R scripts containing boilerplate code
    |  +- template_engine_R-notebook.py       # tool for generating new R-markdown scripts containing boilerplate code
    |
    |- results          # all output from workflows and analyses
    |  |- figures/      # graphs, likely designated for manuscript figures
    |  +- pictures/     # diagrams, images, and other non-graph graphics
    |
    |- .env             # a file to store passwords and usernames needed for the scripts. This will not get synced to GitHub.
    |- notebook.txt     # a lab notebook where activities relating to this project should be entered
    |- requirements.txt # the requirements file for reproducing the analysis environment,
    |                     e.g. generated with `pip freeze > requirements.txt`
    |- scratch/         # temporary files that can be safely deleted or lost
    |
    |- README.md        # the top level description of content
    +- datapackage.json # metadata for the (input and output) data files

How to use
----------

* Create a new directory for your project.
* Download the [latest version] of this repository, and unzip it in the directory you just created.
* **Create a .env file using a basic text editor (keep it empty if you like) and place in the project base directory (where notebook.txt and requirements.txt are). The code templates will not work without it.**


Goals and their status
----------

Collect all BRENDA data -- **Done** \
Split BRENDA file into single EC files -- **Not started** \
Extract all unique substrates from a single EC number -- **Not started** \
Figure out how to filter out co-factors -- **Not started** \
Convert substrate names to Smiles/InChi -- **Not started** \
Calculate fingerprints from Smiles/InChi -- **Not started** \
Implement different algorithms for comparing fingerprints -- **Not started** \
Visualize the (dis)similarity of substrates within an ec class (t-SNE, subgraphs etc.) -- **Not started** \
Match up the inputs and outputs between the different steps to make a complete pipeline -- **Not started** \
Add in code to cache partial files and to look for these before starting a new computation -- **Not started** \
Run pipeline on all of BRENDA -- **Not started** \
Use visualization tools to explore data -- **Not started**


Resources for RDKit
------------

Install rdkit in the anaconda environment (https://www.rdkit.org/docs/Install.html)

```
conda install -c conda-forge rdkit
```

Documentation: https://www.rdkit.org/docs/index.html

Getting started: https://www.rdkit.org/docs/GettingStartedInPython.html




Some starting data
-----------

**Names:** \
2-hydroxydodecanoate \
(S)-lactate \
2-hydroxypalmitate \
Glycolate


**Difficult names:** \
DL-2-hydroxy-4-methylmercaptobutyrate \
an (S)-2-hydroxy carboxylate


**Smiles:** \
[O-]C(=O)C(O)CCCCCCCCCC \
C[C@@H](C(=O)[O-])O \
CCCCCCCCCCCCCCC(C(=O)[O-])O \
C(C(=O)[O-])O


**InChi:** \
InChI=1S/C12H24O3/c1-2-3-4-5-6-7-8-9-10-11(13)12(14)15/h11,13H,2-10H2,1H3,(H,14,15)/p-1 \
InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/p-1/t2-/m0/s1 \
InChI=1S/C16H32O3/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15(17)16(18)19/h15,17H,2-14H2,1H3,(H,18,19)/p-1 \
InChI=1S/C2H4O3/c3-1-2(4)5/h3H,1H2,(H,4,5)/p-1


**Fingerprints:** \
???
