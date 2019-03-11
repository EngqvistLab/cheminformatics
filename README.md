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
    |  |- template_engine_python.py # tool for generating new Python scripts containing boilerplate code
    |  |- template_engine_R.py       # tool for generating new R scripts containing boilerplate code
    |  +- template_engine_Rmd.py       # tool for generating new R-markdown scripts containing boilerplate code
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

Collect all BRENDA data -- Done \
Split BRENDA file into single EC files -- Not started \
Extract all unique substrates from a single EC number -- Not started \
Figure out how to filter out co-factors -- Not started \
Convert substrate names to Smiles/InChi -- Not started \
Calculate fingerprints from Smiles/InChi -- Not started \
Implement different algorithms for comparing fingerprints -- Not started \
Visualize the (dis)similarity of substrates within an ec class (t-SNE, subgraphs etc.) -- Not started \
Match up the inputs and outputs between the different steps to make a complete pipeline -- Not started \
Run pipeline on all of BRENDA -- Not started \
Use visualization tools to explore data -- Not started \
