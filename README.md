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
* You will need the dotenv package for this to work.
```
# do
conda install -c conda-forge python-dotenv

# OR
pip install python-dotenv
```


Initial goals and their status
----------

Collect all BRENDA data -- **Done, Martin** \
Split BRENDA file into single EC files -- **Done, Oskar** \
Extract all unique substrates from a single EC number -- **Done, Martin** \
Figure out how to filter out co-factors -- **Done, Martin** \
Convert substrate names to Smiles/InChi -- **Done, Emma, Rasmus, David** \
Calculate fingerprints from Smiles/InChi -- **Done, Emma, Rasmus, David** \
Implement different algorithms for comparing fingerprints \
......Kernels -- **In progress, David** \
......Graph algorithms to find common subgraphs -- **Done, David, Isabella** \
Visualize the (dis)similarity of substrates within an ec class (t-SNE, subgraphs etc.) -- **In progress, Emma, Rasmus**




More complex goals and their status
----------

To support the final analysis and to have computational flexibility we need functions for the following:



**Single molecule:** \
name -> smiles \
smiles -> fingerprint (should be able to choose which fingerprint) \
smiles -> draw molecule structure



**Two or more molecules in a list:** \
list of names -> list of smiles \
list of smiles -> list of fingerprints (should be able to choose which fingerprint) \
list of smiles -> draw all molecule structures


list of two fingerprints -> compute similarity score between the pair of molecules \
list of many fingerprints -> compute all vs all similarity scores


list of two fingerprints -> compute whether common substructure exists in the pair \
list of many fingerprints -> compute whether common substructure exists between all possible pairs


visualize common substructure of two molecules -- **Done, Jonatan**



**I also want these scripts to support our work in the lab. There we select a subset of substrates from each EC to test**
* Given a set of molecules, with m molecules already selected, select one molecule such that it maximizes the chemical diversity in the selection m + n -- **Done, David**
* Given a set of molecules, make a selection of n molecules that maximizes the chemical diversity in the selection -- **Done, David**



**Before we can analyze EC number substrate diversity we need to remove molecules that are co-factors in the reaction**
* For each EC, filter out substrates matching a pre-defined list. -- **Done, Martin**
* On the remainder calculate all pairwise Tanimoto scores. -- **Done, Emma, Rasmus**
* For each molecule calculate the sum of tanimoto scores, the median tanimoto score, the minimum tanimoto score and the maximum tanimoto score. -- **Done, Emma, Rasmus**
* Use the calculated values to flag substrates that are very different from the others. Manual curation may have to follow. -- **Not started**
* Get some "common sense" descriptors for molecules, such as mw, polarity, etc. -- **In progress, Oskar**


Final goals and their status
------------

Match up the inputs and outputs between the different steps to make a complete pipeline -- **In progress, Martin** \
Add in code to cache partial files and to look for these before starting a new computation -- **In progress, Martin** \
Run pipeline on all of BRENDA -- **Not started** \
Use visualization tools to explore data, make graphs -- **Not started**




Resources for RDKit
------------

Install rdkit in the anaconda environment (https://www.rdkit.org/docs/Install.html)

```
conda install -c conda-forge rdkit
```

You may also need to install gcc

```
conda install libgcc
```

Documentation: https://www.rdkit.org/docs/index.html

Getting started: https://www.rdkit.org/docs/GettingStartedInPython.html


Resources for CIRpy
------------

```
pip install cirpy
```

Getting started: https://cirpy.readthedocs.io/en/latest/guide/gettingstarted.html


https://cirpy.readthedocs.io/en/latest/guide/resolvers.html
