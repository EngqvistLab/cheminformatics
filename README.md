# Description of cheminformatics library
The aim of this package is to easily obtain information the properties of small molecules. **For example usage please refer to the "package_usage.ipynb" file.**


## Installation
The installation instructions below assume that you are using a conda environment as provided by Anaconda (https://www.anaconda.com/). Download repository and unzip (alternatively fork or clone), cd to the project base folder and execute the commands below:

First set up a conda environment and activate.

```bash
conda env create -f environment.yml
conda activate cheminformatics
```

Then install the cheminformatics package.
```bash
pip3 install -e .
```

The library should now be available for loading in all your python scripts.


## Requirements
This library relies on several cheminformatics and machine learning libraries, including: numpy, pandas, matplotlib, PIL, rdkit, cirpy, pubchempy, umap-learn, scikit-learn


## Documentation for cheminformatics packages
For further readiing regarding the installed cheminformatics packages see:


#### RDkit
https://www.rdkit.org/docs/index.html

https://www.rdkit.org/docs/GettingStartedInPython.html


#### CIRpy
https://cirpy.readthedocs.io/en/latest/guide/gettingstarted.html

https://cirpy.readthedocs.io/en/latest/guide/resolvers.html


#### PubChemPy
https://pubchempy.readthedocs.io/en/latest/
