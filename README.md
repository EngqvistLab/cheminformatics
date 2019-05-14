# Description of cheminformatics library
The aim of this package is to easily obtain information the properties of small molecules.


## Installation
Download repository and unzip (alternatively fork or clone), cd to the project base folder and execute the command below:

```
pip3 install -e .
```

If using an anaconda environment you may have to first locate the anaconda pip using whereis.
```
whereis pip
```

Locate the appropriate file path and run the modified command. For example:

```
/home/user/anaconda3/envs/py37/bin/pip install -e .
```

The library should now be available for loading in all your python scripts.


## Requirements
rdkit

cirpy

pubchempy


The installation instructions assume that you are using a conda environment.


Install rdkit (https://www.rdkit.org/docs/Install.html)
```
conda install -c conda-forge rdkit
```

Documentation: https://www.rdkit.org/docs/index.html

Getting started: https://www.rdkit.org/docs/GettingStartedInPython.html



You may also need to install gcc
```
conda install libgcc
```



Install CIRpy
```
pip install cirpy
```

Getting started: https://cirpy.readthedocs.io/en/latest/guide/gettingstarted.html


https://cirpy.readthedocs.io/en/latest/guide/resolvers.html


Install PubChemPy
```
conda install -c bioconda pubchempy
```
