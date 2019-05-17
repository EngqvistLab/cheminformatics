# Description of cheminformatics library
The aim of this package is to easily obtain information the properties of small molecules. **For example usage please refer to the "package_usage.ipynb" file.**


## Installation
Download repository and unzip (alternatively fork or clone), cd to the project base folder and execute the command below:

```
pip3 install -e .
```

If using an anaconda environment you may have to first locate the anaconda pip using whereis.
```
whereis pip
```

Locate the appropriate file path (the one that has anaconda and the correct environment in the filepath) and run the modified command. For example:

```
/home/username/anaconda3/envs/py37/bin/pip install -e .
```

The library should now be available for loading in all your python scripts.


## Requirements
This library relies on several cheminformatics and machine learning libraries

numpy

pandas

matplotlib

PIL

rdkit

cirpy

pubchempy

umap-learn

scikit-learn


The installation instructions below assume that you are using a conda environment as provided by Anaconda (https://www.anaconda.com/).


Install rdkit
```
conda install -c conda-forge rdkit
```


You may also need to install gcc
```
conda install libgcc
```


Install CIRpy (here you might again need to use the "whereis" trick as outlined above)
```
pip install cirpy
```


Install PubChemPy
```
conda install -c bioconda pubchempy
```


Install SciKit-learn
```
conda install scikit-learn
```


Install UMAP-learn
```
conda install -c conda-forge umap-learn
```



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
