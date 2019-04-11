# This will be an installable cheminformatics package that pulls together everything. The "pipeline"


## Install required packages

Install PubChemPy
```
conda install -c bioconda pubchempy
```

Install rdkit in the anaconda environment (https://www.rdkit.org/docs/Install.html)
```
conda install -c conda-forge rdkit
```

You may also need to install gcc
```
conda install libgcc
```


Install CIRpy
```
pip install cirpy
```


## brenda module for loading BRENDA substrate data
