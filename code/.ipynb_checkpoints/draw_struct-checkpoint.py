{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Standard variables loaded, you are good to go!\n"
     ]
    }
   ],
   "source": [
    "import os\n",
    "from dotenv import load_dotenv, find_dotenv\n",
    "from os.path import join, dirname, basename, exists, isdir\n",
    "\n",
    "### Load environmental variables from the project root directory ###\n",
    "# find .env automagically by walking up directories until it's found\n",
    "dotenv_path = find_dotenv()\n",
    "\n",
    "# load up the entries as environment variables\n",
    "load_dotenv(dotenv_path)\n",
    "\n",
    "# now you can get the variables using their names\n",
    "\n",
    "# Check whether a network drive has been specified\n",
    "DATABASE = os.environ.get(\"NETWORK_URL\")\n",
    "if DATABASE == 'None':\n",
    "    pass\n",
    "else:\n",
    "    pass\n",
    "    #mount network drive here\n",
    "\n",
    "# set up directory paths\n",
    "CURRENT_DIR = os.getcwd()\n",
    "PROJ = dirname(dotenv_path) # project root directory\n",
    "\n",
    "DATA = join(PROJ, 'data') #data directory\n",
    "RAW_EXTERNAL = join(DATA, 'raw_external') # external data raw directory\n",
    "RAW_INTERNAL = join(DATA, 'raw_internal') # internal data raw directory\n",
    "INTERMEDIATE = join(DATA, 'intermediate') # intermediate data directory\n",
    "FINAL = join(DATA, 'final') # final data directory\n",
    "\n",
    "RESULTS = join(PROJ, 'results') # output directory\n",
    "FIGURES = join(RESULTS, 'figures') # figure output directory\n",
    "PICTURES = join(RESULTS, 'pictures') # picture output directory\n",
    "\n",
    "\n",
    "# make folders specific for certain data\n",
    "folder_name = 'BRENDA_data_2019_1'\n",
    "if folder_name != '':\n",
    "    #make folders if they don't exist\n",
    "    if not exists(join(RAW_EXTERNAL, folder_name)):\n",
    "        os.makedirs(join(RAW_EXTERNAL, folder_name))\n",
    "\n",
    "    if not exists(join(INTERMEDIATE, folder_name)):\n",
    "        os.makedirs(join(INTERMEDIATE, folder_name))\n",
    "\n",
    "    if not exists(join(FINAL, folder_name)):\n",
    "        os.makedirs(join(FINAL, folder_name))\n",
    "\n",
    "print('Standard variables loaded, you are good to go!')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "metadata": {},
   "outputs": [],
   "source": [
    "from rdkit import Chem\n",
    "from rdkit.Chem import Draw\n",
    "\n",
    "subswithsmiles = {'L-lactate' : \"C([C@@H](O)C)(=O)[O-]\", 'Glucose' : 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO', 'Glutamate' : \"N[C@@H](CCC(=O)[O-])C(=O)[O-]\"}\n",
    "\n",
    "def draw_substrate_structure(SubstratesWSMILES):\n",
    "    \"\"\"Takes a dictionary of substrate names and their SMILES and outputs .png files containing drawings of\n",
    "       their molecular structure.\"\"\"\n",
    "    \n",
    "    for key in SubstratesWSMILES:\n",
    "        struct = Chem.MolFromSmiles(substrateswithSMILES[key])\n",
    "        name = key\n",
    "        Draw.MolToFile(struct, str(name) + \".png\",)\n",
    "        \n",
    "draw_substrate_structure(subswithsmiles)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 39,
   "metadata": {},
   "outputs": [],
   "source": [
    "from rdkit import Chem\n",
    "from rdkit.Chem import Draw\n",
    "\n",
    "# Placeholder dictionary containing sampling substrates.\n",
    "SomeECNum = {'L-lactate' : \"C([C@@H](O)C)(=O)[O-]\", 'Glucose' : 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO', 'Glutamate' : \"N[C@@H](CCC(=O)[O-])C(=O)[O-]\"}\n",
    "\n",
    "# In completed pipeline, suggest to amend this to \"draw_substrate_structures(SubstratesWSMILES, filename)\"\n",
    "def draw_substrate_structures(SubstratesWSMILES):   \n",
    "    \"\"\"Takes as input a dictionary of substrate names and their SMILES representation, and returns a single .png file containing \n",
    "       drawings of all their structures.\"\"\"\n",
    "    \n",
    "    name = []    # Initiate list containing names.\n",
    "    struct = []  # Initiate list containing mols.\n",
    "\n",
    "    for key in SubstratesWSMILES:\n",
    "        name.append(key)\n",
    "        struct.append(Chem.MolFromSmiles(SubstratesWSMILES[key]))  # Calculate mol from SMILES and append to struct-list.\n",
    "    \n",
    "    # Vary row-length according to how many substrates are in dictionary. Modify these values according to preference. \n",
    "    if len(SubstratesWSMILES) < 4:\n",
    "        rowlength = 3\n",
    "    elif len(SubstratesWSMILES) == 4:\n",
    "        rowlength = 2\n",
    "    else:\n",
    "        rowlength = 4\n",
    "    \n",
    "    # Draw grid of molecules and their respective legends. \n",
    "    img=Draw.MolsToGridImage(struct, molsPerRow = rowlength, subImgSize=(300, 300), legends = name)\n",
    "    img.save('SomeECNum.png') # In completed pipeline, suggest to amend this to \"img.save(filename + \"substrates\")\"\n",
    "\n",
    "draw_substrate_structures(SomeECNum)\n",
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
