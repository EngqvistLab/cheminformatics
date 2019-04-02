#!/usr/bin/env python3
"""
Modify this line to briefly discribe the functionality of ./cheminfo/cheminfo.py

Copyright (C) 2017  Martin Engqvist Lab
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
from dotenv import load_dotenv, find_dotenv # do 'pip install python-dotenv'
from os.path import join, dirname, basename, exists, isdir

import rdkit
import cirpy

import json
import re

### Load environmental variables from the project root directory ###
# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)




def remove_dl(substrate):
    '''
    Remove D-, L-, DL- and similar from the beginning of substrates
    '''
    substrate = substrate.lower()
    return re.sub('^d-|^l-|^dl-|^\(r\)-|^\(s\)-', '', substrate)


def get_all_ec():
    '''
    Get a list of all ec numbers
    '''
    return filtered_subst_data.keys()


def get_substrates(ec):
    '''
    Obtain a list of all the natural substrates for an ec number
    '''
    substrate_list = filtered_subst_data[ec]['first_substrate']
    return sorted(list(set([remove_dl(s) for s in substrate_list])))


filepath = join('./data/', 'natural_substrates_filtered.json')
with open(filepath, 'r') as f:
    filtered_subst_data = json.loads(f.read())

x = get_substrates('1.1.3.15')
print(x)


###################


import datetime


currentDT = datetime.datetime.now()
cached_file_dir = join(INTERMEDIATE, 'brenda_data_2019_1')
filename = 'substrate_cache.json'
date = '%s-%s-%s' % (currentDT.year, str(currentDT.month).rjust(2, '0'), str(currentDT.day).rjust(2, '0')) # YY-MM-DD


def todays_filename():
    '''
    Get the filename of todays file.
    '''
    return join(cached_file_dir, '%s_%s' % (date, filename))


def most_recent_filename():
    '''
    Return the filepath of the most recent file.
    '''
    # first make a list of all available files
    file_data = {}
    for fi in os.listdir(cached_file_dir):
        if fi.endswith(filename):
            year, month, day = fi.replace(filename, '').split('-')

            if file_data.get(year) is None:
                file_data[year] = {}

            if file_data[year].get(month) is None:
                file_data[year][month] = {}

            file_data[year][month][day] = join(cached_file_dir, fi)

    # if there is no file
    if file_data == {}:
        return None

    max_year = max(file_data.keys())
    max_month = max(file_data[max_year].keys())
    max_day = max(file_data[max_year][max_month].keys())

    return file_data[max_year][max_month][max_day]


def open_file(filepath):
    '''
    Open up the json file and return data structure
    '''
    with open(filepath, 'r') as f:
        data = json.loads(f.read())
    return data


def save_data(data):
    '''
    Save data to todays file.
    '''
    filepath = todays_filename()
    with open(filepath, 'w') as f:
        f.write(json.dumps(data))


def load_data():
    '''
    Load up existing substrate to smile data
    '''
    # if a file from today exists open it up
    filepath = todays_filename()
    if exists(filepath):
        return open_file(filepath)

    # if no file exists from today then make a copy of the most recent file
    else:
        old_filepath = most_recent_filename()

        if old_filepath is None:
            data = {}

        else:
            data = open_file(old_filepath) # get most recent data

        save_data(data) # save into todays filename
        return data

data = load_data()


###################

import cirpy

def get_brenda_smiles():
	data = load_data()

	counter = 0
	for name in all_substrates:
	    counter += 1

	    if name not in data.keys():
	        data[name] = cirpy.resolve(name, 'smiles')

	    if counter % 100 == 0:
	        print(counter)
	        save_data(data)

	return data


###################


import numpy
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

descriptors = {
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3),
    'morgan5':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,5),
    'rdkit':       lambda m: FingerprintMols.FingerprintMol(m)
}

metrics = {
    'asymmetric':    DataStructs.AsymmetricSimilarity,
    'braunblanquet': DataStructs.BulkBraunBlanquetSimilarity,
    'cosine':        DataStructs.BulkCosineSimilarity,
    'dice':          DataStructs.BulkDiceSimilarity,
    'kulczynski':    DataStructs.BulkKulczynskiSimilarity,
    'mcconnaughey':  DataStructs.BulkMcConnaugheySimilarity,
    'rogotgoldberg': DataStructs.BulkRogotGoldbergSimilarity,
    'russel':        DataStructs.BulkRusselSimilarity,
    'sokal':         DataStructs.BulkSokalSimilarity,
    'tanimoto':      DataStructs.BulkTanimotoSimilarity
}


def diversityPick(names, smiles, ntopick, descriptor='rdkit', metric='tanimoto', firstpicks=[]):
	"""
	Picks a maximally diverse subset of a sert of molecules using the RDKit
	MaxMinPicker. The molecules are given as a list of names and a list with the
	corresponding smiles formats. Optionally, a list of names with already
	chosen molecules can be specified.
	"""
    # make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in smiles]

    # make a list of fingerprints
    fingerprint = descriptors[descriptor] # fingerprint type
    fps = [fingerprint(x) for x in ms]

    ind = []
    for x in firstpicks:
        ind.append(names.index(x)) # indices of picked molecules

    ds = []
    score = metrics[metric] # similarity score
    for i in range(1,len(fps)):
         ds.extend(score(fps[i],fps[:i],returnDistance=True))

    ids = MaxMinPicker().Pick(numpy.array(ds), len(fps), ntopick, ind)
    return ids


import numpy as np, pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols

descriptors = {
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3),
    'morgan5':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,5),
    'rdkit':       lambda m: FingerprintMols.FingerprintMol(m)
}

metrics = {
    'asymmetric':    DataStructs.AsymmetricSimilarity,
    'braunblanquet': DataStructs.BulkBraunBlanquetSimilarity,
    'cosine':        DataStructs.BulkCosineSimilarity,
    'dice':          DataStructs.BulkDiceSimilarity,
    'kulczynski':    DataStructs.BulkKulczynskiSimilarity,
    'mcconnaughey':  DataStructs.BulkMcConnaugheySimilarity,
    'rogotgoldberg': DataStructs.BulkRogotGoldbergSimilarity,
    'russel':        DataStructs.BulkRusselSimilarity,
    'sokal':         DataStructs.BulkSokalSimilarity,
    'tanimoto':      DataStructs.BulkTanimotoSimilarity
}


def similarity(names, smiles, descriptor='rdkit', metric='tanimoto'):
	"""
	Returns a data frame with pairwise similarity scores for a list of molecules,
	specified by a list of names and a list of corresponding smiles formats given
	as input parameters. The fingerprints and similarity coefficients can be chosen
	from the list of descriptors and metrics (default 'rdkit' and 'tanimoto').
	"""
    if descriptor not in descriptors:
        raise ValueError('Invalid descriptor name ' + descriptor)

    if metric not in metrics:
        raise ValueError('Invalid metric ' + metric)

    # make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in smiles]

    # make a list of fingerprints
    fingerprint = descriptors[descriptor] # fingerprint type
    fps = [fingerprint(x) for x in ms]

    S = []
    score = metrics[metric] # similarity score

    # all against all similarity matrix
    for fp in fps: S.append(score(fp, fps))

    return pd.DataFrame(S,index=names,columns=names)


from rdkit.Chem import rdFMCS, Draw, AllChem
from rdkit import Chem





def plotMCS(names, smiles, molRows=5, subSize=200):
	"""
	Plots a given set of molecules in a grid with their maximum common substructure
	highlighted in red. The molecules are given as a list of names and a list with
	the corresponding smiles formats. Optionally, the number of rows in the grid
	and the size of each subfigure can be specified.
	"""
    # make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in smiles]

    # maximum common substructure
    res = rdFMCS.FindMCS(ms)
    mcs = Chem.MolFromSmarts(res.smartsString)

    # align common structure
    AllChem.Compute2DCoords(mcs)
    for m in ms: AllChem.GenerateDepictionMatching2DStructure(m, mcs)

    img = Draw.MolsToGridImage(ms, molsPerRow = molRows, \
            highlightAtomLists = [mol.GetSubstructMatch(mcs) for mol in ms], \
            subImgSize = (subSize, subSize), legends = names)

    img.save('mcs.png')
