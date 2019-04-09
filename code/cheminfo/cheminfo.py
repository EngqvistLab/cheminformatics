#!/usr/bin/env python3
"""
Modify this line to briefly discribe the functionality of ./cheminfo/cheminfo.py

Copyright (C) 2019  Martin Engqvist Lab
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
from os.path import join, dirname, basename, exists, isdir

import cirpy

import json
import re
import datetime

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdFMCS, Draw, AllChem
from rdkit.Chem import Descriptors
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker



class NameToSmile(object):
	'''
	Class that makes use of cirpy to get smiles from substrate names.
	Makes use of caching to save computation.
	Since the chached file is modified with each new substrate I cache the old files,
	just in case something gets corrupted.
	'''
	def __init__(self, names, retest_none=False):
		'''
		Initialize by passing a list of input names.
		These will be converted to lowercase.
		'''
		assert type(names) in [list, set], 'Error, the input must be a list or set'
		assert all([type(s) is str for s in names]), 'Error, each item in the input list must be a string'
		self.input_names = [s.lower() for s in names]
		self.retest_none = retest_none

		# setup variables related to the cached file filename
		self.currentDT = datetime.datetime.now()
		self.cached_file_dir = './data/'
		self.fileend = 'substrate_cache.json'
		self.date = '%s-%s-%s' % (self.currentDT.year, str(self.currentDT.month).rjust(2, '0'), str(self.currentDT.day).rjust(2, '0')) # YY-MM-DD
		self.todays_filename = join(self.cached_file_dir, '%s_%s' % (self.date, self.fileend))

		# get the most recent cached file and load data
		self.most_recent_file = self._get_most_recent_filename()
		self.smile_data = self._load_data()

		# get any missing smiles
		self._get_missing_smiles()


	def _get_most_recent_filename(self):
	    '''
	    Return the filepath of the most recent cached file.
	    '''
	    # first make a list of all available files
	    file_data = {}
	    for fi in os.listdir(self.cached_file_dir):
	        if fi.endswith(self.fileend):
	            year, month, day = fi.replace(self.fileend, '').split('-')

	            if file_data.get(year) is None:
	                file_data[year] = {}

	            if file_data[year].get(month) is None:
	                file_data[year][month] = {}

	            file_data[year][month][day] = join(self.cached_file_dir, fi)

	    # if there is no file
	    if file_data == {}:
	        return None

	    max_year = max(file_data.keys())
	    max_month = max(file_data[max_year].keys())
	    max_day = max(file_data[max_year][max_month].keys())

	    return file_data[max_year][max_month][max_day]


	def _open_file(self, filepath):
	    '''
	    Open up the json file and return data structure
	    '''
	    with open(filepath, 'r') as f:
	        data = json.loads(f.read())
	    return data


	def _save_data(self, data):
	    '''
	    Save data to todays file.
	    '''
	    with open(self.todays_filename, 'w') as f:
	        f.write(json.dumps(data))


	def _load_data(self):
	    '''
	    Load up existing substrate to smile data
	    '''
	    # if a file from today exists open it up
	    if exists(self.todays_filename):
	        return self._open_file(self.todays_filename)

	    # if no file exists from today then make a copy of the most recent file
	    else:

	        if self.most_recent_file is None:
	            data = {}

	        else:
	            data = self._open_file(self.most_recent_file) # get most recent data

	        self._save_data(data) # save into todays filename
	        return data


	def _get_missing_smiles(self):
		'''
		Use cirpy to obtain smiles for molecules that don't have it
		'''
		counter = 0
		for name in self.input_names:
			if name not in self.smile_data.keys():
				counter += 1

			elif self.retest_none is False and name in self.smile_data.keys():
				continue

			elif self.retest_none is True and name in self.smile_data.keys() and self.smile_data.get(name) is not None:
				continue

			elif self.retest_none is True and name in self.smile_data.keys() and self.smile_data.get(name) is None:
				counter += 1

			else:
				raise ValueError

		if counter != 0:
			print('%s name to smile conversions need to be carried out' % counter)


		counter = 0
		for name in self.input_names:

			if name not in self.smile_data.keys():
				counter += 1
				result = cirpy.resolve(name, 'smiles')
				if result is None:
					result =
				self.smile_data[name] = result

			elif self.retest_none is True and name in self.smile_data.keys() and self.smile_data.get(name) is None:
				counter += 1
				result = cirpy.resolve(name, 'smiles')
				if result is None:
					result =
				self.smile_data[name] = result

			if counter % 100 == 0:
				self._save_data(self.smile_data)

		self._save_data(self.smile_data)


	def names(self, exclude_none=False):
		'''
		Return a list of the input molecule names.
		If exclude_none is set to False all names that don't have a smile will be skipped.
		'''
		if exclude_none is True:
			return self.input_names

		elif exclude_none is False:
			return [s for s in self.input_names if self.smile_data[s] is not None]

		else:
			raise ValueError


	def smiles(self, exclude_none=False):
		'''
		Return a list of the output molecule smiles.
		If exclude_none is set to False all names that don't have a smile will be skipped.
		'''
		return [self.smile_data[s] for s in self.names(exclude_none)]


	def data(self, exclude_none=False):
		'''
		Return a dictionary with substrate name keys and smile values.
		If exclude_none is set to False all names that don't have a smile will be skipped.
		'''
		return {k:v for k, v in zip(self.names(exclude_none), self.smiles(exclude_none))}


#subswithsmiles = {'Massa ringar' : 'COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21','L-lactate' : "C([C@@H](O)C)(=O)[O-]", 'Glucose' : 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO', 'Glutamate' : "N[C@@H](CCC(=O)[O-])C(=O)[O-]"}


class SmileToData(object):
	'''
	A class to convert smiles to usable outputs.
	These include similarity measures, molecular properties, scatter plots etc.
	'''
	def __init__(self, names, smiles, descriptor='morgan3', metric='tanimoto'):
		assert type(names) in [list, set], 'Error, the input names must be a list or set'
		assert all([type(s) is str for s in names]), 'Error, each item in the input list names must be a string'
		assert type(smiles) in [list, set], 'Error, the input smiles must be a list or set'
		assert all([type(s) is str for s in smiles]), 'Error, each item in the input list smiles must be a string'

		self.names = names
		self.smiles = smiles
		self.data = {k:v for k, v in zip(self.names, self.smiles)}

		self.descriptors = {
		    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
		    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3),
		    'morgan5':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,5),
    		'rdkit':       lambda m: FingerprintMols.FingerprintMol(m)
			}

		self.metrics = {
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

		assert descriptor in self.descriptors.keys(), 'Error, the argument discriptor must be one of: %s' % ', '.join(self.descriptors.keys())
		assert metric in self.metrics.keys(), 'Error, the argument metri must be one of: %s' % ', '.join(self.metrics.keys())
		self.descriptor = descriptor
		self.metric = metric

		# make a list of mols
		self.mols = [Chem.MolFromSmiles(s) for s in self.smiles]

		# make a list of fingerprints
		self.fingerprints = [self.descriptors[self.descriptor](s) for s in self.mols]

		# compute the similarity matrix
		self.similarity_matrix = self._compute_similarity_matrix()

		# compute the similarity stats
		self.similarity_stats = self._compute_similarity_stats()


	def _compute_similarity_matrix(self):
		"""
		Returns a data frame with pairwise similarity scores for a list of molecules,
		specified by a list of names and a list of corresponding smiles formats given
		as input parameters. The fingerprints and similarity coefficients can be chosen
		from the list of descriptors and metrics (default 'rdkit' and 'tanimoto').
		"""
		# all against all similarity matrix
		S = []
		score = self.metrics[self.metric] # similarity score
		for fp in self.fingerprints:
			S.append(score(fp, self.fingerprints))

		return pd.DataFrame(S, index=self.names, columns=self.names)


	def _compute_similarity_stats(self):
		'''
		'''
		return None


	def valid_descriptors(self):
		'''
		Return a list of all valid values for descriptors
		'''
		return list(self.descriptors.keys())


	def valid_metrics(self):
		'''
		Return a list of all valid values for metrics.
		'''
		return list(self.metrics.keys())


	def diversity_pick(self, n, firstpicks=[]):
		"""
		Picks a maximally diverse subset of a sert of molecules using the RDKit
		MaxMinPicker. Optionally, a list of names with already
		chosen molecules can be specified.
		"""
		assert type(firstpicks) in [list, set], 'Error, the input firstpicks must be a list or set'
		assert all([type(s) is str for s in firstpicks]), 'Error, each item in the input list firstpicks must be a string'
		assert type(n) is int, 'Error, the input n must be an integer'

		firstpicks = [s.lower() for s in firstpicks]
		assert all([s in self.names for s in firstpicks]), 'Error, not all firstpicks are part of the molecule list'
		assert n < len(self.names) - len(firstpicks), 'Error, you have specified an n that is greater or equal to the available molecule number'

		# get indices of already picked molecules
		ind = []
		for x in firstpicks:
			ind.append(self.names.index(x)) # indices of picked molecules

		# compute all pairwise similarity scores
		ds = []
		score = self.metrics[self.metric]
		for i in range(1, len(self.fingerprints)):
			ds.extend(score(self.fingerprints[i], self.fingerprints[:i], returnDistance=True))

		# make the selection (returns indeces)
		ids = MaxMinPicker().Pick(np.array(ds), len(self.fingerprints), n, ind)

		return [self.names[s] for s in ids]


	def pair_similarity(self, molecules):
		'''
		Return similarity measure between a list of specified molecules.
		These have to be part of the original input molecules.
		'''
		assert type(molecules) in [list, set], 'Error, the input molecules must be a list or set'
		assert all([type(s) is str for s in molecules]), 'Error, each item in the input list molecules must be a string'
		assert all([s in self.names for s in molecules]), 'Error, not all molecules are part of the molecule list'

		pass


	def similarity(self):
		'''
		Return similarity matrix with all similarities
		'''
		return self.similarity_matrix


	def similarity_stats(self):
		'''
		Return the similarity statistics
		'''
		return self.similarity_stats


	def plotMCS(self, names, smiles, molRows=5, subSize=200):
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


	def draw_substrate_structures(self, SubstratesWSMILES):
		"""Takes as input a dictionary of substrate names and their SMILES representation, and returns a single .png file containing
		drawings of all their structures.
		"""

		name = []    # Initiate list containing names.
		struct = []  # Initiate list containing mols.

		for key in SubstratesWSMILES:
			name.append(key)
			struct.append(Chem.MolFromSmiles(SubstratesWSMILES[key]))  # Calculate mol from SMILES and append to struct-list.

		# Vary row-length according to how many substrates are in dictionary. Modify these values according to preference.
		if len(SubstratesWSMILES) < 4:
			rowlength = 3
		elif len(SubstratesWSMILES) == 4:
			rowlength = 2
		else:
			rowlength = 4

		# Draw grid of molecules and their respective legends.
		img=Draw.MolsToGridImage(struct, molsPerRow = rowlength, subImgSize=(300, 300), legends = name)
		img.save('SomeECNum.png') # In completed pipeline, suggest to amend this to "img.save(filename + "substrates")"



	def molecular_weight(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES, converts to mol-format and returns a dictionary of substrate names and
	    and their respective molecular weight."""
	    molweight = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        molweight[key] = Descriptors.ExactMolWt(mol)
	    return molweight


	def polar_surface_area(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES and returns a dictionary of substrate names and
	       their respective polar surface area."""
	    TPSA = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        TPSA[key] = Descriptors.TPSA(mol)
	    return TPSA

	# Kanske att föredra att dela upp den här i två funktioner?

	def Hbond_donors_acceptors(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES and returns a list of dictionaries of substrate names and
	       their respective number of hydrogen bond acceptors and donors. Acceptors on index 0, donors on index 1."""
	    acceptors = {}
	    donors = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        acceptors[key] = Descriptors.NumHAcceptors(mol)
	        donors[key] = Descriptors.NumHDonors(mol)
	    return acceptors, donors


	def OHNH_count(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES, converts to mol-format, calculates the number of NH and OH groups, and outputs
	       a dictionary containing the substrate names and their respective counts."""
	    OHNH = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        OHNH[key] = Descriptors.NHOHCount(mol)
	    return OHNH


	def ON_count(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES, converts to mol-format, calculates the number of O and N atoms, and outputs
	        a dictionary containing the substrate names and their respective counts. """
	    ON = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        ON[key] = Descriptors.NOCount(mol)
	    return ON


	def HeteroAtomCount(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES, converts to mol-format, calculates the number of
	        hetero atoms and returns a dictionary containing the substrate names and their respective hetero counts."""
	    HA = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        HA[key] = Descriptors.NumHeteroatoms(mol)
	    return HA


	def RotaBondCount(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES, converts to mol-format, calculates the number of
	        rotatable bonds and returns a dictionary containing the substrate names and their respective rotatable bond counts."""
	    RB = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        RB[key] = Descriptors.NumRotatableBonds(mol)
	    return RB


	def AromRingCount(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES, converts to mol-format, calculates the number of
	        aromatic rings and returns a dictionary containing the substrate names and their respective aromatic ring counts"""
	    RING = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        RING[key] = Descriptors.NumAromaticRings(mol)
	    return RING


	def AliphRingCount(self, substrate_dict):
	    """Takes as input a dictionary of substrate names and their SMILES, converts to mol-format, calculates the number of
	        aliphatic rings and returns a dictionary containing the substrate names and their respective aliphatic ring counts"""
	    RING2 = {}
	    for key in substrate_dict:
	        mol = Chem.MolFromSmiles(substrate_dict[key])
	        RING2[key] = Descriptors.NumAliphaticRings(mol)
	    return RING2


# do MDS on distances
# do PCA on data
# do t-sne on distances
# do t-sne on data



class Plotting(object):
	'''
	An object holding various plotting setups for the molecules.
	'''
	def __init__(self):
		pass
