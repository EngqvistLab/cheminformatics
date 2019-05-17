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
import json
import re
import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

import cirpy
from pubchempy import Compound, get_compounds

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, MACCSkeys, rdFMCS, Draw, Descriptors, rdmolops
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.Chem.Draw import SimilarityMaps, IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults
from rdkit.ML.Cluster import Butina

from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE, MDS

from cheminformatics.helpfunctions import clean_name

from pkg_resources import resource_stream, resource_filename, resource_exists

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
		self.input_names = [clean_name(s) for s in names if s is not '']
		self.retest_none = retest_none

		# setup variables related to the cached file filename
		self.currentDT = datetime.datetime.now()
		self.cached_file_dir = resource_filename(__name__, 'data')
		self.fileend = 'substrate_cache.json'
		self.date = '%s-%s-%s' % (self.currentDT.year, str(self.currentDT.month).rjust(2, '0'), str(self.currentDT.day).rjust(2, '0')) # YY-MM-DD
		self.todays_filename = join(self.cached_file_dir, '%s_%s' % (self.date, self.fileend))

		# get the most recent cached file and load data
		self.most_recent_file = self._get_most_recent_filename()
		self.smile_data = self._load_data()

		# get any missing smiles
		self._get_missing_smiles()

		# remove duplicates
		self._filter_duplicates()


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

			elif name in self.smile_data.keys() and self.retest_none is False:
				continue

			elif name in self.smile_data.keys()  and self.retest_none is True and self.smile_data.get(name) is not None:
				continue

			elif name in self.smile_data.keys()  and self.retest_none is True and self.smile_data.get(name) is None:
				counter += 1

			else:
				raise ValueError

		if counter != 0:
			print('%s name to smile conversions need to be carried out' % counter)


		counter = 1
		for name in self.input_names:

			if name not in self.smile_data.keys():
				counter += 1

				# Check that compound is listed on PubChem, then use cirpy
				result_list = get_compounds(name, 'name')
				if result_list != []:
					result = result_list[0].canonical_smiles
				else:
					result = cirpy.resolve(name, 'smiles')
				self.smile_data[name] = result

			elif name in self.smile_data.keys() and self.retest_none is True and self.smile_data.get(name) is None:
				counter += 1

				# Check that compound is listed on PubChem, then use cirpy
				result_list = get_compounds(name, 'name')
				if result_list != []:
					result = result_list[0].canonical_smiles
				else:
					result = cirpy.resolve(name, 'smiles')
				self.smile_data[name] = result

			else:
				# skip those that I already have
				continue

			if counter % 100 == 0:
				self._save_data(self.smile_data)
				print('%s done' % counter)

		if counter != 1:
			self._save_data(self.smile_data)


	def _filter_duplicates(self):
		'''
		Filter out cases where different molucule names represent the same smiles string.
		For example, gly and glycine.
		'''
		pass


	def names(self, exclude_none=True):
		'''
		Return a list of the input molecule names.
		If exclude_none is set to True all names that don't have a smile will be skipped.
		'''
		if exclude_none is False:
			return self.input_names

		elif exclude_none is True:
			return [s for s in self.input_names if self.smile_data[s] is not None]

		else:
			raise ValueError


	def smiles(self, exclude_none=True):
		'''
		Return a list of the output molecule smiles.
		If exclude_none is set to True all names that don't have a smile will be skipped.
		'''
		return [self.smile_data[s] for s in self.names(exclude_none)]


	def data_dict(self, exclude_none=True):
		'''
		Return a dictionary with substrate name keys and smile values.
		If exclude_none is set to True all names that don't have a smile will be skipped.
		'''
		return {k:v for k, v in zip(self.names(exclude_none), self.smiles(exclude_none))}


	def data_frame(self, exclude_none=True):
		'''
		Return a dictionary with substrate name keys and smile values.
		If exclude_none is set to True all names that don't have a smile will be skipped.
		'''
		new_dict = {'molecule':[], 'smile':[]}
		for mol, smile in zip(self.names(exclude_none), self.smiles(exclude_none)):
			new_dict['molecule'].append(mol)
			new_dict['smile'].append(smile)

		return pd.DataFrame(new_dict)



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

		self.input_names = names
		self.input_smiles = smiles
		self.smile_data = {k:v for k, v in zip(self.input_names, self.input_smiles)}

		self.descriptors = {
		    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
		    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3),
		    'morgan5':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,5),
    		'rdkit':       lambda m: rdmolops.RDKFingerprint(m, fpSize=1024, tgtDensity=0)
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

		self.properties = {
		    'molwt':{'function':Descriptors.ExactMolWt, 'explanation':'molecular weight'},
			'tpsa':{'function':Descriptors.TPSA, 'explanation':'polar surface area'},
			'num_hbond_donors':{'function':Descriptors.NumHDonors, 'explanation':'number of hydrogen bond donors'},
			'num_hbond_acceptors':{'function':Descriptors.NumHAcceptors, 'explanation':'number of hydrogen bond acceptors'},
			'num_nhoh':{'function':Descriptors.NHOHCount, 'explanation':'number of NH and OH groups'},
			'num_no':{'function':Descriptors.NOCount, 'explanation':'number of N and O atoms'},
			'num_hetatom':{'function':Descriptors.NumHeteroatoms, 'explanation':'number of hetero atoms'},
			'num_rotbond':{'function':Descriptors.NumRotatableBonds, 'explanation':'number of rotatable bonds'},
			'num_arorings':{'function':Descriptors.NumAromaticRings, 'explanation':'number of aromatic rings'},
			'num_alirings':{'function':Descriptors.NumAliphaticRings, 'explanation':'number of aliphatic rings'}
		}


		assert descriptor in self.descriptors.keys(), 'Error, the argument discriptor must be one of: %s' % ', '.join(self.descriptors.keys())
		assert metric in self.metrics.keys(), 'Error, the argument metri must be one of: %s' % ', '.join(self.metrics.keys())
		self.descriptor = descriptor
		self.metric = metric

		# make a list of mols
		self.mol_data = [Chem.MolFromSmiles(s) for s in self.input_smiles]

		# make a list of fingerprints
		self.fingerprint_data = [self.descriptors[self.descriptor](s) for s in self.mol_data]

		# compute the similarity matrix
		self.similarity_matrix = self._compute_similarity_matrix()

		# compute the molecule similarity stats
		self.mol_stats = self._compute_mol_stats()

		# compute the global similarity stats
		self.glob_stats = self._compute_glob_stats()


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
		for fp in self.fingerprint_data:
			S.append(score(fp, self.fingerprint_data))

		return pd.DataFrame(S, index=self.input_names, columns=self.input_names)


	def _compute_glob_stats(self):
		'''
		Compute some basic statistics on the similarity matrix
		'''
		sim_data = self.similarity().copy()

		# get rid of half the matrix and the diagonal
		sim_data.values[np.tril_indices_from(sim_data)] = np.nan
		glob_stats = {'min':sim_data.unstack().min(),
								'max':sim_data.unstack().max(),
								'sum':sim_data.unstack().sum(),
								'median':sim_data.unstack().median(),
								'mean':sim_data.unstack().mean(),
								'stdev':sim_data.unstack().std()}
		return glob_stats


	def _compute_mol_stats(self):
		'''
		Compute some basic statistics on the similarity matrix
		'''
		sim_data = self.similarity().copy()

		# remove the diagonal
		sim_data.values[tuple([np.arange(sim_data.shape[0])]*2)] = np.nan

		# calculate stats
		mol_stats = pd.DataFrame()
		mol_stats['min'] = sim_data.min()
		mol_stats['max'] = sim_data.max()
		mol_stats['sum'] = sim_data.sum()
		mol_stats['median'] = sim_data.median()
		mol_stats['mean'] = sim_data.mean()
		mol_stats['stdev'] = sim_data.std()

		return mol_stats


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


	def names(self):
		'''
		Return a list of the input molecule names.
		'''
		return self.input_names


	def smiles(self):
		'''
		Return a list of the output molecule smiles.
		'''
		return [self.smile_data[s] for s in self.names()]


	def molecules(self, return_only=None):
		'''
		Return a list of molecule objects
		'''
		if return_only is None:
			return self.mol_data
		else:
			assert type(return_only) is str, 'Error, the value in return_only has to be a string.'
			return_only = clean_name(return_only)
			assert return_only in self.names(), 'Error, the submitted molecule is not in the original set that was used to initialize the class.'
			ix = self.names().index(return_only)
			return self.molecules()[ix]


	def fingerprints(self):
		'''
		Return a list of the molecule fingerprints as objects.
		'''
		return self.fingerprint_data


	def fingerprints_str(self):
		'''
		Return a list of the molecule fingerprints as bit strings.
		'''
		return [s.ToBitString() for s in self.fingerprint_data]


	def fingerprints_list(self):
		'''
		Return a list of the molecule fingerprints bit lists.
		'''
		return [[int(b) for b in s] for s in self.fingerprints_str()]


	def property(self, property_type):
		'''
		Get a list of a certain chemical property for all molecules
		'''
		assert property_type in self.valid_properties(), 'Error, "property_type" must be one of: %s' % ', '.join(self.valid_properties())
		return [self.properties[property_type]['function'](s) for s in self.molecules()]


	def valid_properties(self):
		'''
		Return a list of valid properties
		'''
		return sorted(list(self.properties.keys()))


	def explain_properties(self):
		'''
		Return a list of valid properties and their explanation
		'''
		out_data = []
		for prop in sorted(list(self.properties.keys())):
			out_data.append('"%s": %s' % (prop, self.properties[prop]['explanation']))
		return out_data


	def data_dict(self):
		'''
		Return a dictionary containing chemical properties for each molecule
		'''
		data_dict = {}
		for i in range(len(self.input_names)):
			name = self.input_names[i]
			smile = self.input_smiles[i]
			mol = self.molecules()[i]
			data_dict[name] = {}
			data_dict[name]['smile'] = smile

			for prop in sorted(list(self.properties.keys())):
				data_dict[self.input_names[i]][prop] = self.properties[prop]['function'](mol)

		return data_dict


	def data_frame(self):
		'''
		Return a data frame containing chemical properties for each molecule
		'''
		new_dict = {'molecule':self.input_names, 'smile':self.input_smiles}

		prop_dict = {}
		for prop in sorted(list(self.properties.keys())):
			prop_dict[prop] = [self.properties[prop]['function'](s) for s in self.molecules()]

		combined_dict = {**new_dict, **prop_dict}

		return pd.DataFrame(combined_dict)


	def similarity(self):
		'''
		Return similarity matrix with all similarities
		'''
		return self.similarity_matrix


	def distance(self):
		'''
		Return distance matrix with all chemnical distances
		'''
		return 1-self.similarity_matrix


	def molecule_similarity_stats(self):
		'''
		Return similarity statistics for each molecule
		'''
		return self.mol_stats


	def global_similarity_stats(self):
		'''
		Return the global similarity statistics for the entire dataset
		'''
		return self.glob_stats


	def diversity_pick(self, n, firstpicks=[]):
		"""
		Picks a maximally diverse subset of a sert of molecules using the RDKit
		MaxMinPicker. Optionally, a list of names with already
		chosen molecules can be specified.
		"""
		assert type(firstpicks) in [list, set], 'Error, the input firstpicks must be a list or set'
		assert all([type(s) is str for s in firstpicks]), 'Error, each item in the input list firstpicks must be a string'

		assert type(n) is int, 'Error, the input n must be an integer'

		firstpicks = [clean_name(s) for s in firstpicks]
		assert all([s in self.input_names for s in firstpicks]), 'Error, not all firstpicks are part of the molecule list'
		assert n < len(self.input_names) - len(firstpicks), 'Error, you have specified an n that is greater or equal to the available molecule number'

		# get indices of already picked molecules
		ind = []
		for x in firstpicks:
			ind.append(self.input_names.index(x)) # indices of picked molecules

		# compute all pairwise similarity scores
		ds = []
		score = self.metrics[self.metric]
		for i in range(1, len(self.fingerprint_data)):
			ds.extend(score(self.fingerprint_data[i], self.fingerprint_data[:i], returnDistance=True))

		# make the selection (returns indeces)
		ids = MaxMinPicker().Pick(np.array(ds), len(self.fingerprint_data), n, ind)

		return [self.input_names[s] for s in ids]


	def draw_structures(self, highlight_substructure=False):
		'''
		Returns a single .png file containing drawings of all substrate structures.
		'''
		mols = self.molecules()
		img_size = (300, 300)
		num_per_row = int(round(len(mols)**0.5))

		# maximum common substructure
		res = rdFMCS.FindMCS(mols)
		mcs = Chem.MolFromSmarts(res.smartsString)

		# align common structure
		# AllChem.Compute2DCoords(mcs)
		# for m in mols:
		# 	AllChem.GenerateDepictionMatching2DStructure(m, mcs)

		if highlight_substructure is False:
			img = Draw.MolsToGridImage(mols, molsPerRow=num_per_row,
										subImgSize=img_size,
										legends=self.names())
			return img

		elif highlight_substructure is True:
			img = Draw.MolsToGridImage(mols, molsPerRow=num_per_row,
							highlightAtomLists = [mol.GetSubstructMatch(mcs) for mol in mols],
							subImgSize=img_size,
							legends=self.names())
			return img

		else:
			raise ValueError


	def draw_mol_comparison(self, refmol, mol):
		'''
		Compare the structure of two molecules.
		The function takes the molecule names as input.
		'''
		fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(self.molecules(refmol), self.molecules(mol),
						SimilarityMaps.GetMorganFingerprint,
						metric=DataStructs.TanimotoSimilarity)
		return fig, maxweight


	# def draw_substructure_grid(self):
	# 	'''
	# 	'''
	# 	import io
	# 	mols = self.molecules()
	# 	img_size = (150, 150)
	#
	# 	new_im = Image.new('RGB', (img_size[0]*len(mols), img_size[1]*len(mols)))
	#
	# 	for i, mol1 in enumerate(mols):
	# 		for j, mol2 in enumerate(mols):
	# 			fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2,
	# 							SimilarityMaps.GetMorganFingerprint,
	# 							metric=DataStructs.TanimotoSimilarity,
	# 							size=(50, 50))
	#
	#
	# 			buf = io.BytesIO()
	# 			fig.savefig(buf, format='png', size=img_size)
	# 			buf.seek(0)
	# 			img = Image.open(buf)
	#
	# 			box = (img_size[0]*i, img_size[1]*j)
	# 			new_im.paste(img, box)
	# 			buf.close()
	# 			img.close()
	#
	# 	return new_im



	def cluster_butina(self, cutoff=0.7):
		'''
		Generate a vector with
		'''
		# make a linear input file
		dists = self.distance().values
		data = []
		for i in range(len(self.names())):
			for j in range(i):
				data.append(dists[i, j])

		# cluster them
		cluster_data = Butina.ClusterData(data, len(self.names()), cutoff, isDistData=True)

		# generate a list with cluster belongings
		cluster = [None] * len(self.names())
		for i, clu in enumerate(cluster_data):
			for member in clu:
				cluster[member] = i

		return cluster


	def tsne(self, include_labels=False, color_categories=None, n_components=2, perplexity=10.0, learning_rate=500.0, n_iter=5000):
		"""
		tsne on fingerprint data down to n_components dimensions.
		"""
		data = self.fingerprints_list()
		tsne = TSNE(n_components=n_components,
				perplexity=perplexity,
				learning_rate=learning_rate,
				n_iter=n_iter).fit_transform(data)

		if color_categories is None:
			plot = plt.scatter(tsne[:, 0], tsne[:, 1])
		else:
			plot = plt.scatter(tsne[:, 0], tsne[:, 1], c=color_categories, cmap='rainbow')

		# add text labels if desired
		if include_labels is True:

			# calculate a reasonable offset for the text
			text_spacing_x = abs(max(tsne[:, 0])-min(tsne[:, 0]))/50
			text_spacing_y = abs(max(tsne[:, 1])-min(tsne[:, 1]))/50

			for i, txt in enumerate(self.names()):
				plt.text(tsne[:, 0][i]+text_spacing_x, tsne[:, 1][i]+text_spacing_y, txt, fontsize=8)

		plt.title('t-SNE plot, %s components' % n_components)
		plt.xlabel('PC1')
		plt.ylabel('PC2')

		return plot


	def pca(self, include_labels=False, color_categories=None, n_components=2):
		'''
		PCA using fingerprint bitstrings down to n_components dimensions.
		'''
		data = self.fingerprints_list()
		pca = PCA(n_components=n_components).fit_transform(data)

		if color_categories is None:
			plot = plt.scatter(pca[:, 0], pca[:, 1])
		else:
			plot = plt.scatter(pca[:, 0], pca[:, 1], c=color_categories, cmap='rainbow')

		# add text labels if desired
		if include_labels is True:

			# calculate a reasonable offset for the text
			text_spacing_x = abs(max(pca[:, 0])-min(pca[:, 0]))/50
			text_spacing_y = abs(max(pca[:, 1])-min(pca[:, 1]))/50

			for i, txt in enumerate(self.names()):
				plt.text(pca[:, 0][i]+text_spacing_x, pca[:, 1][i]+text_spacing_y, txt, fontsize=8)

		plt.title('PCA plot, %s components' % n_components)
		plt.xlabel('PC1')
		plt.ylabel('PC2')

		return plot


	def mds(self, include_labels=False, color_categories=None, n_components=2):
		'''
		MDS using already computed distance values down to n_components dimensions.
		'''
		data = self.distance()
		mds = MDS(n_components=n_components, dissimilarity='precomputed', random_state=42).fit_transform(data)

		if color_categories is None:
			plot = plt.scatter(mds[:, 0], mds[:, 1])
		else:
			plot = plt.scatter(mds[:, 0], mds[:, 1], c=color_categories, cmap='rainbow')

		# add text labels if desired
		if include_labels is True:

			# calculate a reasonable offset for the text
			text_spacing_x = abs(max(mds[:, 0])-min(mds[:, 0]))/50
			text_spacing_y = abs(max(mds[:, 1])-min(mds[:, 1]))/50

			for i, txt in enumerate(self.names()):
				plt.text(mds[:, 0][i]+text_spacing_x, mds[:, 1][i]+text_spacing_y, txt, fontsize=8)

		plt.title('MDS plot, %s components' % n_components)
		plt.xlabel('PC1')
		plt.ylabel('PC2')

		return plot





# make a grid with all pairwise mols with substructure
# to mirror the similarity matrix, but visual
