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

import json
import pandas as pd
from cheminformatics.helpfunctions import clean_name

### Load environmental variables from the project root directory ###
# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)

NATURAL_MOLS_FILEPATH = join('./data/', 'natural_substrates_filtered.json')
MOLS_FILEPATH = join('./data/', 'substrates_filtered.json')



class _BrendaMolBaseClass(object):
    '''
    A class for the generic methods needed for brenda substrates and products.
    Intended for subclassing
    '''
    def __init__(self, data_filepath, typeof):
        assert typeof in ['substrate', 'product', 'substrates', 'products'], 'Error, the argument "typeof" must hold the values "substrate", "product", "substrates" or "products"'
        self.typeof = typeof
        self.data_filepath = data_filepath
        self._open_data()
        self._flatten_data()


    def _open_data(self):
        '''
        Open the json data file holding all the molecule information.
        '''
        with open(self.data_filepath, 'r') as f:
            self.mol_data = json.loads(f.read())


    def _flatten_data(self):
        '''
        Only keep substrates or products, depending on what is chosen.
        '''
        if self.typeof == 'substrate':
            selection = 'first_substrate'
        elif self.typeof == 'product':
            selection = 'first_product'
        elif self.typeof == 'substrates':
            selection = 'other_substrates'
        elif self.typeof == 'products':
            selection = 'other_products'
        else:
            raise ValueError

        self.selected_data = {}
        for ec in self.mol_data.keys():
            self.selected_data[ec] = sorted(list(set([clean_name(s) for s in self.mol_data[ec][selection]])))


    def ec(self):
        '''
        Get a list of all ec numbers
        '''
        return sorted(list(self.selected_data.keys()))


    def names(self):
        '''
        Obtain a list of all the natural substrates for an ec molecules
        '''
        all_mols = set([])
        for ec in self.selected_data.keys():
            all_mols.update(self.selected_data[ec])
        return sorted(list(all_mols))


    def data_dict(self):
        '''
        Obtain a dictionary with ec number keys and lists of substrates as values
        '''
        return self.selected_data


    def data_frame(self):
        '''
        Obtain a data frame with ec numbers and substrates
        '''
        new_dict = {'ec_number':[], 'molecule':[]}
        for ec in sorted(list(self.selected_data.keys())):
            for mol in self.selected_data[ec]:
                new_dict['ec_number'].append(ec)
                new_dict['molecule'].append(mol)

        return pd.DataFrame(new_dict)



class BrendaNaturalMols(_BrendaMolBaseClass):
    '''
    Class for getting natural molecules from BRENDA data.
    '''
    def __init__(self, typeof='substrate'):
    	_BrendaMolBaseClass.__init__(self, data_filepath=NATURAL_MOLS_FILEPATH, typeof=typeof)



class BrendaMols(_BrendaMolBaseClass):
    '''
    Class for getting molecules from BRENDA data.
    '''
    def __init__(self, typeof='substrate'):
    	_BrendaMolBaseClass.__init__(self, data_filepath=MOLS_FILEPATH, typeof=typeof)
