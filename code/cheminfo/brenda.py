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
import re

### Load environmental variables from the project root directory ###
# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)


class BrendaSubstrates():
    '''
    Class for getting, parsing, and obtaining subsets of BRENDA data.
    '''
    def __init__(self):
        pass

    def remove_dl(substrate):
        '''
        Remove D-, L-, DL- and similar from the beginning of substrates
        '''
        substrate = substrate.lower()
        return re.sub('^d-|^l-|^dl-|^\(r\)-|^\(s\)-|^\(\+\)-|^\(\-\)-', '', substrate)


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
###################
