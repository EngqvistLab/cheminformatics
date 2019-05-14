#!/usr/bin/env python3
"""


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

import re


def clean_name(molecule):
    '''
    Remove D-, L-, DL- and similar from the beginning of molecules.
    '''
    # need to make it recursive to capture substrates with multiple patterns
    mol_in = molecule

    molecule = molecule.lower()
    molecule = re.sub('^d-|^l-|^dl-|^\(r\)-|^\(s\)-|^\(\+\)-|^\(\-\)-|^\(\+,\-\)-|^\(\+/-\)-|/in$|/out$', '', molecule)

    # if substrate is unchanged, send it back, otherwise do another round of cleaning
    if mol_in == molecule:
        return molecule
    else:
        return clean_name(molecule)
