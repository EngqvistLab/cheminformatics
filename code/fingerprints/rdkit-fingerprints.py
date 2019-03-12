"""
Converts an input file with smiles-strings to an output file with fingerprints
using RDKit. Three command line arguments specify the input file, the output
file and the type of fingerprint, respectively. The fingerprint formats can be
either 'morgan', 'topological' och 'maccs'.
"""

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys

if len(sys.argv) != 4:
    print("Not enough input arguments")
    print("Specify input file, output file and fingerprint format")
    exit()

inFile = sys.argv[1]    # input file with smiles
outFile = sys.argv[2]   # output file with fingerprints
format = sys.argv[3]    # morgan, topological or maccs

if format not in ['morgan', 'topological', 'maccs']:
    print('Not a valid fingerprint format')
    exit()

# fingerprint function depending on the chosen format
def fingerprint(format):
    if format == 'morgan':
        fp = AllChem.GetMorganFingerprintAsBitVect(m,5)
        # default 2048 number of bits
    elif format == 'topological':
        fp = FingerprintMols.FingerprintMol(m)
    elif format == 'maccs':
        fp = MACCSkeys.GenMACCSKeys(m)
    else:
        print('Not a valid fingerprint format')
        exit()
    return fp

i = open(inFile,'r')
o = open(outFile,'w')

for line in i:
    s = line.rstrip('\n')
    m = Chem.MolFromSmiles(s)
    fp = fingerprint(format) # not a bit string
    o.write(fp.ToBitString() + '\n')

i.close()
o.close()
