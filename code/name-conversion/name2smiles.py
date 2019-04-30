"""
Converts an input file with IUPAC-names to an output file with SMILES-strings
using the CIRpy-package. The input and output files are given as the first and
second command line arguments, respectively.
"""

import cirpy, sys

inFile = sys.argv[1]
outFile = sys.argv[2]

i = open(inFile,'r')
o = open(outFile,'w')

for line in i:
    name = line.rstrip('\n')
    smiles = cirpy.resolve(name,'smiles')
    o.write(smiles + '\n')

i.close()
o.close()
