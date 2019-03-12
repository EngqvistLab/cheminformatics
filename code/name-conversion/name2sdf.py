"""
Converts an input file with IUPAC-names to an sdf file using the CIRpy-package.
The input and output files are given as the first and second command line
arguments, respectively.
"""

import cirpy, sys

inFile = sys.argv[1]
sdfFile = sys.argv[2]

i = open(inFile,'r')
o = open(sdfFile,'w')

for line in i:
    name = line.rstrip('\n')
    sdf = cirpy.resolve(name,'sdf')
    o.write(sdf)

i.close()
o.close()
