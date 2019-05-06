import sys
from rdkit import Chem
from rdkit.Chem import AllChem
#from rdkit.Chem.Fingerprints import FingerprintMols

inFile = sys.argv[1]
outFile = sys.argv[2]

i = open(inFile,'r')
o = open(outFile,'w')

for line in i:
    s = line.rstrip('\n')
    m = Chem.MolFromSmiles(s)

    fp = AllChem.GetMorganFingerprintAsBitVect(m,5)
    #fp = FingerprintMols.FingerprintMol(m)

    o.write(fp.ToBitString() + '\n')

i.close()
o.close()
