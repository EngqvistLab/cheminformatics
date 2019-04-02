from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys

descriptors = {
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3),
    'morgan5':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,5),
    'rdkit':       lambda m: FingerprintMols.FingerprintMol(m)
}

"""
Function that calculates fingerprints for a list of smiles strings. The
fingerprint formats can be chosen from the list of descriptors.
"""
def fingerprint(smiles, fingerprint):
    if fingerprint not in descriptors:
        raise ValueError('Invalid descriptor name ' + descriptor)

    # make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in smiles]

    # make a list of fingerprints
    fingerprint = descriptors[descriptor] # fingerprint type
    fps = [fingerprint(x) for x in ms]
    return fps
