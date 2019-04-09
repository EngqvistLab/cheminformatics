import numpy
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

descriptors = {
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3),
    'morgan5':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,5),
    'rdkit':       lambda m: FingerprintMols.FingerprintMol(m)
}

metrics = {
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

"""
Picks a maximally diverse subset of a sert of molecules using the RDKit
MaxMinPicker. The molecules are given as a list of names and a list with the
corresponding smiles formats. Optionally, a list of names with already
chosen molecules can be specified.
"""
def diversityPick(names, smiles, ntopick, descriptor='rdkit', metric='tanimoto', firstpicks=[]):
    # make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in smiles]

    # make a list of fingerprints
    fingerprint = descriptors[descriptor] # fingerprint type
    fps = [fingerprint(x) for x in ms]

    ind = []
    for x in firstpicks:
        ind.append(names.index(x)) # indices of picked molecules

    ds = []
    score = metrics[metric] # similarity score
    for i in range(1,len(fps)):
         ds.extend(score(fps[i],fps[:i],returnDistance=True))

    ids = MaxMinPicker().Pick(numpy.array(ds), len(fps), ntopick, ind)
    return ids
