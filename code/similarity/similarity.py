import numpy as np, pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols

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
Returns a data frame with pairwise similarity scores for a list of molecules,
specified by a list of names and a list of corresponding smiles formats given
as input parameters. The fingerprints and similarity coefficients can be chosen
from the list of descriptors and metrics (default 'rdkit' and 'tanimoto').
"""
def similarity(names, smiles, descriptor='rdkit', metric='tanimoto'):

    if descriptor not in descriptors:
        raise ValueError('Invalid descriptor name ' + descriptor)

    if metric not in metrics:
        raise ValueError('Invalid metric ' + metric)

    # make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in smiles]

    # make a list of fingerprints
    fingerprint = descriptors[descriptor] # fingerprint type
    fps = [fingerprint(x) for x in ms]

    S = []
    score = metrics[metric] # similarity score

    # all against all similarity matrix
    for fp in fps: S.append(score(fp, fps))

    return pd.DataFrame(S,index=names,columns=names)
