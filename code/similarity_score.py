from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import itertools
import pandas as pd
from rdkit.Chem import MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs



def similarity(smile, fingerprint_algorithm = "morgan"):

    """
    Input: a list containing smiles and fingerprint algorithm: morgan (default), MACCS or topological.
    Returns: A dataframe containing calculated tanimoto, dice, cosine score between every combination of the fingerprints.
    
    """
    
    smile_fp = {}
    
    # Calculates the choosen fingerprint.
    
    for i in smile:
    
        mol=Chem.MolFromSmiles(i)
        
        if fingerprint_algorithm == "morgan" :
            smile_fp[i] = AllChem.GetMorganFingerprintAsBitVect(mol,5)
        elif fingerprint_algorithm == "maccs": 
            smile_fp[i] = MACCSkeys.GenMACCSKeys(mol)
        elif fingerprint_algorithm == "topological":
            smile_fp[i] = FingerprintMols.FingerprintMol(mol)
        else:
            raise NotImplementedError("Not implemented. Available fingerprint algorithms: morgan, maccs, topological.")
            
    # Calculates tanimoto, dice and cosine score for each combination of fingerprints and makes a dataframe.
    
    index = pd.MultiIndex.from_tuples(list(itertools.combinations(smile_fp, 2)), names=['Molecule x', 'Molecule Y'])
    similarity_df = pd.DataFrame(index=index, columns=['Tanimoto','Dice','Cosine'])
    
    for x,y in itertools.combinations(smile_fp, 2):
        similarity_df.loc[(x,y), 'Tanimoto'] = DataStructs.FingerprintSimilarity(smile_fp[x],smile_fp[y])
        similarity_df.loc[(x,y), 'Dice']= DataStructs.FingerprintSimilarity(smile_fp[x],smile_fp[y], metric=DataStructs.DiceSimilarity)
        similarity_df.loc[(x,y), 'Cosine']= DataStructs.FingerprintSimilarity(smile_fp[x],smile_fp[y], metric=DataStructs.CosineSimilarity)
    
    similarity_df = similarity_df.sort_values('Tanimoto', ascending=False)
    
    
    return similarity_df
