from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import itertools
import pandas as pd



smile = ["CCCCCCCCCCCCCCC(C(=O)[O-])O", "[O-]C(=O)C(O)CCCCCCCCCC", "C[C@@H](C(=O)[O-])O", "C(C(=O)[O-])O","CCCCCCCCCCCCCCC(C(=O)[O-])O"]

def tanimoto(smile):

    """
    Input: a list containing smiles.
    Returns: A dataframe containing calculated tanimoto, dice, cosine score between every combination of the fingerprints.
    
    """
    
    smile_fp = {}
    
    for i in smile:
    
        mol=Chem.MolFromSmiles(i)
        smile_fp[i] = AllChem.GetMorganFingerprintAsBitVect(mol,5)
    
    index = pd.MultiIndex.from_tuples(list(itertools.combinations(smile_fp, 2)), names=['Molecule x', 'Molecule Y'])
    similarity_df = pd.DataFrame(index=index, columns=['Tanimoto','Dice','Cosine'])
    
    for x,y in itertools.combinations(smile_fp, 2):
        similarity_df.loc[(x,y), 'Tanimoto'] = DataStructs.FingerprintSimilarity(smile_fp[x],smile_fp[y])
        similarity_df.loc[(x,y), 'Dice']= DataStructs.FingerprintSimilarity(smile_fp[x],smile_fp[y], metric=DataStructs.DiceSimilarity)
        similarity_df.loc[(x,y), 'Cosine']= DataStructs.FingerprintSimilarity(smile_fp[x],smile_fp[y], metric=DataStructs.CosineSimilarity)
    
    similarity_df = similarity_df.sort_values('Tanimoto', ascending=False)
    
    
    return similarity_df