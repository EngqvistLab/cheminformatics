from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import itertools

smile = ["CCCCCCCCCCCCCCC(C(=O)[O-])O", "[O-]C(=O)C(O)CCCCCCCCCC", "C[C@@H](C(=O)[O-])O", "C(C(=O)[O-])O","CCCCCCCCCCCCCCC(C(=O)[O-])O"]

def tanimoto(smile):

    """
    Input: a list containing smiles.
    Returns: A list containing calculated tanimoto score between every combination of the fingerprint.
    
    """
    fp = []
    for i in smile:
    
        mol=Chem.MolFromSmiles(i)
        fp.append(AllChem.GetMorganFingerprintAsBitVect(mol,5))
    
    for x,y in itertools.combinations(fp, 2):
        
        tanimoto_score = []
        tanimoto_score.append(DataStructs.FingerprintSimilarity(fp[fp.index(x)],fp[fp.index(y)]))
    
    return tanimoto_score