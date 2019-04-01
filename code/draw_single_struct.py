from rdkit import Chem
from rdkit.Chem import Draw

subswithsmiles = {'L-lactate' : "C([C@@H](O)C)(=O)[O-]", 'Glucose' : 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO', 'Glutamate' : "N[C@@H](CCC(=O)[O-])C(=O)[O-]"}

def draw_substrate_structure(SubstratesWSMILES):
    """Takes a dictionary of substrate names and their SMILES and outputs .png files containing drawings of
       their molecular structure."""
    
    for key in SubstratesWSMILES:
        struct = Chem.MolFromSmiles(substrateswithSMILES[key])
        name = key
        Draw.MolToFile(struct, str(name) + ".png",)
        
draw_substrate_structure(subswithsmiles)