from rdkit import Chem
from rdkit.Chem import Draw

# Placeholder dictionary containing sampling substrates.
SomeECNum = {'L-lactate' : "C([C@@H](O)C)(=O)[O-]", 'Glucose' : 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO', 'Glutamate' : "N[C@@H](CCC(=O)[O-])C(=O)[O-]"}

# In completed pipeline, suggest to amend this to "draw_substrate_structures(SubstratesWSMILES, filename)"
def draw_substrate_structures(SubstratesWSMILES):   
    """Takes as input a dictionary of substrate names and their SMILES representation, and returns a single .png file containing 
       drawings of all their structures."""
    
    name = []    # Initiate list containing names.
    struct = []  # Initiate list containing mols.

    for key in SubstratesWSMILES:
        name.append(key)
        struct.append(Chem.MolFromSmiles(SubstratesWSMILES[key]))  # Calculate mol from SMILES and append to struct-list.
    
    # Vary row-length according to how many substrates are in dictionary. Modify these values according to preference. 
    if len(SubstratesWSMILES) < 4:
        rowlength = 3
    elif len(SubstratesWSMILES) == 4:
        rowlength = 2
    else:
        rowlength = 4
    
    # Draw grid of molecules and their respective legends. 
    img=Draw.MolsToGridImage(struct, molsPerRow = rowlength, subImgSize=(300, 300), legends = name)
    img.save('SomeECNum.png') # In completed pipeline, suggest to amend this to "img.save(filename + "substrates")"

draw_substrate_structures(SomeECNum)