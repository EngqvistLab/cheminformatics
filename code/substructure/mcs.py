from rdkit.Chem import rdFMCS, Draw, AllChem
from rdkit import Chem

"""
Plots a given set of molecules in a grid with their maximum common substructure
highlighted in red. The molecules are given as a list of names and a list with
the corresponding smiles formats. Optionally, the number of rows in the grid 
and the size of each subfigure can be specified.
"""
def plotMCS(names, smiles, molRows=5, subSize=200):
    # make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in smiles]

    # maximum common substructure
    res = rdFMCS.FindMCS(ms)
    mcs = Chem.MolFromSmarts(res.smartsString)

    # align common structure
    AllChem.Compute2DCoords(mcs)
    for m in ms: AllChem.GenerateDepictionMatching2DStructure(m, mcs)

    img = Draw.MolsToGridImage(ms, molsPerRow = molRows, \
            highlightAtomLists = [mol.GetSubstructMatch(mcs) for mol in ms], \
            subImgSize = (subSize, subSize), legends = names)

    img.save('mcs.png')
