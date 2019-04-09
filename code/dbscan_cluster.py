from rdkit import Chem
from rdkit.Chem import MACCSkeys
from sklearn.cluster import DBSCAN
import numpy as np

def dbscan_clustering(smile, eps=0.3, min_samples=2):
    """
    Input: a list containing smiles.
    Output: DBSCAN cluster for the molecules, and a list with the smiles of each cluster.
    DBSCAN clustering: DBSCAN clustering algorithm from scikit learn is used on the 
    fingerprint vectors, with Jaccard as the pairwise distance metric. 
    The different clusters are visualised with different colours in the plot. 
    The function returns a list,called smile_cluster, containing lists with all the smiles in each cluster.
    To achive a good clustering parameters eps and min_sample can be changed to fit the dataset.
    Input to DBSCAN: 
    eps: The maximum distance between two points for them to be considered to belong in the same neighborhood.
    min_sample: The number of samples in a neighborhood for a point to be considered as a core point (including the point itself).
    """
    
    # Calculates fingerprints for all smiles.
    smiles_fps = {}
    for molecule in smile:
        try:
            mol = Chem.MolFromSmiles(molecule)
            smiles_fps[molecule] = MACCSkeys.GenMACCSKeys(mol)
        except: 
            print(f"The fingerprint for {molecule} could not be generated.")
        
    # Creates one vector (list) with all the bits in a fingerprint and one list with all the corresponding smiles.
    bit_list = [] 
    smile_list = []
    for key, value in smiles_fps.items():
        bits = []
        bit_fp = value.ToBitString()
        for bit in bit_fp:
            bits.append(int(bit))  
        bit_list.append(bits)
        smile_list.append(key)

    #DBSCAN, another clustering algoritm on the same data.
    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='jaccard').fit(np.array(bit_list))
    print(f"The clustering labels are: {clustering.labels_}. \nThe index of this list corresponds to the index of smile_list.")
    
    #Generates a list that contains lists with all smiles that are in a cluster.
    dbscan_cluster = []
    for number in set(clustering.labels_):
        index = 0
        lista = []
        for label in clustering.labels_:
            if label == number:
                lista.append(smile_list[index])
            index += 1
        dbscan_cluster.append(lista)     
    
    return dbscan_cluster
    