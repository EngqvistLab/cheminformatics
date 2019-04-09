from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys
from rdkit.ML.Cluster import Butina
import numpy as np

def butina_clustering(smile, cutoff=0.8):
    
    """
    Input: a list containing smiles.
    Output: Butina cluster for the molecules, visualised with colours in a 2D plot and 
    a list with the smiles of each cluster.
    Butina clustering: Butina clustering algorithm from RDKit is used on the fingerprint vectors.
    The Tanimoto score is used for the clustering. The function returns the list: smile_cluster, 
    that contains lists with all the smiles in each of the clusters. To achive a good clustering, 
    parameter cutoff can be changed to fit the dataset.
    Butina parameters:
    cutoff: The cluster centroid is at least this similar to all other cluster members. 
    """

    # Calculates fingerprints for all smiles.
    smiles_fps = {}
    for molecule in smile:
        try:
            mol = Chem.MolFromSmiles(molecule)
            smiles_fps[molecule] = MACCSkeys.GenMACCSKeys(mol)
        except: 
            print(f"The fingerprint for {molecule} could not be generated.")
            
    # Makes 2 lists, one with all fingerprints and one with all smiles.        
    fps_ = []
    for value in smiles_fps.values():
        fps_.append(value)
        
    # Creates one vector with the bits in the fingerprint and one list with the corresponding smiles.
    bit_list = [] 
    smile_list = []
    for key, value in smiles_fps.items():
        bits = []
        bit_fp = value.ToBitString()
        for bit in bit_fp:
            bits.append(int(bit)) 
        bit_list.append(bits)
        smile_list.append(key)

    #Butina cluster: Generate the distance matrix (metric=tanimoto):
    dists = []
    len_fps = len(smile_list) 
    for fp in range(1, len_fps):
        similarities = DataStructs.BulkTanimotoSimilarity(fps_[fp], fps_[: fp])
        dists.extend([1-x for x in similarities])
        
    #Cluster the data:
    cluster = Butina.ClusterData(dists, len_fps, cutoff, isDistData=True)
    print(f"Clusters: {cluster} \nEach tuple represent a cluster, the numbers represent the index in smile_list.")
       
    #Generates smile_cluster: a list with lists of all the smiles in one cluster and the cluster labels.
    counter = 0
    smile_cluster = []
    cluster_label = ['#']*len_fps
    for molecule in cluster:
        clusters = []
        for numb in molecule:
            cluster_label[numb] = counter
            clusters.append(smile_list[numb])
        counter += 1
        smile_cluster.append(clusters)

    return smile_cluster  
    