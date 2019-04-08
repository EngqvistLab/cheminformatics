from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys
from rdkit.ML.Cluster import Butina
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np

def butina_cluster(smile, cutoff=0.8, perplexity=3.0, learning_rate=400.0,n_iter=5000):
    
    """
    Input: a list containing smiles.
    Output: Butina cluster for the molecules, visualised with colours in a 2D plot and 
    a list with the smiles of each cluster.
    Data processing: The data have first been processed first with PCA down to 50 dimensions 
    (or less if the number of smiles is less) and then with tSNE, down to 2 dimensions. 
    The result after tSNE will be plotted. To achive a nice looking plot, the parameters 
    perplexity and learning rate can be changed to fit the dataset.
    Input to tSNE:
    perplexity: A approximation on how many close neighbours for each point.
    n_iter: maximum numbers of iterations
    learning_rate: Usually in the rang [10:1000] 
    Butina clustering: Butina clustering algorithm from RDKit is used on the fingerprint vectors.
    The Tanimoto score is used for the clustering. The different clusters are visualised with 
    different colours in the plot. The function returns the list: smile_cluster, that contains lists
    with all the smiles in each of the clusters. To achive a good clustering parameter cutoff can be 
    changed to fit the dataset.
    Input to Butina:
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

    #PCA on fingerprint data down to 50 dimensions or if that number is less than 50, to the number of fps.
    data = np.array(bit_list)
    if len(smile_list) >= 50:
        n_components = 50
    else:
        n_components = len(bit_list)  
    pca = PCA(n_components=n_components).fit_transform(data)

    #tSNE on the resulting PCA data
    data_embedded = TSNE(n_components=2, perplexity=perplexity, learning_rate=learning_rate,n_iter=n_iter).fit_transform(pca)
    
    #Butina cluster: Generate the distance matrix (metric=tanimoto):
    dists = []
    len_fps = len(smile_list) 
    for fp in range(1, len_fps):
        similarities = DataStructs.BulkTanimotoSimilarity(fps_[fp], fps_[: fp])
        dists.extend([1-x for x in similarities])

    #Cluster the data:
    cluster = Butina.ClusterData(dists, len_fps, cutoff, isDistData=True)
    print(f"Clusters: {cluster} \nEach tuple represent a cluster, visualised with a specific colour.")
    
    
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
    #Generates a plot of the molecules, coloured according to cluster. 
    butina_plot = plt.scatter(data_embedded[:, 0], data_embedded[:, 1], c=cluster_label, cmap='rainbow')

    return butina_plot, smile_cluster
    