from rdkit import Chem
from rdkit.Chem import MACCSkeys
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import numpy as np

def dbscan_cluster(smile, eps=0.4, min_samples=2, perplexity=5.0, learning_rate=500.0,n_iter=5000):
    """
    Input: a list containing smiles.
    Output: DBSCAN cluster for the molecules, visualised with colours in a 2D plot and 
    a list with the smiles of each cluster.
    Data processing: The data have first been processed first with PCA down to 50 dimensions 
    (or less if the number of smiles is less) and then with tSNE, down to 2 dimensions. 
    The result after tSNE will be plotted. To achive a nice looking plot, 
    the parameters perplexity and learning rate can be changed to fit the dataset.
    Input to tSNE:
    perplexity: A approximation on how many close neighbours for each point.
    n_iter: maximum numbers of iterations
    learning_rate: Usually in the rang [10:1000]
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

    #PCA on fingerprint data down to 50 dimensions or if that number is less than 50, to the number of fps.
    data = np.array(bit_list)
    if len(smile_list) >= 50:
        n_components = 50
    else:
        n_components = len(bit_list)  
    pca = PCA(n_components=n_components).fit_transform(data)

    #tSNE on the resulting PCA data
    data_embedded = TSNE(n_components=2, perplexity=perplexity, learning_rate=learning_rate,n_iter=n_iter).fit_transform(pca)
    
    #DBSCAN, another clustering algoritm on the same data.
    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='jaccard').fit(np.array(bit_list))
    print(f"The clustering labels are: {clustering.labels_}. \n Each color in the plot correspond to a specific cluster.")
    dbscan_plot = plt.scatter(data_embedded[:, 0], data_embedded[:, 1], c=clustering.labels_, cmap='rainbow')

    #Generates a list that contains lists with all smiles that are in a cluster.
    smile_cluster = []
    for number in set(clustering.labels_):
        index = 0
        lista = []
        for label in clustering.labels_:
            if label == number:
                lista.append(smile_list[index])
            index += 1
        smile_cluster.append(lista)     
    
    return dbscan_plot, smile_cluster