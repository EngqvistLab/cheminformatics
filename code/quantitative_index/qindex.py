import numpy as np
from rdkit import DataStructs

metrics = {
    'asymmetric':    DataStructs.AsymmetricSimilarity,
    'braunblanquet': DataStructs.BulkBraunBlanquetSimilarity,
    'cosine':        DataStructs.BulkCosineSimilarity,
    'dice':          DataStructs.BulkDiceSimilarity,
    'kulczynski':    DataStructs.BulkKulczynskiSimilarity,
    'mcconnaughey':  DataStructs.BulkMcConnaugheySimilarity,
    'rogotgoldberg': DataStructs.BulkRogotGoldbergSimilarity,
    'russel':        DataStructs.BulkRusselSimilarity,
    'sokal':         DataStructs.BulkSokalSimilarity,
    'tanimoto':      DataStructs.BulkTanimotoSimilarity
}

def mean_dist_set(fps):
    """
    Returns an overall set dissimilarity for a list of fingerprints from RDKit, as defined in
    "Nath A, Atkins WM. A quantitative index of substrate promiscuity. Biochemistry. 2008."
    The fingerprints have to be the same length.
    """

    # convert to numpy array and add vectors
    s = np.array(list(fps[0]))
    for i in range(1,len(fps)):
        s = s + np.array(list(fps[i]))

    l = (s == len(fps)).sum() # bits common to all fps
    k = np.count_nonzero(s) - l # on bits for at least one but not all fps

    return k/(k+l)

def mean_dist(fps,metric):
    """
    Returns the mean distances between molecules represented by a list of
    fingerprints from RDKit, with a specified metric.
    """
    score = metrics[metric] # similarity score

    n = len(fps)
    md = []

    for fp in fps:
        s = sum(score(fp,fps,returnDistance=1))-1
        md.append(s/(n-1))

    return md

def qindex(fps,e,metric):
    """
    Returns a quantitative index of substrate promiscuity (J) as defined in
    "Nath A, Atkins WM. A quantitative index of substrate promiscuity. Biochemistry. 2008."
    for a list of fingerprints from RDKit, with corresponding values for e = kcat/KM.
    Metric specifies the chosen similarity score.
    """
    N = len(fps)

    md = np.array(mean_dist(fps,metric)) # mean distances
    nmd = md/mean_dist_set(fps) # normalized mean distances
    e = np.array(e) # kcat/KM
    p = e/sum(e)
    A = -N/(sum(nmd)*np.log(N))

    return A * sum(nmd * p * np.log(p))
