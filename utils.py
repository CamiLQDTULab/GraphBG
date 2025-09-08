# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause
## Preprocessing is very important
# sc.pp.highly_variable_genes(adata_bk, flavor="seurat_v3", n_top_genes=3000)
# sc.pp.normalize_total(adata_bk, target_sum=1e4)
# sc.pp.log1p(adata_bk)
# sc.pp.scale(adata_bk, zero_center=False, max_value=10)

import time
import warnings
from itertools import cycle, islice
import matplotlib.pyplot as plt
import numpy as np
from sklearn import cluster, datasets, mixture
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize, StandardScaler 
import pandas as pd
from sklearn import metrics
import scanpy as sc
import ot
from sklearn.decomposition import PCA
# from GraphST.utils import mclust_R
# from GraphST.utils import refine_label
import scipy as sp
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import mode
from sklearn.cluster import (KMeans, MiniBatchKMeans, AgglomerativeClustering, SpectralClustering, Birch)
from sklearn.metrics.cluster import adjusted_rand_score, normalized_mutual_info_score, homogeneity_score, completeness_score

####################################### Leiden clustering ###############################################################
import scanpy as sc
from anndata import AnnData
from sklearn import preprocessing
def eval_leiden(adata, resol):
    sc.tl.leiden(adata, key_added='louvain', resolution = resol)
    louv_labels = np.array(adata.obs['louvain'].tolist())

    le = preprocessing.LabelEncoder().fit(louv_labels)
    cell_labels_pred  = le.transform(louv_labels)
    return int(np.unique(cell_labels_pred).shape[0])

def eval_bisect_leiden(adata, minresolution, maxresolution, n_clusters):
    M = -1
    k1 = n_clusters + 1
    k2 = 0
    iterations = 0 
    #adata = AnnData(X = data)
    #sc.pp.neighbors(adata, n_neighbors=20, use_rep='X')

    # find minresolution and maxresolution s.t k1 < n_clusters < k2 
    while k1 > n_clusters and iterations < 8:
        minresolution = minresolution/2 
        k1 = eval_leiden(adata, minresolution)
        if k1 == n_clusters:
            return minresolution
        iterations = iterations + 1
    while k2 < n_clusters and iterations < 8:
        maxresolution = 2*maxresolution
        k2 = eval_leiden(adata, maxresolution)
        if k2 == n_clusters:
            return maxresolution
        iterations = iterations + 1

    # bisection main
    while iterations < 40 and abs(maxresolution - minresolution) > 1e-5:
        M = (minresolution + maxresolution)/2 
        iterations = iterations + 1
        k3 = eval_leiden(adata, M)
        if k3 == n_clusters:
            return M 
        elif k3 < n_clusters:
            minresolution = M 
        else:
            maxresolution = M

    if iterations >= 40:
        print("bisection algorithm could not find the right resolution")

def leiden_exact_K(X_dimred, n_clusters):
    adata = AnnData(X=X_dimred)
    sc.pp.neighbors(adata,n_neighbors=10, use_rep='X')

    resol = eval_bisect_leiden(adata, 0.5, 1.7, n_clusters)

    sc.tl.leiden(adata, key_added='louvain', resolution = resol)
    louv_labels = np.array(adata.obs['louvain'].tolist())
    le = preprocessing.LabelEncoder().fit(louv_labels)
    cell_labels_pred  = le.transform(louv_labels)
    return cell_labels_pred

############################################## Refinement ################################################
############################################## Refinement ################################################
from scipy.spatial import cKDTree
from scipy.stats import mode
def cKD_refine_label(coords, labels, k):
    # Step 1: Build KD-Tree
    tree = cKDTree(coords.copy())

    # Step 2: Find k-nearest neighbors for each spot
    # k+1 because the closest point is itself
    distances, neighbors = tree.query(coords, k=k+1)

    # Exclude self-neighbor (first column)
    neighbors = neighbors[:, 1:]

    # Step 3: Reassign labels
    new_labels = labels.copy()
    for i, nbrs in enumerate(neighbors):
        # Get the labels of neighboring spots
        neighbor_labels = labels[nbrs]
        # Find the most common label among neighbors
        # most_common_label = mode(neighbor_labels).mode[0]
        most_common_label = np.atleast_1d(mode(neighbor_labels).mode)[0]
        # Reassign the label
        new_labels[i] = most_common_label
    return (new_labels)

################################################## TSNE ##################################################
if False:
    ## TSNE plot
    from sklearn.manifold import TSNE
    import matplotlib.pyplot as plt
    import seaborn as sns
    tsne = TSNE(n_components=2, perplexity=30, random_state=42)
    # Transform the data
    X_embedded = tsne.fit_transform(X)
    # Convert to DataFrame for visualization
    df_tsne = pd.DataFrame(X_embedded, columns=['TSNE1', 'TSNE2'])
    plt.figure(figsize=(10, 7))
    sns.scatterplot(x='TSNE1', y='TSNE2', palette='tab10', data=df_tsne, alpha=0.7)
    plt.title("t-SNE Visualization of Digits Dataset")
    plt.legend(title="Digit Label", bbox_to_anchor=(1, 1))
    plt.show()

    
def mclust_R_raw(X_data, num_cluster, modelNames='EEE', random_seed=2020):
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    
    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")

    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']
    
    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(X_data), num_cluster, modelNames)
    mclust_res = np.array(res[-2])

    mclust_res = mclust_res.astype('int')
    # adata.obs['mclust'] = adata.obs['mclust'].astype('category')
    return mclust_res


def clustering(adata, n_clusters=7, radius=50, key='emb', refinement=False):
    pca = PCA(n_components=20, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    adata.obsm['emb_pca'] = embedding 
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 


############################################## Clustering (unimodal) ###########################################################
############################################## Clustering (unimodal) ###########################################################
# SOTA method
from CustomVBG import CustomBayesianGaussianMixture
def runGraphBG(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, covariance_type='full', n_init = 5, max_iter = 1000, linkage='ward'):
    print(covariance_type, linkage)
    print("Data input shape n: ", adata.obsm[key].shape)
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    connectivity = 0.5 * (connectivity + connectivity.T)
    if 'pca' not in key:
        print("Run PCA")
        pca = PCA(n_components=n_components, random_state=42) 
        embedding = pca.fit_transform(np.asarray(adata.obsm[key].copy()))
        adata.obsm['raw_pca'] = embedding
        embedding = connectivity.dot(embedding)
        # A = connectivity.dot(embedding)
        # B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
        # embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    else:
        embedding = connectivity.dot(adata.obsm[key])

    gmm = CustomBayesianGaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42, init_params = 'random_from_data', n_init = n_init, max_iter = max_iter) # 5, 1000
    adata.obs['domain']  = gmm.fit_predict(embedding)
    
    if refinement:  
        # new_type = refine_label(adata, radius, key='domain')
        # adata.obs['domain'] = new_type 
        adata.obs['domain'] = cKD_refine_label(np.array(adata.obsm['spatial']), adata.obs['domain'], k = radius)


from sklearn.cluster import AgglomerativeClustering   
from sklearn.mixture import GaussianMixture
def clustering_v2(adata, n_clusters=7, radius=50, key='emb', refinement=False):
    pca = PCA(n_components=20, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    # adata.obsm['emb_pca'] = embedding     
    gmm = GaussianMixture(n_components=n_clusters, covariance_type='full', random_state=42, n_init = 5, max_iter = 50)
    gmm.fit(embedding)
    
    adata.obs['domain'] = gmm.predict(embedding)
    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
        
from numpy import dot, array
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
def clustering_v3(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20):
    print("Data input shape: ", adata.obsm[key].shape)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    
    adata.obsm['raw_pca'] = embedding

    # sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # sc.tl.pca(adata, use_highly_variable=True, n_comps = 10)
    # embedding = adata.obsm['X_pca']
    
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # print(A.shape, B.shape)
    embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    
    adata.obsm['emb_pca'] = np.asarray(embedding) 
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
        
def clustering_v3b(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20):
    print("Data input shape: ", adata.obsm[key].shape)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    
    adata.obsm['raw_pca'] = embedding

    # sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # sc.tl.pca(adata, use_highly_variable=True, n_comps = 10)
    # embedding = adata.obsm['X_pca']
    
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # print(A.shape, B.shape)
    embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    
    adata.obsm['emb_pca'] = np.array(embedding) 
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 

def clustering_v3c(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20):
    print("Data input shape: ", adata.obsm[key].shape)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    
    adata.obsm['raw_pca'] = embedding
    
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
  
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    # B = np.diag(np.sqrt(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1)))))
    # embedding = B.dot(embedding)
    # A = connectivity.dot(embedding)
    # embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # embedding = B.dot(embedding)
    # embedding = connectivity.dot(embedding)
    embedding = dot(B, A)# + 0.01*adata.obsm['raw_pca']
    
    adata.obsm['emb_pca'] = np.array(embedding) 
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
        
def clustering_v3d(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20):
    print("Data input shape: ", adata.obsm[key].shape)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    
    adata.obsm['raw_pca'] = embedding

    # sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # sc.tl.pca(adata, use_highly_variable=True, n_comps = 10)
    # embedding = adata.obsm['X_pca']
    
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # print(A.shape, B.shape)
    embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    
    adata.obsm['emb_pca'] = np.asarray(embedding) 
    
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
           
def clustering_v3e(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, modelName='EEE'):
    # Concat features
    print("Data input shape: ", adata.obsm[key].shape, 'modelName: ', modelName)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())  
    adata.obsm['raw_pca'] = np.asarray(embedding)  
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # print(A.shape, B.shape)
    embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    
    adata.obsm['emb_pca'] = np.asarray(embedding) 
    # adata.obsm['emb_pca'] = np.concatenate((adata.obsm['emb_pca'], adata.obsm['raw_pca']), axis=1)
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters, modelNames=modelName)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
        
# clustering_v4(adata, n_clusters=n_clusters, radius=50, key='raw', refinement=True, inner_neighbors = 2, outer_neighbors = 4)
def clustering_v4(adata, n_clusters=7, radius=50, key='emb', refinement=True, inner_neighbors = 3, outer_neighbors = 6, n_components = 20, eps = 0.1):
    print("Data input shape: ", adata.obsm[key].shape)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    
    adata.obsm['raw_pca'] = embedding

    # sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # sc.tl.pca(adata, use_highly_variable=True, n_comps = 10)
    # embedding = adata.obsm['X_pca']
    
    connectivity1 = kneighbors_graph(adata.obsm['spatial'], n_neighbors=inner_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity1 = 0.5 * (connectivity1 + connectivity1.T)
    
    connectivity2 = kneighbors_graph(adata.obsm['spatial'], n_neighbors=outer_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity2 = 0.5 * (connectivity2 + connectivity2.T)
    connectivity = connectivity2 - connectivity1
    
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # print(A.shape, B.shape)
    embedding = np.asarray(dot(B, A)) #+ 0.5*embedding #dot(connectivity, embedding)
    
    # Ain = connectivity1.dot(adata.obsm['raw_pca'])
    # connectivity1 = connectivity2 - connectivity1
    
    Ain = connectivity1.dot(embedding)
    Bin = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity1.todense(), axis=1))))
    # print(A.shape, B.shape)
    embedding_in = np.asarray(dot(Bin, Ain)) #+ 0.5*embedding #dot(connectivity, embedding)
    
    # adata.obsm['emb_pca'] = eps*embedding/np.linalg.norm(embedding) + embedding_in/np.linalg.norm(embedding_in)
    # adata.obsm['emb_pca'] = np.concatenate((np.asarray(embedding_in), eps*np.asarray(embedding) ), axis=1)
    adata.obsm['emb_pca'] = embedding_in
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 

def clustering_v5(adata, n_clusters=7, radius=50, key='emb', refinement=True, inner_neighbors = 3, outer_neighbors = 6, n_components = 20, eps = 0.1):
    print("Data input shape: ", adata.obsm[key].shape)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    
    adata.obsm['raw_pca'] = embedding

    # sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # sc.tl.pca(adata, use_highly_variable=True, n_comps = 10)
    # embedding = adata.obsm['X_pca']
    
    connectivity1 = kneighbors_graph(adata.obsm['spatial'], n_neighbors=inner_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity1 = 0.5 * (connectivity1 + connectivity1.T)
    
    connectivity2 = kneighbors_graph(adata.obsm['spatial'], n_neighbors=outer_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity2 = 0.5 * (connectivity2 + connectivity2.T)
    connectivity = connectivity1
    
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # print(A.shape, B.shape)
    embedding = np.asarray(dot(B, A)) #+ 0.5*embedding #dot(connectivity, embedding)
    
    # Ain = connectivity1.dot(adata.obsm['raw_pca'])
    Ain = connectivity1.dot(embedding)
    Bin = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity1.todense(), axis=1))))
    # print(A.shape, B.shape)
    embedding_in = np.asarray(dot(Bin, Ain)) #+ 0.5*embedding #dot(connectivity, embedding)
    
    # adata.obsm['emb_pca'] = eps*embedding/np.linalg.norm(embedding) + embedding_in/np.linalg.norm(embedding_in)
    # adata.obsm['emb_pca'] = np.concatenate((np.asarray(embedding_in), eps*np.asarray(embedding) ), axis=1)
    adata.obsm['emb_pca'] = embedding_in
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
        
def clustering_v3a(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, covariance_type='full', linkage='ward'):
    print(covariance_type, linkage)
    print("Data input shape: ", adata.obsm[key].shape)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(np.asarray(adata.obsm[key].copy()))
    adata.obsm['raw_pca'] = embedding
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    
    hierarchical = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
    initial_labels = hierarchical.fit_predict(embedding)
    
    gmm = GaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42, init_params = 'random_from_data', n_init = 5, max_iter = 1000)
    gmm.means_init = np.array([embedding[initial_labels == i].mean(axis=0) for i in range(n_clusters)])
    # print("dddd")
    # Step 4: Fit GMM
    gmm.fit(embedding)
    cluster_labels = gmm.predict(embedding)
    adata.obs['domain'] = cluster_labels
    
    print('BIC = ', gmm.bic(embedding), 'AIC = ', gmm.aic(embedding))
    
    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
 

import numpy as np
from sklearn.neighbors import NearestNeighbors
def compute_connectivity_matrix(data1, data2, k1, k2):
    """
    Compute the connectivity matrix based on two-step nearest neighbor selection.
    
    Parameters:
        data1 (ndarray): Data for first modality (N x D1).
        data2 (ndarray): Data for second modality (N x D2).
        k1 (int): Number of nearest neighbors in first modality.
        k2 (int): Number of nearest neighbors in second modality.
    
    Returns:
        connectivity_matrix (ndarray): Binary matrix (N x N), where M[i, j] = 1 if j is a selected neighbor of i.
    """
    N = data1.shape[0]
    connectivity_matrix = np.zeros((N, N), dtype=int)
    
    # Step 1: Find k1 nearest neighbors for all points using modality 1
    nbrs1 = NearestNeighbors(n_neighbors=k1 + 1, algorithm='auto').fit(data1)
    _, indices1 = nbrs1.kneighbors(data1)

    for i in range(N):
        k1_neighbors = indices1[i][1:]  # Exclude self (first element)

        # Step 2: Find k2 nearest neighbors among k1_neighbors using modality 2
        nbrs2 = NearestNeighbors(n_neighbors=min(k2, len(k1_neighbors)), algorithm='auto').fit(data2[k1_neighbors])
        _, indices2 = nbrs2.kneighbors(data2[i].reshape(1, -1))
        
        # Map indices back to original dataset
        final_neighbors = [k1_neighbors[j] for j in indices2[0]]
        
        # Update connectivity matrix
        connectivity_matrix[i, final_neighbors] = 1

    return connectivity_matrix

def clustering_bigraph(adata, n_clusters=7, radius=50, key='emb', refinement=False, k1 = 15, k2 = 5, n_components = 20, covariance_type='full', n_init = 5, max_iter = 1000, linkage='ward'):
    print(covariance_type, linkage)
    print("Data input shape n: ", adata.obsm[key].shape)
    # connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # connectivity = 0.5 * (connectivity + connectivity.T)
    if 'pca' not in key:
        print("Run PCA")
        pca = PCA(n_components=n_components, random_state=42) 
        embedding = pca.fit_transform(np.asarray(adata.obsm[key].copy()))
        adata.obsm['raw_pca'] = embedding
        connectivity = compute_connectivity_matrix(adata.obsm['spatial'], embedding, k1, k2)
        connectivity = 0.5 * (connectivity + connectivity.T)
        embedding = connectivity.dot(embedding)
        # A = connectivity.dot(embedding)
        # B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
        # embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    else:
        connectivity = compute_connectivity_matrix(adata.obsm['spatial'], adata.obsm[key], k1, k2)
        # connectivity = compute_connectivity_matrix(adata.obsm[key], adata.obsm['spatial'],  k1, k2)
        connectivity = 0.5 * (connectivity + connectivity.T)
        embedding = connectivity.dot(adata.obsm[key])
    
    # hierarchical = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
    # initial_labels = hierarchical.fit_predict(embedding)
    
    #gmm = GaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42)
    gmm = BayesianGaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42, init_params = 'random_from_data', n_init = n_init, max_iter = max_iter) # 5, 1000
    # gmm.means_init = np.array([embedding[initial_labels == i].mean(axis=0) for i in range(n_clusters)])
    # gmm.fit(embedding)
    # cluster_labels = gmm.predict(embedding)
    # adata.obs['domain'] = cluster_labels
    adata.obs['domain']  = gmm.fit_predict(embedding)
    
    # Calculate ELBO Score
    # print("ELBO Score = ", gmm.lower_bound_)
    if refinement:  
        # new_type = refine_label(adata, radius, key='domain')
        # adata.obs['domain'] = new_type 
        adata.obs['domain'] = cKD_refine_label(np.array(adata.obsm['spatial']), adata.obs['domain'], k = radius)

def clustering_metacells(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, covariance_type='full', n_init = 5, max_iter = 1000, n_metacells = 100):
    print("Data input shape n: ", adata.obsm[key].shape, covariance_type)
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    connectivity = 0.5 * (connectivity + connectivity.T)
    if 'pca' not in key:
        # print("Run PCA")
        pca = PCA(n_components=n_components, random_state=42) 
        embedding = pca.fit_transform(np.asarray(adata.obsm[key].copy()))
        adata.obsm['raw_pca'] = embedding
        embedding = connectivity.dot(embedding)
        # A = connectivity.dot(embedding)
        # B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
        # embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    else:
        embedding = connectivity.dot(adata.obsm[key])
    
    # Compute meta-cells
    clusteringobj = MiniBatchKMeans(n_clusters=n_metacells, random_state=42, tol=1e-3, max_iter=10000)
    cluster_ = clusteringobj.fit_predict(embedding)
    # cluster_ = clusteringobj.fit_predict(adata.obsm['spatial'])
    # cluster_ = clusteringobj.fit_predict(adata.obsm['raw_pca'])
    cluster_vec = np.unique(cluster_)
    clusterIDvec = []
    clusterValuevec = []
    clusterSpatialvec = []
    cluster2Cellsdict = {}
    for clusterID in cluster_vec:
        sliceClusterID = 'ClusterID' + str(clusterID)
        sliceClusterValue = np.mean(embedding[cluster_ == clusterID], axis=0)
        clusterIDvec.append(sliceClusterID)
        clusterValuevec.append(sliceClusterValue)
        clusterSpatialvec.append(np.mean(np.array(adata.obsm['spatial'])[cluster_ == clusterID], axis=0))
        cluster2Cellsdict[sliceClusterID] = adata.obs_names[cluster_ == clusterID]
        
    meta_emb = np.array(clusterValuevec)
    meta_spatial = np.array(clusterSpatialvec)
    
    # connectivity = kneighbors_graph(meta_spatial, n_neighbors=3, include_self=False)
    # connectivity = 0.5 * (connectivity + connectivity.T)
    # meta_emb = connectivity.dot(meta_emb)
    # print(meta_emb.shape)
    
    gmm = BayesianGaussianMixture(n_components=n_clusters, covariance_type='full', random_state=42, init_params = 'random_from_data', n_init = n_init, max_iter = max_iter) # 5, 1000
    meta_labels  = gmm.fit_predict(meta_emb)
    # if refinement:  
    #     # print(meta_spatial.shape, meta_labels.shape)
    #     meta_labels  = cKD_refine_label(meta_spatial, meta_labels, k = radius)
    # Extrapolation
    adata.obs['domain'] = 0
    for i in range(len(clusterIDvec)):
        metacell = clusterIDvec[i]
        adata.obs.loc[cluster2Cellsdict[metacell], 'domain'] = meta_labels[i]
    if refinement:  
        adata.obs['domain-1'] = cKD_refine_label(np.array(adata.obsm['spatial']), adata.obs['domain'], k = radius)
    true_labels = pd.factorize(adata.obs['Region'])[0]
    score_1 = normalized_mutual_info_score(true_labels[true_labels>=0], adata.obs['domain-1'][true_labels>=0])
    
    # Tied
    gmm = BayesianGaussianMixture(n_components=n_clusters, covariance_type='tied', random_state=42, init_params = 'random_from_data', n_init = n_init, max_iter = max_iter) # 5, 1000
    meta_labels  = gmm.fit_predict(meta_emb)
    # if refinement:  
    #     # print(meta_spatial.shape, meta_labels.shape)
    #     meta_labels  = cKD_refine_label(meta_spatial, meta_labels, k = radius)
    # Extrapolation
    adata.obs['domain'] = 0
    for i in range(len(clusterIDvec)):
        metacell = clusterIDvec[i]
        adata.obs.loc[cluster2Cellsdict[metacell], 'domain'] = meta_labels[i]
    if refinement:  
        adata.obs['domain'] = cKD_refine_label(np.array(adata.obsm['spatial']), adata.obs['domain'], k = radius)
    
    score_2 = normalized_mutual_info_score(true_labels[true_labels>=0], adata.obs['domain'][true_labels>=0])
    print(score_1, score_2)
    if score_1 > score_2:
        print('full is better than tied in this case')
        adata.obs['domain'] = adata.obs['domain-1']
        
def clustering_v3a3(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, covariance_type='full', n_init = 5, max_iter = 1000, linkage='ward'):
    print(covariance_type, linkage)
    print("Data input shape nx: ", adata.obsm[key].shape)
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    connectivity = 0.5 * (connectivity + connectivity.T)
    if 'pca' not in key:
        print("Run PCA")
        pca = PCA(n_components=n_components, random_state=42) 
        embedding = pca.fit_transform(np.asarray(adata.obsm[key].copy()))
        # adata.obsm['raw_pca'] = embedding
        embedding = connectivity.dot(embedding)
        # A = connectivity.dot(embedding)
        # B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
        # embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    else:
        embedding = connectivity.dot(adata.obsm[key])
    
    # hierarchical = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
    # initial_labels = hierarchical.fit_predict(embedding)
    
#     scaler = StandardScaler() 
#     X_scaled = scaler.fit_transform(embedding) 

#     # Normalizing the data so that the data 
#     # approximately follows a Gaussian distribution 
#     embedding = normalize(X_scaled)
    
    #gmm = GaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42)
    gmm = BayesianGaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42, init_params = 'random_from_data', n_init = n_init, max_iter = max_iter) # 5, 1000
    
    # gmm.means_init = np.array([embedding[initial_labels == i].mean(axis=0) for i in range(n_clusters)])
    gmm.fit(embedding)
    cluster_labels = gmm.predict(embedding)
    adata.obs['domain'] = cluster_labels
    
    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 


def clustering_v3a2Light(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, covariance_type='full', linkage='ward'):
    # print(covariance_type, linkage)
    # print("Data input shape: ", adata.obsm[key].shape)
    # pca = PCA(n_components=n_components, random_state=42) 
    # embedding = pca.fit_transform(adata.obsm[key].copy())
    # adata.obsm['raw_pca'] = embedding
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    # A = connectivity.dot(adata.obsm[key])
    embedding = connectivity.dot(adata.obsm[key])
    # n_cells = int(adata.shape[0])
    # B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    #B = sp.sparse.dia_matrix(1.0 / np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))), shape=(n_cells,n_cells));
    # B = sp.sparse.dia_matrix(1.0 / np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))), shape=(n_cells,n_cells));
    # embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    # embedding = A/connectivity.sum(axis=0)[0:10].mean()
    
    # hierarchical = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
    # initial_labels = hierarchical.fit_predict(embedding)
    
    #gmm = GaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42)
    gmm = BayesianGaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=0, init_params = 'random_from_data', n_init = 10, max_iter = 1000)
    # gmm.means_init = np.array([embedding[initial_labels == i].mean(axis=0) for i in range(n_clusters)])
    # print("dddd")
    # Step 4: Fit GMM
    gmm.fit(embedding)
    cluster_labels = gmm.predict(embedding)
    adata.obs['domain'] = cluster_labels
    
    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
        
######################################## Clustering Multimodal data #################################################
######################################## Clustering Multimodal data #################################################
### Good for multimodal dataset
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from sklearn.cluster import AgglomerativeClustering
def clustering_v3amulti(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, covariance_type='full', linkage='ward'):
    print(covariance_type, linkage)
    print("Data input shape: ", adata.obsm[key].shape)
    pca = PCA(n_components=n_components, random_state=42) 
    embedding = pca.fit_transform(adata.obsm[key].copy())
    adata.obsm['raw_pca'] = embedding
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    embedding = dot(B, A) #+ 0.5*embedding #dot(connectivity, embedding)
    
    hierarchical = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
    initial_labels = hierarchical.fit_predict(embedding)
    
    gmm = GaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42)
    # gmm = BayesianGaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42)
    gmm.means_init = np.array([embedding[initial_labels == i].mean(axis=0) for i in range(n_clusters)])
    # print("dddd")
    # Step 4: Fit GMM
    gmm.fit(embedding)
    cluster_labels = gmm.predict(embedding)
    adata.obs['domain'] = cluster_labels
    
    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type
        

from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import kneighbors_graph
def clustering_multi(adata, n_clusters=7, radius=50, key='emb', refinement=True, n_neighbors = 3, n_components = 20):
    # print("Data input shape: ", adata.obsm[key].shape)
#     pca = PCA(n_components=n_components, random_state=42) 
#     embedding = pca.fit_transform(adata.obsm[key].copy())
#     adata.obsm['raw_pca'] = embedding

    embedding = adata.obsm['RNA_feat']
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    A = connectivity.dot(embedding)
    B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    embedding_RNA = np.asarray(dot(B, A))
    
    embedding = adata.obsm['Pro_feat']
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    A = connectivity.dot(embedding)
    # B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    n_cells = int(adata.shape[0])
    B = sp.sparse.spdiags(1.0/connectivity.sum(axis=0), 0, n_cells, n_cells)
    embedding_Pro = np.asarray(dot(B, A))
    
    pca = PCA(n_components=2, random_state=42) 
    RNA_pca_embedding = pca.fit_transform(embedding_RNA)
    
    pca = PCA(n_components=2, random_state=42) 
    Pro_pca_embedding = pca.fit_transform(embedding_Pro)
    
    # Standardize the data
    scaler_a = StandardScaler()
    scaler_b = StandardScaler()

    data_a_train = scaler_a.fit_transform(embedding_RNA)
    data_b_train = scaler_b.fit_transform(embedding_Pro)

    # Define and train the CCA model
    n_components = 5  # Number of canonical components
    cca = CCA(n_components=n_components)
    cca.fit(data_a_train, data_b_train)

    # Transform the data into canonical space
    data_a_train_cca, data_b_train_cca = cca.transform(data_a_train, data_b_train)
    # adata.obsm['emb_pca'] = np.concatenate((data_a_train_cca, data_b_train_cca), axis=1)
    adata.obsm['emb_pca'] = np.concatenate((embedding_RNA, data_b_train_cca), axis=1)
    # adata.obsm['emb_pca'] = np.concatenate((data_a_train, data_b_train), axis=1)
    # adata.obsm['emb_pca'] = np.concatenate((embedding_RNA/np.linalg.norm(embedding_RNA), embedding_Pro/np.linalg.norm(embedding_Pro)), axis=1)
    # adata.obsm['emb_pca'] = np.asarray(embedding) 
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['domain'] = adata.obs['mclust']

    if refinement:  
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type 
        
def clustering_multiv2(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, covariance_type='full', linkage='ward'):
    embedding = adata.obsm['RNA_feat']
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    embedding_RNA = connectivity.dot(embedding)
    # A = connectivity.dot(embedding)
    # B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # n_cells = int(adata.shape[0])
    # B = sp.sparse.spdiags(1.0/connectivity.sum(axis=0), 0, n_cells, n_cells)
    # embedding_RNA = np.asarray(dot(B, A))
     
    embedding = adata.obsm['Pro_feat']
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    embedding_Pro = connectivity.dot(embedding)
    # A = connectivity.dot(embedding)
    # B = np.diag(1.0/np.squeeze(np.asarray(np.sum(connectivity.todense(), axis=1))))
    # n_cells = int(adata.shape[0])
    # B = sp.sparse.spdiags(1.0/connectivity.sum(axis=0), 0, n_cells, n_cells)
    # embedding_Pro = np.asarray(dot(B, A))
    
    pca = PCA(n_components=2, random_state=42) 
    RNA_pca_embedding = pca.fit_transform(embedding_RNA)
    
    pca = PCA(n_components=2, random_state=42) 
    Pro_pca_embedding = pca.fit_transform(embedding_Pro)
    
    # Standardize the data
    scaler_a = StandardScaler()
    scaler_b = StandardScaler()

    data_a_train = scaler_a.fit_transform(embedding_RNA)
    data_b_train = scaler_b.fit_transform(embedding_Pro)

    # Define and train the CCA model
    n_components = 5  # Number of canonical components
    cca = CCA(n_components=n_components)
    cca.fit(data_a_train, data_b_train)

    # Transform the data into canonical space
    data_a_train_cca, data_b_train_cca = cca.transform(data_a_train, data_b_train)
    adata.obsm['emb_pca'] = np.concatenate((data_a_train_cca, data_b_train_cca), axis=1)
    # adata.obsm['emb_pca'] = np.concatenate((embedding_RNA, data_b_train_cca), axis=1)
    
    # adata.obsm['emb_pca'] = np.concatenate((data_a_train, data_b_train), axis=1)
    # 0.38532702122083684
    # 0.41865301896421203
    
    # adata.obsm['emb_pca'] = np.concatenate((data_a_train/np.linalg.norm(data_a_train), 100*data_b_train), axis=1)
    
    # 0.3498349520666313, 0.44629947822423305
    # 0.6485536904822723, 0.4483355465783212

    # adata.obsm['emb_pca'] = np.concatenate((embedding_RNA/np.linalg.norm(embedding_RNA), embedding_Pro/np.linalg.norm(embedding_Pro)), axis=1)
    adata.obsm['emb_pca'] = np.concatenate((embedding_RNA/np.linalg.norm(embedding_RNA), embedding_Pro), axis=1)

    # adata.obsm['emb_pca'] = np.asarray(embedding) 
    
    gmm = BayesianGaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42, init_params = 'random_from_data', n_init = 5, max_iter = 1000)
    # gmm.means_init = np.array([embedding[initial_labels == i].mean(axis=0) for i in range(n_clusters)])
    # Step 4: Fit GMM
    gmm.fit(adata.obsm['emb_pca'])
    cluster_labels = gmm.predict(adata.obsm['emb_pca'])
    adata.obs['domain'] = cluster_labels

import mvlearn
from mvlearn.decomposition import GroupPCA
def clustering_multimvlearn(adata, n_clusters=7, radius=50, key='emb', refinement=False, n_neighbors = 3, n_components = 20, covariance_type='full', linkage='ward'):
    embedding = adata.obsm['RNA_feat']
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    embedding_RNA = connectivity.dot(embedding)
  
    embedding = adata.obsm['Pro_feat']
    connectivity = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    embedding_Pro = connectivity.dot(embedding)
    
    pca = PCA(n_components=2, random_state=42) 
    RNA_pca_embedding = pca.fit_transform(embedding_RNA)
    
    pca = PCA(n_components=2, random_state=42) 
    Pro_pca_embedding = pca.fit_transform(embedding_Pro)
    
    # Standardize the data
    scaler_a = StandardScaler()
    scaler_b = StandardScaler()

    data_a_train = scaler_a.fit_transform(embedding_RNA)
    data_b_train = scaler_b.fit_transform(embedding_Pro)

    # Define and train the CCA model
    n_components = 5  # Number of canonical components
    cca = CCA(n_components=n_components)
    cca.fit(data_a_train, data_b_train)

    #adata.obsm['emb_pca'] = np.concatenate((embedding_RNA/np.linalg.norm(embedding_RNA), embedding_Pro), axis=1)
    
    Xs = [data_a_train, data_b_train] # multiview data
    Xs_components = GroupPCA(multiview_output = False).fit_transform(Xs)
    print(Xs_components.shape)
    adata.obsm['emb_pca'] = Xs_components

    gmm = BayesianGaussianMixture(n_components=n_clusters, covariance_type=covariance_type, random_state=42, init_params = 'random_from_data', n_init = 5, max_iter = 1000)
    # gmm.means_init = np.array([embedding[initial_labels == i].mean(axis=0) for i in range(n_clusters)])
    # Step 4: Fit GMM
    gmm.fit(adata.obsm['emb_pca'])
    cluster_labels = gmm.predict(adata.obsm['emb_pca'])
    adata.obs['domain'] = cluster_labels
