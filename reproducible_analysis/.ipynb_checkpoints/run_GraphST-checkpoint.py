#!/usr/bin/env python
# coding: utf-8

from numpy.random import default_rng
import scanpy as sc
# import squidpy as sq
from anndata import AnnData
import scipy
sc.logging.print_header()
from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np
import pandas as pd
import seaborn as sns
import os
import torch
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphST import GraphST

import glob
data_dir_full = []
for data_dir in glob.glob("/data/hoan/spatial_transcriptomics/data/SDMBench/Data/*.h5ad", recursive=True):
    # print(data_dir)
    data_dir_full.append(data_dir)
for data_dir in glob.glob("/data/hoan/spatial_transcriptomics/data/SDMBench/Data/*/*.h5ad", recursive=True):
    if '378k' not in data_dir:
        # print(data_dir)
        data_dir_full.append(data_dir)

for data_dir in data_dir_full:
    print(data_dir)
    data_name = data_dir.split('/')[-2]+data_dir.split('/')[-1][:-5]
    file_pickle = '/data/hoan/spatial_transcriptomics/GraphBGM/Data/Pickle/'+data_name+'.pickle'
    adata=sc.read_h5ad(data_dir)
    true_labels = pd.factorize(adata.obs['Region'])[0]


    # Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    os.environ['R_HOME'] = '/data/hoan/mybin/miniconda3/lib/R'


    if np.min(true_labels) >= 0:
        n_clusters = len(np.unique(true_labels))
    else:
        n_clusters = len(np.unique(true_labels)) - 1

    import pickle
    import os.path
    if os.path.exists(file_pickle):
        # Reload the file
        adata = pickle.load(open(file_pickle, "rb"))
    else:
        # define model
        model = GraphST.GraphST(adata, device=device)
        # train model
        adata = model.train()
        # Save the file
        pickle.dump(adata, file = open(file_pickle, "wb"))


    if scipy.sparse.issparse(adata.X):
        adata.obsm['raw'] = adata.X.todense() 
    else:
        adata.obsm['raw'] = adata.X

    from sklearn.decomposition import PCA
    from GraphST.utils import refine_label
    from sklearn.preprocessing import StandardScaler
    from GraphST.utils import mclust_R

    # Run GraphST
    from GraphST.utils import clustering
    clustering(adata, n_clusters=n_clusters, radius=50, key='emb', refinement=True)
    adata.obs['GraphST'] = adata.obs['domain']
    print('GraphST: ', normalized_mutual_info_score(true_labels[true_labels>=0], adata.obs['domain'][true_labels>=0]))


    
