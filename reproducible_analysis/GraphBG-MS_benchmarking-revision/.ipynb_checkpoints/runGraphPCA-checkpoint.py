# python env biopy311
import GraphPCA as sg
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances as pair
from sklearn.metrics import adjusted_rand_score as ari_score

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

data_dir = '/data/hoan/spatial_transcriptomics/data/SDMBench/Data/378k/BrainAgingSpatialAtlas_MERFISH.h5ad'
adata=sc.read_h5ad(data_dir)
adata.obs['slice_id'] = adata.obs['donor_id'].astype(str) +'Slice'+ adata.obs['slice'].astype(str)
adata.obsm['spatial'] = adata.obsm['spatial_coords']
adata.obs['batch'] = adata.obs['donor_id'].astype(str) +'Slice'+ adata.obs['slice'].astype(str)
adata.obs['Region'] = adata.obs['tissue']
output_folder = '/data/hoan/spatial_transcriptomics/GraphBG-analysis/output/graphpca/'

import time
start_time = time.perf_counter()

from GraphPCA import Run_Hierarchical_Multi_GPCA
# adata = sc.read_h5ad('/path/to/MERFISH_Aging.h5ad')
sections = [adata[adata.obs['slice_id'].eq(name)].copy() for name in adata.obs['slice_id'].unique()]
locations = [one.obsm['spatial'][:, :2] for one in sections]
Zs, W0, Ws, info = Run_Hierarchical_Multi_GPCA(
    sections, locations, n_components=30, lambdas=0.5, rhos=0.14872615,
    sample_weights='equal_slice', max_iter=50, random_seed=666, return_info=True, n_neighbors = 7)

Z = np.vstack(Zs)
true_labels = adata.obs['Region'].values
estimator = KMeans(n_clusters=len(np.unique(true_labels)))
res = estimator.fit(Z)
lable_pred=res.labels_
print(normalized_mutual_info_score(true_labels, lable_pred))

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.4f} seconds")

df_out = pd.DataFrame({'ground_truth': true_labels, 'GraphPCA': lable_pred})
df_out.to_csv('output/GraphPCA_BrainAgingSpatialAtlas_v2.csv', index=False)