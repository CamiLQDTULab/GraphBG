#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scanpy as sc
import umap
import INSPIRE
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import anndata as ad
import scipy.sparse
import os
from sklearn.metrics import pairwise_distances
from sklearn.metrics import silhouette_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from annoy import AnnoyIndex

data_dir = '/data/hoan/spatial_transcriptomics/data/SDMBench/Data/378k/BrainAgingSpatialAtlas_MERFISH.h5ad'
adata=sc.read_h5ad(data_dir)

adata.obs["slice_stage"] = [adata.obs["slice"][i] + "_" + adata.obs["development_stage"][i] for i in range(adata.shape[0])]
adata.obs["slice_stage_donor"] = [adata.obs["slice_stage"][i] + "_" + adata.obs["donor_id"][i] for i in range(adata.shape[0])]
adata.obs['Region'] = adata.obs['tissue']
adata.X = adata.raw.X
adata.obsm["spatial"] = adata.obsm["X_spatial_coords"]
set(adata.obs["slice_stage_donor"].values)


import time
start_time = time.perf_counter()

adata_st_list = []
section_id_list = list(sorted(set(adata.obs["slice_stage_donor"].values)))
for i in range(len(section_id_list)):
    print(section_id_list[i])
    adata_st_list.append(adata[adata.obs["slice_stage_donor"] == section_id_list[i], :].copy())

adata_st_list, adata_full = INSPIRE.utils.preprocess(adata_st_list=adata_st_list,
                                                     num_hvgs=1000,
                                                     min_genes_qc=1,
                                                     min_cells_qc=1,
                                                     spot_size=20,
                                                     limit_num_genes=True)

adata_st_list = INSPIRE.utils.build_graph_LGCN(adata_st_list=adata_st_list,
                                               rad_cutoff_list=[60]*len(adata_st_list))


model = INSPIRE.model.Model_LGCN(adata_st_list=adata_st_list,
                                 n_spatial_factors=30, #may need to choose a suitable number by users
                                 n_training_steps=10000,
                                 batch_size=2048,
                                 different_platforms=False
                                )
model.train(adata_st_list)

adata_full, basis_df = model.eval(adata_st_list, adata_full)
basis = np.array(basis_df.values)

## Note from INSPIRE: https://inspire-tutorial.readthedocs.io/en/latest/examples/INSPIRE_LGCN.html#
# The learned spot representations enable the identification of spatial domains. 
#When the number of spatial domains is known a priori, clustering is performed using a Gaussian Mixture Model implemented in 
#the sklearn package. In the absence of such prior knowledge, spatial domains are identified using the Louvain clustering 
#algorithm, implemented through the sc.pp.neighbors and sc.tl.louvain functions in the Scanpy package


import scanpy as sc
from sklearn.mixture import GaussianMixture

n_domains = 8

gmm = GaussianMixture(n_components=n_domains, random_state=42)
gmm_labels = gmm.fit_predict(adata_full.obsm["latent"])

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.4f} seconds")

true_labels = adata_full.obs['Region'].values
print(normalized_mutual_info_score(true_labels, gmm_labels))
df_out = pd.DataFrame({'ground_truth': true_labels, 'INSPIRE': gmm_labels})
df_out.to_csv('/data/hoan/spatial_transcriptomics/GraphBG-analysis/output/INSPIRE_BrainAgingSpatialAtlas.csv', index=False)

