#!/usr/bin/env python
# coding: utf-8
import scanpy as sc
from numpy.random import default_rng
from anndata import AnnData
import scipy
from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np
import pandas as pd
import seaborn as sns
import os
from sklearn import metrics
import multiprocessing as mp
from sklearn.metrics.cluster import normalized_mutual_info_score
import pickle
from sklearn.neighbors import kneighbors_graph
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from numpy import genfromtxt
import time

data_dir = '/data/hoan/spatial_transcriptomics/data/SDMBench/Data/378k/BrainAgingSpatialAtlas_MERFISH.h5ad'
adata_big=sc.read_h5ad(data_dir)
donor_id = np.unique(adata_big.obs['donor_id'])

# ## ----------------------------------- SpaceFlow Training runtime ----------------------------------------------------
# from SpaceFlow import SpaceFlow
# import squidpy as sq
# sf_embedding_all = np.array([], dtype=np.int64).reshape(0,50)
# batch_labels = []
# batch_index = 0
# sf_meta = []
# for sliceID in ['0', '1', '2']:
#     for donorID in donor_id:
#         print(sliceID, donorID)
#         adata = adata_big[adata_big.obs['slice']==sliceID]#[:5000]
#         adata = adata[adata.obs['donor_id']==donorID]
#         adata.obsm['spatial'] = adata.obsm['spatial_coords']
#         true_labels = pd.factorize(adata.obs['tissue'])[0]
#         if adata.shape[0] > 10:
#             time_st = time.time()
#             data_name = donorID + 'Slice' + sliceID
#             adataNew = AnnData(1.0 + adata.obsm['X_pca'] - adata.obsm['X_pca'].min())
#             adataNew.obsm['spatial'] = adata.obsm['spatial']
#             sf = SpaceFlow.SpaceFlow(adata=adataNew)
#             sf.preprocessing_data(n_top_genes=40)
#             embed_filepath = 'output/'+data_name+'SpaceFlowembedding.tsv'
#             domain_label = 'output/'+data_name+'SpaceFlowdomains.tsv'
#             sf.train(embedding_save_filepath=embed_filepath, spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=1000, max_patience=50, min_stop=100, random_seed=42, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)
            
#             # sf.segmentation(domain_label_save_filepath=domain_label, n_neighbors=50, resolution=1.0)
#             sf_embedding = genfromtxt(embed_filepath, delimiter='\t')
#             # embedding_adata = AnnData(sf_embedding)
#             # sc.pp.neighbors(embedding_adata, n_neighbors=50, use_rep='X')
#             # sc.tl.leiden(embedding_adata, resolution=1.0)
#             # domains = embedding_adata.obs["leiden"].cat.codes
#             time_ed = time.time()
#             time_spaceflow = time_ed-time_st

#             # cluster_label = pd.read_csv(domain_label, header = None)
#             # cluster_SpaceFlow = cluster_label.values.flatten()
#             # cluster_SpaceFlow = np.array(domains)
#             # df_out = pd.DataFrame()
#             # df_out['SpaceFLow'] = np.array(cluster_SpaceFlow[true_labels>=0]).astype(int)
#             # df_out['ground_truth'] = np.array(true_labels[true_labels>=0]).astype(int)
#             # out_data = '/data/hoan/spatial_transcriptomics/GraphST/output/'+data_name+'SpaceFlow.csv'
#             # df_out.to_csv(out_data, index=False)
#             # time_st = time.time()
#             # sf_embedding = genfromtxt(embed_filepath, delimiter='\t')
#             # sf_embedding_all = np.vstack([sf_embedding_all, sf_embedding])
#             # batch_labels += [batch_index]*sf_embedding.shape[0]
#             # batch_index += 1
#             # time_ed = time.time()
#             # time_combine_preprocess = time_ed-time_st
#             # sf_info = [data_name, sf_embedding.shape[0], normalized_mutual_info_score(true_labels, cluster_SpaceFlow), time_spaceflow, time_combine_preprocess]
#             sf_info = [data_name, sf_embedding.shape[0], time_spaceflow]
#             sf_meta.append(sf_info)
#             # break;

# sf_meta_df = pd.DataFrame(sf_meta, columns = ['metaSlice', 'n', 'SpaceFlowTrainRuntime'])
# sf_meta_df.to_csv('sf_meta_df_runtime_round1.csv')

# # sf_meta_df = pd.DataFrame(sf_meta, columns = ['metaSlice', 'n', 'NMI', 'SpaceFlowRuntime', 'Hamony_prepro_runtime'])
# # sf_meta_df.to_csv('sf_meta_df_runtime_round1.csv')

# # import scanpy.external as sce
# # time_st = time.time()
# # adata_combine = sc.AnnData(sf_embedding_all)
# # adata_combine.obsm['X_pca'] = sf_embedding_all
# # adata_combine.obs["batch"] = np.array(batch_labels)
# # adata_combine.obs["batch"] = adata_combine.obs["batch"].astype("category")
# # # sc.pp.combat(adata_combine, key="batch", inplace=True)
# # # sce.pp.scanorama_integrate(adata, key="batch")
# # # sc.pp.pca(adata_combine, n_comps=30)
# # sce.pp.harmony_integrate(adata_combine, key="batch")
# # time_ed = time.time()
# # time_harmony = time_ed-time_st

# # time_st = time.time()
# # sc.pp.neighbors(adata_combine,use_rep='X_pca_harmony')
# # sc.tl.leiden(adata_combine,resolution=0.5)
# # time_ed = time.time()
# # time_leiden = time_ed-time_st

# # sf_meta_df['Harmony_combine_Time'] = time_harmony
# # sf_meta_df['time_leiden'] =time_leiden
# # sf_meta_df.to_csv('sf_meta_df_runtime.csv')



## ----------------------------------- SpaceFlow clustering runtime ----------------------------------------------------
sf_embedding_all = np.array([], dtype=np.int64).reshape(0,50)
batch_labels = []
batch_index = 0
sf_meta = []
for sliceID in ['0', '1', '2']:
    for donorID in donor_id:
        print(sliceID, donorID)
        adata = adata_big[adata_big.obs['slice']==sliceID]#[:5000]
        adata = adata[adata.obs['donor_id']==donorID]
        adata.obsm['spatial'] = adata.obsm['spatial_coords']
        true_labels = pd.factorize(adata.obs['tissue'])[0]
        if adata.shape[0] > 10:
            data_name = donorID + 'Slice' + sliceID
#             adataNew = AnnData(1.0 + adata.obsm['X_pca'] - adata.obsm['X_pca'].min())
#             adataNew.obsm['spatial'] = adata.obsm['spatial']
#             sf = SpaceFlow.SpaceFlow(adata=adataNew)
#             sf.preprocessing_data(n_top_genes=40)
            embed_filepath = 'output/'+data_name+'SpaceFlowembedding.tsv'
#             domain_label = 'output/'+data_name+'SpaceFlowdomains.tsv'
#             sf.train(embedding_save_filepath=embed_filepath, spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=1000, max_patience=50, min_stop=100, random_seed=42, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)
            
            # sf.segmentation(domain_label_save_filepath=domain_label, n_neighbors=50, resolution=1.0)
            sf_embedding = genfromtxt(embed_filepath, delimiter='\t')
            time_st = time.time()
            embedding_adata = AnnData(sf_embedding)
            sc.pp.neighbors(embedding_adata, n_neighbors=50, use_rep='X')
            sc.tl.leiden(embedding_adata, resolution=1.0)
            domains = embedding_adata.obs["leiden"].cat.codes
            time_ed = time.time()
            time_spaceflow = time_ed-time_st

            cluster_SpaceFlow = np.array(domains)
            df_out = pd.DataFrame()
            df_out['SpaceFLow'] = np.array(cluster_SpaceFlow[true_labels>=0]).astype(int)
            df_out['ground_truth'] = np.array(true_labels[true_labels>=0]).astype(int)
            out_data = '/data/hoan/spatial_transcriptomics/GraphST/output/'+data_name+'SpaceFlow.csv'
            df_out.to_csv(out_data, index=False)
            time_st = time.time()
            # sf_embedding = genfromtxt(embed_filepath, delimiter='\t')
            sf_embedding_all = np.vstack([sf_embedding_all, sf_embedding])
            batch_labels += [batch_index]*sf_embedding.shape[0]
            batch_index += 1
            time_ed = time.time()
            time_combine_preprocess = time_ed-time_st
            sf_info = [data_name, sf_embedding.shape[0], time_spaceflow, time_combine_preprocess,
                       normalized_mutual_info_score(true_labels, cluster_SpaceFlow)]
            # sf_info = [data_name, sf_embedding.shape[0], time_spaceflow]
            sf_meta.append(sf_info)
            print(np.array(sf_info))

# sf_meta_df = pd.DataFrame(sf_meta, columns = ['metaSlice', 'n', 'SpaceFlowTrainRuntime'])
# sf_meta_df.to_csv('sf_meta_df_runtime_round2_clustering.csv')

sf_meta_df = pd.DataFrame(sf_meta, columns = ['metaSlice', 'n', 'SpaceFlowRuntime', 'Harmony_prepro_runtime', 'NMI'])
sf_meta_df.to_csv('sf_meta_df_runtime_round2_clustering.csv')

import scanpy.external as sce
time_st = time.time()
adata_combine = sc.AnnData(sf_embedding_all)
adata_combine.obsm['X_pca'] = sf_embedding_all
adata_combine.obs["batch"] = np.array(batch_labels)
adata_combine.obs["batch"] = adata_combine.obs["batch"].astype("category")
# sc.pp.combat(adata_combine, key="batch", inplace=True)
# sce.pp.scanorama_integrate(adata, key="batch")
# sc.pp.pca(adata_combine, n_comps=30)
sce.pp.harmony_integrate(adata_combine, key="batch")
time_ed = time.time()
time_harmony = time_ed-time_st

time_st = time.time()
sc.pp.neighbors(adata_combine,use_rep='X_pca_harmony')
sc.tl.leiden(adata_combine,resolution=0.5)
time_ed = time.time()
time_leiden = time_ed-time_st

sf_meta_df['Harmony_combine_Time'] = time_harmony
sf_meta_df['time_leiden'] =time_leiden
sf_meta_df.to_csv('sf_meta_df_runtime_clusteringstep_final.csv')




