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

import scanpy as sc
from banksy.initialize_banksy import initialize_banksy
from banksy.run_banksy import run_banksy_multiparam
from banksy_utils.color_lists import spagcn_color

import glob
data_dir_full = []
for data_dir in glob.glob("/data/hoan/spatial_transcriptomics/data/SDMBench/Data/*.h5ad", recursive=True):
    # print(data_dir)
    data_dir_full.append(data_dir)
for data_dir in glob.glob("/data/hoan/spatial_transcriptomics/data/SDMBench/Data/*/*.h5ad", recursive=True):
    if '378k' not in data_dir:
        # print(data_dir)
        data_dir_full.append(data_dir)


import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.metrics import normalized_mutual_info_score

output_folder = '/data/hoan/spatial_transcriptomics/GraphBG-analysis/output/banksy/'

for data_dir in data_dir_full:
    print(data_dir)
    data_name = data_dir.split('/')[-2] + data_dir.split('/')[-1][:-5]
    out_data = output_folder + data_name + '.csv'

    # --- Check if out_data file exists ---
    if os.path.exists(out_data):
        print(f"File already exists: {out_data}. Skipping...")
    else:
        print(f"File not found. Proceeding with analysis: {out_data}")
        file_pickle = '/data/hoan/spatial_transcriptomics/GraphBGM/Data/Pickle/' + data_name + '.pickle'
        
        # Load and filter out NaN labels
        adata = sc.read_h5ad(data_dir)
        adata = adata[~adata.obs['Region'].isna()].copy()
        
        true_labels = pd.factorize(adata.obs['Region'])[0]
        adata.obs['x'] = adata.obsm['spatial'][:, 0]
        adata.obs['y'] = adata.obsm['spatial'][:, 1]
        
        coord_keys = ('x', 'y', 'spatial')
        pca_dims = [20] # number of dimensions to keep after PCA
        lambda_list = [.8] # lambda
        k_geom = 15 # spatial neighbours
        max_m = 1 # use AGF
        nbr_weight_decay = "scaled_gaussian" # can also be "reciprocal", "uniform" or "ranked"
        
        banksy_dict = initialize_banksy(
            adata,
            coord_keys,
            k_geom,
            nbr_weight_decay=nbr_weight_decay,
            max_m=max_m,
            plt_edge_hist=False,
            plt_nbr_weights=False,
            plt_agf_angles=False,
            plt_theta=False,
        )
        
        annotation_key = 'Region' # Key to access annotations in adata.obs[annotation_keys]
        num_clusters = len(set(true_labels))
        cluster_algorithm = 'leiden' 
        
        # --- Binary Search for Resolution ---
        low_res = 0.01
        high_res = 3.0
        best_res = 0.8
        best_labels = None
        max_iterations = 10
        
        print(f"Target number of clusters: {num_clusters}")
        
        for iteration in range(max_iterations):
            mid_res = (low_res + high_res) / 2.0
            resolutions = [mid_res]
            
            # Dynamically generate a large enough color list using mcolors
            dynamic_color_list = [plt.cm.tab20(i % 20) for i in range(100)]
            dynamic_color_list = [mcolors.to_hex(c) for c in dynamic_color_list]
            
            results_df = run_banksy_multiparam(
                adata,
                banksy_dict,
                lambda_list,
                resolutions,
                color_list=dynamic_color_list,
                max_m=max_m,
                filepath=output_folder,
                key=coord_keys,
                pca_dims=pca_dims,
                annotation_key=annotation_key,
                max_labels=num_clusters,
                cluster_algorithm=cluster_algorithm,
                match_labels=False,
                savefig=False,
                add_nonspatial=False,
                variance_balance=False,
            )
            
            params_name = results_df.index[0]
            tmp_label = results_df.loc[params_name, "labels"].dense
            n_clusters_found = len(np.unique(tmp_label))
            
            print(f"Iteration {iteration+1}: Resolution = {mid_res:.4f} -> Clusters found = {n_clusters_found}")
            
            best_res = mid_res
            best_labels = tmp_label
            
            if n_clusters_found == num_clusters:
                print("Exact cluster match achieved!")
                break
            elif n_clusters_found < num_clusters:
                low_res = mid_res
            else:
                high_res = mid_res
        
        # Final evaluation using the best found labels
        tmp_label = best_labels
        print("NMI = ", normalized_mutual_info_score(tmp_label, true_labels))
        print(num_clusters, len(np.unique(tmp_label)))
        print('=='*15)
        
        df_out = pd.DataFrame({'ground_truth': true_labels, 'BANKSY': tmp_label})
        df_out.to_csv(out_data, index=False)