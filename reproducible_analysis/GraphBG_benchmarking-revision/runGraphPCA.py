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

import glob
data_dir_full = []
for data_dir in glob.glob("/data/hoan/spatial_transcriptomics/data/SDMBench/Data/*.h5ad", recursive=True):
    # print(data_dir)
    data_dir_full.append(data_dir)
for data_dir in glob.glob("/data/hoan/spatial_transcriptomics/data/SDMBench/Data/*/*.h5ad", recursive=True):
    if '378k' not in data_dir:
        # print(data_dir)
        data_dir_full.append(data_dir)

output_folder = '/data/hoan/spatial_transcriptomics/GraphBG-analysis/output/graphpca/'
for data_dir in data_dir_full:
    # data_dir = data_dir_full[1]
    print(data_dir)
    data_name = data_dir.split('/')[-2] + data_dir.split('/')[-1][:-5]
    out_data = output_folder + data_name + '.csv'

    if 'Slice_3_removed' in data_dir:
        continue;
    
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
    
        adata.var_names_make_unique()
        sc.pp.filter_genes(adata, min_cells=20)
        sc.experimental.pp.normalize_pearson_residuals(adata)
        sc.pp.scale(adata)
        x_array=adata.obs["x"].tolist()
        y_array=adata.obs["y"].tolist()
        location=np.array([x_array, y_array]).T.astype(np.float32)
        Z,W = sg.Run_GPCA(adata, location=location, n_components = 30, method = "knn",  n_neighbors = 7, _lambda =0.5,
                        save_reconstruction=True)
        # cu la 20, ma_iter la 50
        adata.obsm["GraphPCA"] = Z
        
        estimator = KMeans(n_clusters=len(np.unique(true_labels)))
        res = estimator.fit(Z[:,:])
        lable_pred=res.labels_
        adata.obs["GPCA_pred"]= lable_pred
        adata.obs["GPCA_pred"] = adata.obs["GPCA_pred"].astype('category')
        refined_pred=sg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["GPCA_pred"].tolist(), dis= pair(location),
                               shape="generic",neighbor_num=7)
        adata.obs["GPCA_pred"]= refined_pred
        adata.obs["GPCA_pred"] = adata.obs["GPCA_pred"].astype('category')
        print(ari_score(adata.obs.GPCA_pred,adata.obs.Region))
        
        df_out = pd.DataFrame({'ground_truth': true_labels, 'GraphPCA': adata.obs.GPCA_pred})
        df_out.to_csv(out_data, index=False)