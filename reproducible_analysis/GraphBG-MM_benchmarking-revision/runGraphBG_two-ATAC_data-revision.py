# ============================================================
# GraphBG  multimodal spatial domain analysis
#
# Dataset:
#   Mouse brain
#   RNA-seq + CUT&Tag H3K4me3
#
# Main pipeline:
#
# RNA:
#   QC -> HVG -> normalization -> log1p -> scaling -> PCA
#
# CUT&Tag:
#   TF-IDF -> LSI

import os
import itertools
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc

from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import kneighbors_graph
from sklearn.mixture import BayesianGaussianMixture
from sklearn.decomposition import TruncatedSVD

from esda.moran import Moran
from libpysal.weights import W

# IMPORTANT:
# KMCCA comes from mvlearn, NOT sklearn.
from mvlearn.embed import KMCCA

from SpatialGlue.preprocess import fix_seed
from utils import cKD_refine_label


warnings.filterwarnings("ignore")

RANDOM_SEED = 2022

fix_seed(RANDOM_SEED)

np.random.seed(RANDOM_SEED)


# ============================================================
# INPUT FILES
# ============================================================

# file_fold = (
#     "/data/hoan/spatial_transcriptomics/GraphBG-analysis/"
#     "Data/SpatialGlue/Data_SpatialGlue/"
#     "Dataset8_Mouse_Brain_H3K4me3"
# )

file_fold = (
    "/data/hoan/spatial_transcriptomics/GraphBG-analysis/"
    "Data/SpatialGlue/Data_SpatialGlue/"
    "Dataset7_Mouse_Brain_ATAC"
)

RNA_FILE = os.path.join(
    file_fold,
    "adata_RNA.h5ad"
)

CUTTAG_FILE = os.path.join(
    file_fold,
    "adata_peaks_normalized.h5ad"
)

adata_rna = sc.read_h5ad(
    RNA_FILE
)

adata_h3k4 = sc.read_h5ad(
    CUTTAG_FILE
)


adata_rna.var_names_make_unique()
adata_h3k4.var_names_make_unique()

# ============================================================
# MATCH CELLS BEFORE QC
# ============================================================

common_cells = (
    adata_rna.obs_names
    .intersection(
        adata_h3k4.obs_names
    )
)

print(
    "Common cells before QC:",
    len(common_cells)
)


adata_rna = (
    adata_rna[
        common_cells
    ].copy()
)

adata_h3k4 = (
    adata_h3k4[
        common_cells
    ].copy()
)


# Make absolutely sure the ordering is identical
adata_h3k4 = (
    adata_h3k4[
        adata_rna.obs_names
    ].copy()
)


assert np.array_equal(
    adata_rna.obs_names,
    adata_h3k4.obs_names
), "RNA and CUT&Tag cell ordering is not identical."


print(
    "Matched RNA:",
    adata_rna.shape
)

print(
    "Matched CUT&Tag:",
    adata_h3k4.shape
)


# ============================================================
# RNA PREPROCESSING
# ============================================================

print("\n" + "=" * 70)
print("RNA preprocessing")
print("=" * 70)


# ------------------------------------------------------------
# 4.1 Filter genes
# ------------------------------------------------------------

sc.pp.filter_genes(
    adata_rna,
    min_cells=10
)


# ------------------------------------------------------------
# 4.2 Filter cells
# ------------------------------------------------------------

sc.pp.filter_cells(
    adata_rna,
    min_genes=200
)

common_cells_after_qc = (
    adata_rna.obs_names
    .intersection(
        adata_h3k4.obs_names
    )
)

adata_rna = (
    adata_rna[
        common_cells_after_qc
    ].copy()
)

adata_h3k4 = (
    adata_h3k4[
        common_cells_after_qc
    ].copy()
)


# Force identical order again
adata_h3k4 = (
    adata_h3k4[
        adata_rna.obs_names
    ].copy()
)


assert (
    adata_rna.n_obs
    ==
    adata_h3k4.n_obs
), (
    "RNA and CUT&Tag have different "
    "numbers of cells after QC."
)


assert np.array_equal(
    adata_rna.obs_names,
    adata_h3k4.obs_names
), (
    "RNA and CUT&Tag cell ordering "
    "differs after QC."
)


print(
    "Cells after RNA QC:",
    adata_rna.n_obs
)


# ------------------------------------------------------------
# 4.3 Highly variable genes
# ------------------------------------------------------------

sc.pp.highly_variable_genes(
    adata_rna,
    flavor="seurat_v3",
    n_top_genes=3000,
    subset=True
)


print(
    "RNA genes after HVG:",
    adata_rna.n_vars
)


# ------------------------------------------------------------
# 4.4 Normalize
# ------------------------------------------------------------

sc.pp.normalize_total(
    adata_rna,
    target_sum=1e4
)


sc.pp.log1p(
    adata_rna
)


# ------------------------------------------------------------
# 4.5 Scale
# ------------------------------------------------------------

sc.pp.scale(
    adata_rna,
    max_value=10
)


# ------------------------------------------------------------
# 4.6 PCA
# ------------------------------------------------------------

sc.tl.pca(
    adata_rna,
    n_comps=50,
    svd_solver="arpack",
    random_state=RANDOM_SEED
)


rna_feat = (
    adata_rna
    .obsm["X_pca"]
    .copy()
)


print(
    "RNA feature:",
    rna_feat.shape
)


# ============================================================
# 5. CUT&TAG PREPROCESSING
# ============================================================

print("\n" + "=" * 70)
print("CUT&Tag H3K4me3 preprocessing")
print("=" * 70)


X = adata_h3k4.X


if sp.issparse(X):

    X = X.tocsr()

else:

    X = sp.csr_matrix(
        X
    )


print(
    "CUT&Tag matrix:",
    X.shape
)


# ------------------------------------------------------------
# 5.1 TF-IDF
# ------------------------------------------------------------

# Term frequency
cell_sum = (
    np.asarray(
        X.sum(axis=1)
    )
    .ravel()
)

cell_sum[
    cell_sum == 0
] = 1.0


tf = X.multiply(
    1.0 /
    cell_sum[:, None]
)


# ------------------------------------------------------------
# 5.2 Inverse document frequency
# ------------------------------------------------------------

n_cells = X.shape[0]


binary_X = X.copy()

binary_X.data = np.ones_like(
    binary_X.data
)


n_cells_feature = (
    np.asarray(
        binary_X.sum(axis=0)
    )
    .ravel()
)


idf = np.log(
    1.0
    +
    n_cells /
    (
        n_cells_feature
        +
        1.0
    )
)


tfidf = tf.multiply(
    idf
)


# ------------------------------------------------------------
# 5.3 LSI
#
# We calculate 51 components and remove the first component,
# giving 50 usable dimensions.
# ------------------------------------------------------------

n_lsi_total = 51


svd = TruncatedSVD(
    n_components=n_lsi_total,
    random_state=RANDOM_SEED
)


lsi_all = (
    svd.fit_transform(
        tfidf
    )
)


# Remove first LSI component
lsi_feat = (
    lsi_all[:, 1:51]
)


print(
    "CUT&Tag LSI:",
    lsi_feat.shape
)


# ============================================================
# 6. FINAL ALIGNMENT CHECK
# ============================================================

print("\n" + "=" * 70)
print("Final alignment check")
print("=" * 70)


assert (
    rna_feat.shape[0]
    ==
    lsi_feat.shape[0]
), (
    "RNA and CUT&Tag feature matrices "
    "have different numbers of cells."
)


assert (
    rna_feat.shape[0]
    ==
    adata_rna.n_obs
), (
    "RNA feature matrix does not match "
    "AnnData cell count."
)


assert (
    lsi_feat.shape[0]
    ==
    adata_h3k4.n_obs
), (
    "CUT&Tag feature matrix does not match "
    "AnnData cell count."
)


print(
    "Number of cells:",
    rna_feat.shape[0]
)

print(
    "RNA dimensions:",
    rna_feat.shape[1]
)

print(
    "CUT&Tag dimensions:",
    lsi_feat.shape[1]
)


# ============================================================
# 7. SPATIAL COORDINATES
# ============================================================

coords = np.asarray(
    adata_rna.obsm["spatial"]
)[:, :2]


assert (
    coords.shape[0]
    ==
    rna_feat.shape[0]
)


print(
    "Spatial coordinates:",
    coords.shape
)


# ============================================================
# 8. SPATIAL GRAPH
# ============================================================

def build_spatial_graph(
    coords,
    n_neighbors
):
    """
    Build a symmetric row-normalized kNN spatial graph.

    Parameters
    ----------
    coords : ndarray
        Spatial coordinates.

    n_neighbors : int
        Number of spatial neighbors.

    Returns
    -------
    A : scipy sparse matrix
        Row-normalized spatial graph.
    """

    n_cells = coords.shape[0]

    if n_neighbors >= n_cells:

        raise ValueError(
            "n_neighbors must be smaller "
            "than the number of cells."
        )


    A = kneighbors_graph(
        coords,
        n_neighbors=n_neighbors,
        mode="connectivity",
        include_self=False
    )


    # Symmetrize
    A = (
        0.5
        *
        (
            A
            +
            A.T
        )
    )


    # Row normalization
    row_sum = np.asarray(
        A.sum(axis=1)
    ).ravel()


    row_sum[
        row_sum == 0
    ] = 1.0


    D_inv = sp.diags(
        1.0 /
        row_sum
    )


    A = (
        D_inv
        @
        A
    )


    return A.tocsr()


# ============================================================
# 9. SPATIAL SMOOTHING
# ============================================================

def spatial_smooth(
    X,
    A,
    alpha
):
    """
    Spatial smoothing:

        X_smooth =
            (1-alpha) * X
            +
            alpha * A X

    alpha = 0:
        no spatial smoothing

    alpha = 1:
        completely neighbor-based
    """

    X_smooth = (
        (1.0 - alpha)
        *
        X
        +
        alpha
        *
        A.dot(X)
    )


    return np.asarray(
        X_smooth
    )


# ============================================================
# 10. STANDARDIZATION
# ============================================================

def standardize_features(
    X
):
    """
    Standardize features and remove
    zero-variance dimensions.
    """

    X = np.asarray(
        X,
        dtype=np.float64
    )


    X = np.nan_to_num(
        X,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )


    variance = X.var(
        axis=0
    )


    keep = (
        variance
        >
        1e-8
    )


    if keep.sum() < 2:

        raise ValueError(
            "Too few non-zero variance "
            "features after smoothing."
        )


    X = X[:, keep]


    scaler = StandardScaler()


    X = scaler.fit_transform(
        X
    )


    X = np.nan_to_num(
        X,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )


    return X


# ============================================================
# 11. KMCCA
# ============================================================

def run_kmcca(
    embedding_rna,
    embedding_cuttag,
    n_components,
    kernel="poly",
    regs=1e-3
):
    """
    Robust KMCCA wrapper.

    Important:
    - no X - X.min()
    - zero-variance features removed
    - effective rank checked
    - n_components automatically restricted
    """

    X_rna = standardize_features(
        embedding_rna
    )

    X_cuttag = standardize_features(
        embedding_cuttag
    )


    # --------------------------------------------------------
    # Check dimensions
    # --------------------------------------------------------

    if X_rna.shape[0] != X_cuttag.shape[0]:

        raise ValueError(
            "KMCCA input sample mismatch: "
            f"RNA={X_rna.shape}, "
            f"CUT&Tag={X_cuttag.shape}"
        )


    # --------------------------------------------------------
    # Numerical check
    # --------------------------------------------------------

    if not np.isfinite(
        X_rna
    ).all():

        raise ValueError(
            "RNA contains non-finite values."
        )


    if not np.isfinite(
        X_cuttag
    ).all():

        raise ValueError(
            "CUT&Tag contains non-finite values."
        )


    # --------------------------------------------------------
    # Effective rank
    #
    # Smoothing can reduce effective rank.
    # --------------------------------------------------------

    rank_rna = np.linalg.matrix_rank(
        X_rna
    )

    rank_cuttag = np.linalg.matrix_rank(
        X_cuttag
    )


    max_components = min(
        rank_rna,
        rank_cuttag,
        X_rna.shape[1],
        X_cuttag.shape[1]
    )


    if max_components < 2:

        raise ValueError(
            "Insufficient rank for KMCCA: "
            f"RNA rank={rank_rna}, "
            f"CUT&Tag rank={rank_cuttag}"
        )


    n_components_used = min(
        n_components,
        max_components
    )


    # --------------------------------------------------------
    # KMCCA
    # --------------------------------------------------------

    Xs = [
        X_rna,
        X_cuttag
    ]


    kmcca = KMCCA(
        n_components=n_components_used,
        kernel=kernel,
        regs=regs
    )


    kmcca.fit(
        Xs
    )


    Xs_latents = kmcca.transform(
        Xs
    )


    # --------------------------------------------------------
    # Validate KMCCA output
    # --------------------------------------------------------

    if len(Xs_latents) != 2:

        raise ValueError(
            "KMCCA returned "
            f"{len(Xs_latents)} views instead of 2."
        )


    Z_rna = np.asarray(
        Xs_latents[0]
    )


    Z_cuttag = np.asarray(
        Xs_latents[1]
    )


    if (
        Z_rna.shape[0]
        !=
        Z_cuttag.shape[0]
    ):

        raise ValueError(
            "KMCCA output dimension mismatch: "
            f"RNA={Z_rna.shape}, "
            f"CUT&Tag={Z_cuttag.shape}"
        )


    # --------------------------------------------------------
    # Combine KMCCA views
    # --------------------------------------------------------

    embedding = np.concatenate(
        [
            Z_rna,
            Z_cuttag
        ],
        axis=1
    )


    embedding = np.nan_to_num(
        embedding,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )


    embedding = StandardScaler().fit_transform(
        embedding
    )


    return (
        embedding,
        rank_rna,
        rank_cuttag,
        n_components_used
    )


# ============================================================
# 12. PY-SAL MORAN WEIGHTS
# ============================================================

def scipy_to_pysal(
    A
):
    """
    Convert scipy sparse spatial graph
    to libpysal W.
    """

    A = A.tocsr()


    neighbors = {}
    weights = {}


    for i in range(
        A.shape[0]
    ):

        start = A.indptr[i]
        end = A.indptr[i + 1]


        idx = A.indices[
            start:end
        ]

        val = A.data[
            start:end
        ]


        neighbors[i] = (
            idx.tolist()
        )

        weights[i] = (
            val.tolist()
        )


    return W(
        neighbors,
        weights
    )


# ============================================================
# 13. GMM
# ============================================================

def run_gmm(
    embedding,
    gmm_components,
    covariance_type
):

    gmm = BayesianGaussianMixture(

        n_components=gmm_components,

        covariance_type=covariance_type,

        random_state=RANDOM_SEED,

        init_params="kmeans",

        n_init=3,

        max_iter=1000,

        reg_covar=1e-6
    )


    gmm.fit(
        embedding
    )


    labels = gmm.predict(
        embedding
    )


    return labels

### Run the main method

n_neighbors = 4
alpha_rna = 0.2
alpha_cuttag = 0.2
kmcca_regs = 1e-2

GMM_COMPONENTS = 18
COVARIANCE = "diag"
REFINE_K = 30
kmcca_components = 30

A = build_spatial_graph(
    coords,
    n_neighbors
)

rna_spatial = spatial_smooth(
    rna_feat,
    A,
    alpha_rna
)


cuttag_spatial = spatial_smooth(
    lsi_feat,
    A,
    alpha_cuttag
)

embedding, rank_rna, rank_cuttag, n_used = (
    run_kmcca(
        rna_spatial,
        cuttag_spatial,
        n_components=kmcca_components,
        kernel="poly",
        regs=kmcca_regs
    )
)

labels = run_gmm(
    embedding,
    GMM_COMPONENTS,
    COVARIANCE
)

refined_labels = cKD_refine_label(
    coords,
    labels,
    k=REFINE_K
)

refined_labels = np.asarray(
    refined_labels
)
