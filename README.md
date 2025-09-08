# GraphBG

**Fast Bayesian Domain Detection via Spectral Graph Convolutions for Multi-slice and Multimodal Spatial Transcriptomics**

GraphBG is a unified and scalable framework for spatial domain detection in spatial transcriptomics (ST) and spatial multi-omics data.
It integrates **spectral graph convolutions** with a **variational Bayesian Gaussian mixture model (VB-GMM)** to enable accurate, robust, and computationally efficient clustering of spatially coherent domains.

GraphBG extends to:

* **GraphBG-MS** → for multi-slice analysis with metacell aggregation, batch correction, and joint clustering.
* **GraphBG-MM** → for multimodal analysis with modality-specific graph embeddings and kernel canonical correlation analysis.

---

##  Key Features

* Scalable to datasets with **hundreds of thousands of cells**.
* Supports **unimodal, multi-slice, and multimodal** ST data.
* Achieves state-of-the-art accuracy across multiple benchmarks.
* Runs large datasets (e.g., **>370,000 cells in 5 minutes**) on standard hardware.
* Open-source and easy to integrate with existing ST analysis workflows.

---

## Installation

Clone the repository and install dependencies:

```bash
git clone git@github.com:CamiLQDTULab/GraphBG.git
cd GraphBG
pip install -r requirements.txt
```

---

## Usage

### 1. Unimodal ST analysis

```python
runGraphBG.ipynb
```

### 2. Multi-slice analysis (GraphBG-MS)

```python
runGraphBG-MS.ipynb
```

### 3. Multimodal ST analysis (GraphBG-MM)

```python
runGraphBG-MM.ipynb
```

## Reproducibility
```python
GraphBG/reproducible_analysis
```



## Citation

If you use GraphBG in your research, please cite:

> Do, V.H., Tran, T.P.L., & Canzar, S.
> **GraphBG: Fast Bayesian Domain Detection via Spectral Graph Convolutions for Multi-slice and Multimodal Spatial Transcriptomics**
> *Preprint*, 2025.

---

## License

GraphBG is released under the **BSD 3-Clause License**.

---
