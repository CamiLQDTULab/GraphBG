# GraphBG

**Fast Bayesian Domain Detection via Spectral Graph Convolutions for Multi-slice and Multimodal Spatial Transcriptomics**

GraphBG is a unified and scalable framework for spatial domain detection in spatial transcriptomics (ST) and spatial multi-omics data. It combines **spectral graph convolutions** with a **variational Bayesian Gaussian mixture model (VB-GMM)** to identify spatially coherent domains accurately and efficiently.

GraphBG supports three analysis settings:

- **GraphBG:** Spatial domain detection for unimodal ST data.
- **GraphBG-MS:** Multi-slice analysis using metacell aggregation, batch correction, and joint clustering.
- **GraphBG-MM:** Multimodal analysis using modality-specific graph embeddings and kernel canonical correlation analysis (KCCA).

## Key Features

- **Scalable:** Handles datasets containing hundreds of thousands of cells.
- **Versatile:** Supports unimodal, multi-slice, and multimodal ST data.
- **Accurate:** Achieves competitive, state-of-the-art performance across multiple benchmarks.
- **Efficient:** Processes datasets with more than 370,000 cells in approximately 5 minutes on standard hardware.
- **Open source:** Integrates with existing spatial transcriptomics analysis workflows.

## Installation

GraphBG requires **Python 3.8.20**. We recommend using Conda to create an isolated environment and manage dependencies.

### 1. Clone the repository

```bash
git clone git@github.com:CamiLQDTULab/GraphBG.git
cd GraphBG
```

### 2. Create and activate the Conda environment

```bash
conda create -n graphbg python=3.8.20 -y
conda activate graphbg
```

### 3. Install dependencies

Upgrade pip and install the packages listed in `requirements.txt`:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

**Dependency notes:**

- **Scanpy and Matplotlib:** If `scanpy==1.9.8` causes a Matplotlib-related import error, try installing Scanpy 1.9.8 without reinstalling its dependencies:

  ```bash
  python -m pip install --no-deps "scanpy==1.9.8"
  ```

  This assumes the required dependencies are already installed and compatible.

- **Mvlearn:** This package may introduce dependency conflicts with the Matplotlib version required by Scanpy. Install it only if it is needed for your workflow.

### 4. Verify the installation

Check dependency consistency and verify that Scanpy and POT can be imported:

```bash
python -m pip check
python -c "import scanpy as sc; print('Scanpy:', sc.__version__)"
python -c "import ot; print('POT:', ot.__version__)"
```

Note that `pip check` reports dependency conflicts but does not guarantee that every package can be imported successfully.

## Usage

### 1. Unimodal spatial transcriptomics

Run the Jupyter notebook:

```text
runGraphBG.ipynb
```

Alternatively, execute the Python script:

```bash
python runGraphBG.py
```

### 2. Multi-slice analysis (GraphBG-MS)

Run the notebook:

```text
runGraphBG-MS.ipynb
```

### 3. Multimodal analysis (GraphBG-MM)

Run the notebook:

```text
runGraphBG-MM.ipynb
```

**Note:** Run the notebooks in an environment where all required dependencies are installed. Update the input data paths and analysis parameters as needed for your dataset.

## Reproducibility

Scripts, notebooks, and supporting materials for reproducing the analyses are available in:

```text
reproducible_analysis/
```

Refer to the files in this directory for the corresponding datasets, configurations, and analysis procedures.

## Citation

If you use GraphBG in your research, please cite:

> Do, V. H., Tran, T. P. L., & Canzar, S.  
> *GraphBG: Fast Bayesian Domain Detection via Spectral Graph Convolutions for Multi-slice and Multimodal Spatial Transcriptomics.*  
> Preprint, 2025.

## License

GraphBG is released under the **BSD 3-Clause License**. See the `LICENSE` file for details.
