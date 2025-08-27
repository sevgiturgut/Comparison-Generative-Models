# 🧬 Flow-Based and Tabular Generative Models for Single-cell RNA-seq Data

This repository contains the full implementation of the generative model experiments presented in our study, focusing on synthetic data generation for **single-cell RNA-seq** datasets using **flow-based models** and **tabular deep generative models**.

---

## 📁 Datasets

We evaluate models on the following single-cell datasets:

- **PBMC3K** – data downloaded from [Link](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) [Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
- **PBMC68K** – processed data downloaded from the [ACTIVA repository](https://github.com/SindiLab/ACTIVA)
- **HCA-BM10K** – subset of Human Cell Atlas bone marrow data (Tutorial)[https://bioconductor.org/books/3.12/OSCA/hca-human-bone-marrow-10x-genomics.html]
- **Integrated Pancreatic Dataset** – Custom merged dataset

---

## 🧪 Models

### 1. **Flow-based Models**
- **MAF-FB**: Baseline model using *Masked Affine Autoregressive Transform*
- **MOE-FB**: Our proposed model extending MAF-FB with:
  - **Learnable feature masking**
  - **Mixture-of-Experts attention** (10 heads, 4 experts)
  - **ActNorm** layer (optional via CLI)

Training params (shared):
- Epochs: `100`
- Batch size: `128`
- Hidden features: `1024`
- Learning rate: `1e-6`

### 2. **SDV Tabular Models**
- **CTGAN** – Conditional GAN for tabular data
- **TVAE** – Variational Autoencoder for tabular data
- **GC** - Gaussian Copula for tabular data
---

## 🔁 Data Splitting Strategies

- **PBMC3K / PBMC68K**: Simple `train-test` split. Generative models trained **only on train**, synthetic samples generated to match **test set size**.
- **HCA-BM10K / Pancreatic**: `5-fold cross-validation` setup.
  - Synthetic samples are generated per fold.
  - If class labels (`cell_type`) are available, synthesis is done **per-class** up to the **Q3 count**.
  - **Only training data is used** in generation; validation/test are untouched.

---


