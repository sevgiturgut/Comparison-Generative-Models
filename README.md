# 🧬 Flow-Based and Tabular Generative Models for Single-cell RNA-seq Data

This repository contains the full implementation of the generative model experiments presented in our study, focusing on synthetic data generation for **single-cell RNA-seq** datasets using **flow-based models** and **tabular deep generative models**.

---

## 📁 Datasets

We evaluate models on the following single-cell datasets:

- **PBMC3K** – data downloaded from [Link](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)
- **PBMC68K** – processed data downloaded from the [ACTIVA repository](https://github.com/SindiLab/ACTIVA)
- **HCA-BM10K** – subset of Human Cell Atlas bone marrow data [Package](https://bioconductor.org/packages/release/data/experiment/html/HCAData.html)
- **Integrated Pancreatic Dataset** – Custom merged dataset ([GSE84133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133), [GSE85241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241), [E-MTAB-5061](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/), [GSE83139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83139]))

---

## Dataset Preprocessing and Tutorials

### 🧬 Integrated Pancreatic Dataset
- **Individual dataset processing scripts:**
  - `xin.R`
  - `seg.R`
  - `muraro.R`
  - `baron.R`
- **Merging script:**  
  - `sceProcess.R`
- **Batch effect removal visualization:**  
  - `BatchEffectRemoval.ipynb`

---

### 🩸 PBMC3K
- **Preprocessing script:** `PBMC3K.R`  
- **Tutorial:** [Seurat PBMC3K Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

---

### 🩸 PBMC68K
- **Note:** Preprocessed version was downloaded directly.

---

### 🧪 HCA-BM10K
- **Notebook:** `HCA-BM10K.ipynb`  
- **Tutorial:** [Bioconductor OSCA: HCA Human Bone Marrow 10x Genomics](https://bioconductor.org/books/3.12/OSCA/hca-human-bone-marrow-10x-genomics.html)


## 🔁 Data Splitting Strategies
- **Notebook:** `Main.ipynb` 
- **PBMC3K / PBMC68K**: Simple `train-test` split. Generative models trained **only on train**, synthetic samples generated to match **test set size**.
- **HCA-BM10K / Pancreatic**: `5-fold cross-validation` setup.
  - Synthetic samples are generated per fold.
  - If class labels (`cell_type`) are available, synthesis is done **per-class** up to the **Q3 count**.
  - **Only training data is used** in generation; validation/test are untouched.

---

## 🧪 Models

### 1. **Flow-based Models**
- **MAF-FB**: Baseline model using *Masked Affine Autoregressive Transform* (**Notebook:** `MAF_FB.ipynb`)
- **MOE-FB**: Our proposed model extending MAF-FB with: (**Notebook:** `MOE_FB.ipynb`)
  - **Learnable feature masking**
  - **Mixture-of-Experts attention** (10 heads, 4 experts)
  - **ActNorm** layer (optional via CLI)

Training params (shared):
- Epochs: `100`
- Batch size: `128`
- Hidden features: `1024`
- Learning rate: `1e-6`

### 2. **SDV Tabular Models**
- **CTGAN** – Conditional GAN for tabular data (**Notebook:** `SDV_CTGAN.ipynb`)
- **TVAE** – Variational Autoencoder for tabular data (**Notebook:** `SDV_TVAE.ipynb`)
- **GC** - Gaussian Copula for tabular data (**Notebook:** `SDV_GC.ipynb`)
---

## Analysis and Evaluation

### 🔗 Cell-Cell Interaction
- **Notebook:** `CellPhoneDB.ipynb`  
- **Includes:**  
  - Statistical analysis of cell–cell interactions  
  - Heatmaps  
  - Dot plots  

---

### 🧬 Differential Expression Analysis (DEG)
- **Notebook:** `DEG.ipynb`  
- **Criteria:**  
  - Adjusted *p*-value < 0.05  
  - log2 Fold Change > 1  
- **Notes:**  
  - Top 20 mutually exclusive genes were examined  

---

### 📊 Evaluation
- **Notebook:** `Evaluation.ipynb`  
- **Metrics:**  
  - **CD** (Correlation Discrepancy)  
  - **MMD** (Maximum Mean Discrepancy)  
  - **WD** (Wasserstein Distance)  
  - **Classification:** Precision, Recall, F1  


