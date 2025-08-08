# Transcriptomic Profiling of *Plasmodium*-Induced Lung Pathology

This repository accompanies a study investigating the immune mechanisms underlying **Plasmodium-induced lung injury**, with a particular focus on the role of **IFN-γ signaling in T cells** and its impact on **T cell–monocyte interactions** and pulmonary pathology.

Key approaches:  
- **Bulk RNA-seq** and **spatial transcriptomics** performed on a murine *Plasmodium* lung infection model to map the temporal and spatial dynamics of immune responses.  
- **scRNA-seq** from a MA-ARDS model used to validate immune cell subsets and regulatory interactions.

The datasets provided enable a comprehensive, multi-dimensional exploration of the host immune landscape during *Plasmodium* infection.

---

## 📂 Script Guidance

This repository contains raw count matrices and R scripts for analyzing bulk RNA-seq, single-cell RNA-seq, and spatial transcriptomics data.  
Data include multiple time points, as well as pre- and post-infection conditions, allowing in-depth investigation of immune responses in **Plasmodium-induced lung injury**.

---

## 1️⃣ Data

- **01_Bulk_DEGs.Rds**, **01_Bulk_Matrix.Rds**, **01_Bulk_SampleTable.Rds** — Bulk RNA-seq data, including differential expression results and sample information.  
- **01_biomart_hom_mus_gene.Rds** — Gene annotation from Biomart, including homologous gene mappings for *Mus musculus*.  
- **02_01_scedata.Rds**, **02_01_scedata_Normaldata.Rds**, **02_01_scedata_PCA.Rds**, **02_01_scedata_UMAP.Rds** — Single-cell RNA-seq data with normalized counts and dimensionality reduction results (PCA, UMAP).  
- **03_01_Infect.Rds**, **03_01_Infect_Ifngr1KO.Rds**, **03_01_Uninfect.Rds** — Spatial transcriptomics data for infected, Ifngr1KO, and uninfected samples for comparative analysis.  

---

## 2️⃣ Scripts

**Bulk RNA-seq Analysis**  
- `01_Bulk_RNAseq_01_SampleDistance.Rmd` — Sample distance analysis.  
- `01_Bulk_RNAseq_02_DEG.Rmd` — Differential expression analysis ([DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)).  
- `01_Bulk_RNAseq_03_Mfuzz.Rmd` — Time-series clustering using [Mfuzz](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2139991/).  
- `01_Bulk_RNAseq_04_Xcell.Rmd` — Cell type enrichment with [Xcell](https://xcell.ucsf.edu/).  

**Single-cell RNA-seq Analysis**  
- `02_ScRNAseq_01_Quantity.Rmd` — Immune cell quantification.  
- `02_ScRNAseq_02_DEG.Rmd` — Differential expression analysis.  
- `02_ScRNAseq_03_NicheNet.Rmd` — Cell–cell communication analysis ([NicheNet](https://www.nature.com/articles/s41592-019-0667-5)).  

**Spatial Transcriptomics Analysis**  
- `03_SpatialData_01_CCA.Rmd` — Batch integration using Canonical Correlation Analysis (CCA).  
- `03_SpatialData_02_SubRegion.Rmd` — Sub-region analysis.  
- `03_SpatialData_03_Deconvolution.Rmd` — Cell type deconvolution.  

---

## 📄 Citation

If you use the **data** or **code** from this repository, please cite as:

> This repository accompanies a manuscript currently under peer review. Citation details will be updated upon publication.

