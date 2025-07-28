# Transcriptomic profiling of Plasmodium-induced lung pathology

This repository supports a study investigating the immune mechanisms underlying **Plasmodium-induced lung injury**, focusing on the role of **IFN-γ signaling in T cells** and its impact on **T cell–monocyte interactions** and pulmonary pathology.

- **Bulk RNA-seq** and **spatial transcriptomics** were performed on a murine *Plasmodium* lung infection model to characterize temporal and spatial dynamics of immune responses.
- **scRNA-seq data** from a MA-ARDS model were used as a complementary approach to validate immune cell subsets and regulatory interactions.

The datasets provided here enable comprehensive, multi-dimensional exploration of the host immune landscape during *Plasmodium* infection.

---

# Script guidance

The project includes raw count matrices and R scripts to analyze bulk RNA-seq, single-cell RNA-seq, and spatial transcriptomics data. We provide data from multiple time points, as well as pre- and post-infection conditions, to allow comprehensive exploration of immune responses in **Plasmodium-induced lung injury**.

---

## 1) Data

- **01_Bulk_DEGs.Rds, 01_Bulk_Matrix.Rds, 01_Bulk_SampleTable.Rds**: Bulk RNA-seq data matrices, including differential gene expression analysis and sample information.
- **01_biomart_hom_mus_gene.Rds**: Gene annotation data from Biomart, including homologous gene information for _Mus musculus_.
- **02_01_scedata.Rds, 02_01_scedata_Normaldata.Rds, 02_01_scedata_PCA.Rds, 02_01_scedata_UMAP.Rds**: Single-cell RNA-seq data, including normalized data and dimensionality reduction (PCA, UMAP) results.
- **03_01_Infect.Rds, 03_01_Infect_Ifngr1KO.Rds, 03_01_Uninfect.Rds**: Spatial transcriptomics data of infected, Ifngr1KO, and uninfected samples for comparative analysis.

---

## 2) Script

- **01_Bulk_RNAseq_01_SampleDistance.Rmd**: Analysis of sample distance in bulk RNA-seq data.
- **01_Bulk_RNAseq_02_DEG.Rmd**: Differential expression analysis using [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).
- **01_Bulk_RNAseq_03_Mfuzz.Rmd**: [Mfuzz](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2139991/) clustering analysis for temporal gene expression patterns.
- **01_Bulk_RNAseq_04_Xcell.Rmd**: Cell type enrichment analysis using [Xcell](https://xcell.ucsf.edu/).
- **02_ScRNAseq_01_Quantity.Rmd**: Quantification of immune cell types in scRNA-seq data.
- **02_ScRNAseq_02_DEG.Rmd**: scRNA-seq differential expression analysis.
- **02_ScRNAseq_03_NicheNet.Rmd**: Cell-cell communication analysis using [NicheNet](https://www.nature.com/articles/s41592-019-0667-5).
- **03_SpatialData_01_CCA.Rmd**: Batch Integration Using Canonical Correlation Analysis (CCA) for Spatial Transcriptomics Data.
- **03_SpatialData_02_SubRegion.Rmd**: Sub-region analysis of spatial transcriptomics data.
- **03_SpatialData_03_Deconvolution.Rmd**: Deconvolution of spatial transcriptomics data.

---
