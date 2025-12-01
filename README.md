# Single-nuclei RNA-sequencing reveals cellular contributions to progressive circuit dysfunction in Huntington’s Disease mice

## Abstract

Huntington’s disease (HD) involves progressive corticostriatal dysfunction; however, the timing of region-specific transcriptional changes remains unresolved at cellular resolution. Here, we provide a temporal single-nuclei transcriptomic atlas of the striatum and motor cortex in the zQ175 knock-in mouse model at early (6 months) and late (18 months) stages. Striatal projection neurons show extensive early dysregulation with progressive striosomal identity loss, whereas cortical pyramidal neurons exhibit layer-specific changes coinciding with motor symptom onset. By modeling genotype-dependent effects, we distinguished region- and cell type-specific gene signatures of core disease mechanisms from those of aging and compensatory adaptations. Integrating these findings with gene co-expression and transcription factor regulatory networks, we identified candidate regulators of compensatory disease signatures also noted in recently identified CAG repeat length-dependent genes from studies on HD patient brain tissues. Further, this resource reports cell type-specific gene expression programs with conserved dynamics in human HD brain and other mice datasets along with coordinated cross-regional transcriptional responses, establishing a comprehensive temporal framework to guide further work.

## Links to additional related repositories and protocols

[SPLiT-seq demultiplexing repository](https://github.com/DavidsonLabCHOP/SPLiT-seq)

[SPLiT-seq library preparation protocol (optimized for this study)](https://www.protocols.io/view/protocol-for-split-seq-4r3l27o8jg1y/v1)


## About repository

This repository includes all data processing and analysis code from Robbins et al. 2025 [link]. This README is structured to align with the manuscript organization and provides links to the relevant scripts. All necessary code is included to reproduce rsults from either raw sequencing reads or processed data objects. 

## Data Availability

All raw sequencing data generated as a part of this study has been deposited to the NIH NCBI Sequence Read Archive (SRA) and will be made available under the Bioproject accession number PRJNA1337572 upon publication. Processed data and results output will be made available on [Synapse](https://www.synapse.org) under the Project SynID: syn69978129. 

## SPLiT-seq Demultiplexing for pre-processing of raw sequencing data

To pre-process raw FASTQ files following Illumina sequencing, we developed a custom de-multiplexing pipeline that is compatible with the original chemistry and primers published in Rosenberg et al. 2018 
## Clustering and annotation of single-nuclei RNA-seq data (Figure 1 and S1)

After importing counts matrices into [**Seurat**](https://github.com/satijalab/seurat), we separately clustered nuclei from each region (motor cortex and striatum). Biological origin was identified by the SPLiT-seq round 1 barcode. Samples were integrated by batch using [**Harmony**](https://github.com/immunogenomics/harmony), representing two separate SPLiT-seq plates that were evenly balanced by experimental condition. Doublets were identified at a predicted rate of 7.5% and removed with [**scDblFinder**](https://github.com/plger/scDblFinder).

* [Importing counts into Seurat and demultiplexing by metadata](https://github.com/ashley-brooke/Robbins_HDsnRNAseq_2025/blob/main/Figure%201/01_DataImport_Seurat.Rmd)

#### Mouse Motor Cortex

* [Clustering and annotation of the motor cortex dataset](https://github.com/ashley-brooke/Robbins_HDsnRNAseq_2025/blob/main/Figure%201/02_QC_Preprocessing_Seurat_Ctx.Rmd)

  Nuclei were annotated using [**Azimuth**](https://github.com/satijalab/azimuth) and the mouse motor cortex reference dataset, [Yao et al., Nature 2021]  (https://www.nature.com/articles/s41586-021-03500-8). 

#### Mouse Striatum

* [Clustering and annotation of the striatum dataset](https://github.com/ashley-brooke/Robbins_HDsnRNAseq_2025/blob/main/Figure%201/02_QC_Preprocessing_Seurat_Str.Rmd)

  Nuclei were annotated by canonical marker expression and using [MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells), available from   the Allen Brain Atlas, with the [Mammalian Basal Ganglia Consensus Cell Type Atlas](https://alleninstitute.github.io/HMBA_BasalGanglia_Consensus_Taxonomy/).

## snRNA-seq Differential Gene Expression (Figure 1, 2, and S2)

First, we sought to identify cell type-specific gene dysregulation at both timepoints. Our main approach was to create "pseudocells" by aggregating UMI matrices for every 35 cells within each individual and cell type prior to differential gene expression (DGE) analysis using [**limma-voom**](https://bioconductor.org/packages/devel/bioc/html/limma.html) with sample weights. First described in [Gazestani and Kamath et al.,Cell 2023](https://www.sciencedirect.com/science/article/pii/S0092867423009790/pdfft?md5=6a5f57356e90b95afd827b9589de44dc&pid=1-s2.0-S0092867423009790-main.pdf), this pseudocell method addresses data sparsity, as well as pseudoreplication by accounting for individual as a random effect.

* [Cell type-specific Pseudocell DGE analysis](https://github.com/ashley-brooke/Robbins_HDsnRNAseq_2025/blob/main/DGE%20Analysis/05_DEGAnalysis_Pseudocell.Rmd)

For comparison, we also report our results using "pseudobulks", aggregating raw UMI counts for all cells within a cell type and individual. We report the correlations between both methods for log<sub>2</sub>(fold change) and *p*-value.  

* [Cell type-specific Pseudobulk DGE analysis](https://github.com/ashley-brooke/Robbins_HDsnRNAseq_2025/blob/main/DGE%20Analysis/06_DEGAnalysis_Pseudobulk.Rmd)
* [Comparison between DGE methods](https://github.com/ashley-brooke/Robbins_HDsnRNAseq_2025/blob/main/DGE%20Analysis/DEGAnalysis_methodcomparison.R)

DGE analysis code was adapted from the [original source code](https://app.terra.bio/#workspaces/fbrihuman/sconline_integrative_analysis/analysis/launch/Pseudocell_differential_expression_analysis.ipynb) provided in Gazestani and Kamath et al., Cell 2023.

## Assessing transcriptional vulnerabilities of Striatal Medium Spiny Neuron (MSN) subtypes (Figure 3 and S3)

To characterize MSN subtypes, stratified by compartment and projection, we subsetted and re-clustered MSN clusters (D1 and D2 MSNs) from the striatum dataset. 

## MSN-specific gene co-expression network dynamics in HD (Figure 4)

## Transcription Factor Expression and Regulon Activity in HD MSNs (Figure 5)
## HD disease signatures in layer-specific cortical neurons (Figure 6)

