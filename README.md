# scRNAseq_E_intestinalis

Single-cell RNA sequencing analysis of human macrophages infected with microsporidian *Encephalitozoon intestinalis*

**Table of contents**
- [Introduction](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/edit/main/README.md#introduction)
- [Workflow summary](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/tree/main?tab=readme-ov-file#workflow-summary)
- [Documentation](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/tree/main?tab=readme-ov-file#documentation)
- [Credits](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/tree/main?tab=readme-ov-file#credits)
- [Citation](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/tree/main?tab=readme-ov-file#citation)

## Introduction
In this study, we use single-cell RNA sequencing to investigate transcriptional changes in human macrophages when they are infected with microsporidian E. intestinalis. Primary human PBMCs were collected from 4 healthy donors, differentiated into macrophages, and infected with E. intestinalis for 3, 12, 24, 48, or 72 hours. Uninfected cells were used as negative control. 

Cell hashing was used to combine cells from different infection timepoints and the library preparation was prepared using Chromium Single-Cell 3' Reagent Kits v2 Chemistry (10x Genomics).

Three libraries were prepared and sequenced, including
1. Donor 1 library contains 10,100 cells (uninfected control, 24h, 48h, 72h)
2. Donor 2 library contains 12,279 cells (uninfected control, 3h, 12h, 24h)
3. Donor 3 and 4 library contains 30,355 cells (uninfected control, 24h, 48h)


## Workflow summary
### Processing of raw sequencing reads

- Generation of the combined reference genome between human and E. intestinalis
- Mapping raw reads with the combined reference genome and generation of gene expression matrix

### scRNA-seq data processing
- Demultiplexing cells from different infection timepoints
- Filtering of bad quality cells
- Normalization and dimensionality reduction
- Cell clustering
- Identification of marker genes and differentially expressed genes


## Documentation
Raw sequencing reads were mapped with the combined reference genome and the gene expression matrices were generated using Cell Ranger. Codes are available in the [Processing_of_raw_sequencing_reads.md](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/blob/main/Processing_of_raw_sequencing_reads.md) file.

After obtaining the gene expression matrices, they were processed using Seurat in R. We processed these datasets using three strategies
1. [Total transcriptome processing](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/blob/main/PJSC22_KM_Total_Transcriptome.R)
2. [Human only transcripts processing](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/blob/main/PJSC20_KM_Human_Only_Transcripts.R)
3. [Parasite only transcripts processing](https://github.com/pjaroenlak/scRNAseq_E_intestinalis/blob/main/PJSC23_KM_Parasite_Only_Transcripts.R)

R scripts that were used for data processing are uploaded in separate .R file for each data processing strategy

## Credits
- `Cell Ranger` from 10X Genomics was used to process raw sequencing reads
- `Seurat` in R was used to process scRNA-seq data

## Citation
Jaroenlak, P.*, McCarty, K.L.*, Xia, B., Lam, C., Zwack, E.E., Yanai, I., Bhabha, G. and Ekiert, D.C., 2024. scRNA-seq reveals transcriptional dynamics of Encephalitozoon intestinalis parasites in human macrophages. bioRxiv, pp.2024-05. (*equally contributed). *[Link.](https://www.biorxiv.org/content/10.1101/2024.05.30.596468v1.abstract "Link.")*****
