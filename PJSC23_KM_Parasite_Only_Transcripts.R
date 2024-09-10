#Set the working directory
setwd("/gpfs/data/bhabhaekiertlabs/home/km5313/scrnaseq/RStudio/PJSC23_KM/PJSC23_KM_3")

#Load necessary libraries
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

#Processing strategies: we will separate E. intestinalis transcripts from each donor.Then normalize the reads prior to integration of 4 donors

#Donor 1:
raw.data <- Read10X(data.dir="/gpfs/data/bhabhaekiertlabs/home/jaroep01/MichaelPJ/scRNA-seq/mapping_v3/Macrophages_Ei_v3/outs/filtered_feature_bc_matrix")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features, 10,100 cells
raw.HTO.count <- raw.data$`Antibody Capture` #4 features, 10,100 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs] 
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
row.names(raw.HTO.count) <- c('ctrl','1dpi','2dpi','3dpi')
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
percent.mt <- PercentageFeatureSet(gene.hashtag, pattern = "^MT-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.hashtag, pattern = "^Eint-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.microsporidia, col.name = "percent.microsporidia")

#separate E. intestinalis transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor1.microsporidian <- subset(gene.hashtag[36602:38612,])
rownames(Donor1.microsporidian)
Donor1.microsporidian <- NormalizeData(Donor1.microsporidian, normalization.method = "LogNormalize", scale.factor = 10000)
Donor1.microsporidian <- FindVariableFeatures(Donor1.microsporidian, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor1.microsporidian)
Donor1.microsporidian <- ScaleData(Donor1.microsporidian, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor1.microsporidian[["HTO"]] <- CreateAssayObject(counts=raw.HTO.count)
Donor1.microsporidian <- MULTIseqDemux(Donor1.microsporidian, assay = "HTO")
table(Donor1.microsporidian$MULTI_ID)
#1dpi     2dpi     3dpi     ctrl  Doublet Negative 
#1559     1631     2136     2870      878     1026 
Donor1.microsporidian <- subset(Donor1.microsporidian, idents = c("ctrl","1dpi","2dpi","3dpi"))
dim(Donor1.microsporidian) #8,196 cells

VlnPlot(Donor1.microsporidian, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.microsporidia"))

PJSC20_infected_cell_barcodes_and_clusters_D1 <- subset(PJSC20_infected_cell_barcodes_and_clusters_updated, subset = PJSC20_infected_cell_barcodes_and_clusters_updated$Exp == 'Donor1' & PJSC20_infected_cell_barcodes_and_clusters_updated$percent.microsporidia > 2)
dim(PJSC20_infected_cell_barcodes_and_clusters_D1)
#1295 cells

PJSC20_infected_cell_barcodes_and_clusters_D2 <- subset(PJSC20_infected_cell_barcodes_and_clusters_updated, subset = PJSC20_infected_cell_barcodes_and_clusters_updated$Exp == 'Donor2' & PJSC20_infected_cell_barcodes_and_clusters_updated$percent.microsporidia > 2)
dim(PJSC20_infected_cell_barcodes_and_clusters_D2)
#1108

PJSC20_infected_cell_barcodes_and_clusters_D3 <- subset(PJSC20_infected_cell_barcodes_and_clusters_updated, subset = PJSC20_infected_cell_barcodes_and_clusters_updated$Exp == 'Donor3' & PJSC20_infected_cell_barcodes_and_clusters_updated$percent.microsporidia > 2)
dim(PJSC20_infected_cell_barcodes_and_clusters_D3)
#1925

PJSC20_infected_cell_barcodes_and_clusters_D4 <- subset(PJSC20_infected_cell_barcodes_and_clusters_updated, subset = PJSC20_infected_cell_barcodes_and_clusters_updated$Exp == 'Donor4' & PJSC20_infected_cell_barcodes_and_clusters_updated$percent.microsporidia > 2)
dim(PJSC20_infected_cell_barcodes_and_clusters_D4)
#2580

Donor1.microsporidia_infected <- subset(Donor1.microsporidian@meta.data, subset = Donor1.microsporidian@meta.data$percent.microsporidia > 2)

cell.use <- PJSC20_infected_cell_barcodes_and_clusters_D1$...1
Donor1.microsporidian <- subset(Donor1.microsporidian, cells = cell.use)
dim(Donor1.microsporidian)
#1295 cells, 2011 features 

PJSC20_clusters_D1 <- PJSC20_infected_cell_barcodes_and_clusters_D1$integrated_snn_res.0.2
Donor1.microsporidian <- AddMetaData(Donor1.microsporidian, metadata = PJSC20_clusters_D1, col.name = "PJSC20_clusters")


remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)

#Donor 2
raw.data <- Read10X(data.dir = "/gpfs/data/bhabhaekiertlabs/home/jaroep01/MichaelPJ/scRNA-seq/Macrophage-Ei-rep2/processing/Macrophage-Ei-rep2/outs/filtered_feature_bc_matrix")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features and 12,279 cells
raw.HTO.count <- raw.data$`Antibody Capture`  #4 features and 12,279 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs]
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
row.names(raw.HTO.count) <- c("ctrl","3hpi","12hpi","1dpi")
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
percent.mt <- PercentageFeatureSet(gene.hashtag, pattern = "^MT-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.hashtag, pattern = "^Eint-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.microsporidia, col.name = "percent.microsporidia")

#separate E. intestinalis transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor2.microsporidian <- subset(gene.hashtag[36602:38612,])
rownames(Donor2.microsporidian)
Donor2.microsporidian <- NormalizeData(Donor2.microsporidian, normalization.method = "LogNormalize", scale.factor = 10000)
Donor2.microsporidian <- FindVariableFeatures(Donor2.microsporidian, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor2.microsporidian)
Donor2.microsporidian <- ScaleData(Donor2.microsporidian, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor2.microsporidian[["HTO"]] <- CreateAssayObject(counts = raw.HTO.count)
Donor2.microsporidian <- MULTIseqDemux(Donor2.microsporidian, assay = "HTO")
table(Donor2.microsporidian$MULTI_ID)
#12hpi     1dpi     3hpi     ctrl  Doublet Negative 
#3483     3219     3434      415      628     1100
remove(gene.hashtag)
Donor2.microsporidian <- subset(Donor2.microsporidian, idents = c("ctrl","3hpi","12hpi","1dpi"))
dim(Donor2.microsporidian) #10,551 cells

#Subset only infected cells that contain more than 2% of microsporidian transcripts from PJSC015
#By doing this, we don't have to perform filtering out
cell.use <- PJSC15_infected_cell_barcodes_and_clusters_D2$...1
Donor2.microsporidian <- subset(Donor2.microsporidian, cells = cell.use)
dim(Donor2.microsporidian)
#1,108 cells and 2011 features

#Add human clusters
view(Donor2.microsporidian)
PJSC20_clusters_D2 <- PJSC20_infected_cell_barcodes_and_clusters_D2$integrated_snn_res.0.2
Donor2.microsporidian <- AddMetaData(Donor2.microsporidian, metadata = PJSC20_clusters_D2, col.name = "PJSC20_clusters")

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)

#Donor3 and Donor4
raw.data <- Read10X(data.dir="/gpfs/data/bhabhaekiertlabs/home/jaroep01/MichaelPJ/scRNA-seq/GTC-PJ-6725/processing/GTC-PJ-6725/outs/filtered_feature_bc_matrix")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features and 30,355 cells
raw.HTO.count <- raw.data$`Antibody Capture`  #6 features and 30,355 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs]
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
rownames(raw.HTO.count) <- c("D3_ctrl","D3_1dpi","D3_2dpi","D4_ctrl","D4_1dpi","D4_2dpi")
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
percent.mt <- PercentageFeatureSet(gene.hashtag, pattern = "^MT-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.hashtag, pattern = "^Eint-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.microsporidia, col.name = "percent.microsporidia")

#separate E. intestinalis transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor3_4.microsporidian <- subset(gene.hashtag[36602:38612,])
rownames(Donor3_4.microsporidian)
Donor3_4.microsporidian <- NormalizeData(Donor3_4.microsporidian, normalization.method = "LogNormalize", scale.factor = 10000)
Donor3_4.microsporidian <- FindVariableFeatures(Donor3_4.microsporidian, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor3_4.microsporidian)
Donor3_4.microsporidian <- ScaleData(Donor3_4.microsporidian, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor3_4.microsporidian[["HTO"]] <- CreateAssayObject(counts = raw.HTO.count)
Donor3_4.microsporidian <- MULTIseqDemux(Donor3_4.microsporidian, assay = "HTO")
table(Donor3_4.microsporidian$MULTI_ID)
#D3-1dpi  D3-2dpi  D3-ctrl  D4-1dpi  D4-2dpi  D4-ctrl  Doublet Negative 
#3500     5984      830     5583     4431     1572     6535     1920 

Donor3.microsporidian <- subset(Donor3_4.microsporidian, idents = c("D3-ctrl","D3-1dpi","D3-2dpi"))
Donor4.microsporidian <- subset(Donor3_4.microsporidian, idents = c("D4-ctrl","D4-1dpi","D4-2dpi"))

#Subset Donor 3: subset only infected cells that contain more than 2% of microsporidian transcripts from PJSC015
#By doing this, we don't have to perform filtering out
cell.use <- PJSC15_infected_cell_barcodes_and_clusters_D3$...1
Donor3.microsporidian <- subset(Donor3.microsporidian, cells = cell.use)
dim(Donor3.microsporidian)
#1,925 cells and 2011 features

#Add human clusters
view(Donor3.microsporidian)
PJSC20_clusters_D3 <- PJSC20_infected_cell_barcodes_and_clusters_D3$integrated_snn_res.0.2
Donor3.microsporidian <- AddMetaData(Donor3.microsporidian, metadata = PJSC20_clusters_D3, col.name = "PJSC20_clusters")
view(Donor3.microsporidian)

#Subset Donor 4: subset only infected cells that contain more than 2% of microsporidian transcripts from PJSC015
#By doing this, we don't have to perform filtering out
cell.use <- PJSC15_infected_cell_barcodes_and_clusters_D4$...1
Donor4.microsporidian <- subset(Donor4.microsporidian, cells = cell.use)
dim(Donor4.microsporidian)
#2,580 cells and 2011 features

#Add human clusters
view(Donor4.microsporidian)
PJSC20_clusters_D4 <- PJSC20_infected_cell_barcodes_and_clusters_D4$integrated_snn_res.0.2
Donor4.microsporidian <- AddMetaData(Donor4.microsporidian, metadata = PJSC20_clusters_D4, col.name = "PJSC20_clusters")
view(Donor4.microsporidian)

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)
remove(Donor3_4.microsporidian)
remove(PJSC20_clusters_D1)
remove(PJSC20_clusters_D3)
remove(PJSC20_clusters_D4)
remove(cell.use)

#Assign donor name into each cell and store the data in 'Exp' slot
Donor1.microsporidian$Exp <- "Donor1"
Donor2.microsporidian$Exp <- "Donor2"
Donor3.microsporidian$Exp <- "Donor3"
Donor4.microsporidian$Exp <- "Donor4"

saveRDS(Donor1.microsporidian, "Donor1.microsporidian.rds")
saveRDS(Donor2.microsporidian, "Donor2.microsporidian.rds")
saveRDS(Donor3.microsporidian, "Donor3.microsporidian.rds")
saveRDS(Donor4.microsporidian, "Donor4.microsporidian.rds")

#2022-05-20
#Perform integration and batch correction
#First, select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(Donor1.microsporidian, Donor2.microsporidian, Donor3.microsporidian, Donor4.microsporidian), nfeatures = 3000)
options(max.print = 100000)
features

#Identify anchors which will be used for integration
combined.anchors <- FindIntegrationAnchors(object.list = list(Donor1.microsporidian, Donor2.microsporidian, Donor3.microsporidian, Donor4.microsporidian), anchor.features = features)

remove(PJSC15_infected_cell_barcodes_and_clusters_D1)
remove(PJSC15_infected_cell_barcodes_and_clusters_D2)
remove(PJSC15_infected_cell_barcodes_and_clusters_D3)
remove(PJSC15_infected_cell_barcodes_and_clusters_D4)

remove(combined.dataset)

#Combined datasets
combined.dataset <- IntegrateData(anchorset = combined.anchors)
dim(combined.dataset)
#6,908 cells with 1,990 features

#Next, let's see the integration results
DefaultAssay(combined.dataset) <- "integrated"
combined.dataset <- ScaleData(combined.dataset, verbose = FALSE)
combined.dataset <- RunPCA(combined.dataset)
##The following 3 features requested have zero variance (running reduction without them): Eint-010010, Eint-050010, Eint-070830
percent.pc <- combined.dataset@reductions$pca@stdev / sum(combined.dataset@reductions$pca@stdev) * 100
cum.percent.pc <- cumsum(percent.pc)
cum.percent.pc
#[1]  11.71819  18.42599  23.05443  26.48881  29.75891  32.29309  34.73331  36.74845  38.69800  40.43277  42.08944  43.71448  45.32395  46.91712
#[15]  48.47919  50.02434  51.54678  53.05754  54.56615  56.06987  57.56857  59.05778  60.54122  62.02402  63.50400  64.98250  66.45899  67.93477
#[29]  69.40718  70.87786  72.34730  73.81530  75.28040  76.74429  78.20676  79.66846  81.12952  82.59016  84.04760  85.50437  86.95985  88.41404
#[43]  89.86706  91.31896  92.76827  94.21630  95.66380  97.11060  98.55603 100.00000
#PC1-44

#Run UMAP
set.seed(2022)
combined.dataset <- RunUMAP(combined.dataset, dims = 1:44)
Idents(combined.dataset) <- "Exp"

png("PJSC23_UMAP_microsporidia_all_donors.png", width = 12, height = 8, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", pt.size = 0.1, cols = "Set1")
dev.off()

png("PJSC23_UMAP_microsporidia_all_donors_split.png", width = 40, height = 10, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", pt.size = 0.1, split.by = "Exp", ncol = 4, cols = "Set1")
dev.off()
#Donor 2 have only few cells on the island on the right

#Let's see how %microsporidia are distributed
png("PJSC23_UMAP_percent.microsporidia_NoLabel.png", width = 12, height = 8, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = "percent.microsporidia")
dev.off()
#The %microsporidia does not increase in 1 direction. The highest percentage is in the middle of the plot

#Split according to the infection timepoints
#First, we need to change the name of MULTI_ID for D3 and D4 from "D3-ctrl" to only "ctrl"
combined.dataset$MULTI_ID <- recode_factor(combined.dataset$MULTI_ID, "D3-ctrl" = "ctrl", "D3-1dpi" = "1dpi", "D3-2dpi" = "2dpi", "D4-ctrl"  = "ctrl", "D4-1dpi" = "1dpi", "D4-2dpi" = "2dpi")
#Then, order the timepoints
combined.dataset$MULTI_ID <- factor(combined.dataset$MULTI_ID, levels = c("ctrl","3hpi","12hpi","1dpi","2dpi","3dpi"))

png("PJSC23_UMAP_microsporidia_all_donors_split2.png", width = 32, height = 20, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", cols = "Set1", pt.size = 0.1) + facet_wrap(~combined.dataset$MULTI_ID + combined.dataset$Exp, ncol = 6)
dev.off()
#we can see how infection goes from early infection timepoint to late infection
#This shows asynchronistic infection of E. intestinalis

png("PJSC16_UMAP_microsporidia_all_donors_split3.png", width = 24, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", split.by = "MULTI_ID", ncol = 3, cols = "Set1", pt.size = 0.1) 
dev.off()

#Perform the clustering
set.seed(2022)
Idents(combined.dataset) <- "orig.ident"
combined.dataset <- FindNeighbors(combined.dataset, dims = 1:44)


#Resolution 0.35 --> DECIDED TO USE 0.35
combined.dataset <- FindClusters(combined.dataset, resolution = 0.35)
png("PJSC23_UMAP_res0.35_No_label.png", width = 12, height = 8, units = "cm", res = 300)
DimPlot(combined.dataset, pt.size = 0.1, label = T, label.size = 4)
dev.off()
#5 clusters


DimPlot(combined.dataset, split.by = "MULTI_ID", pt.size = 0.1) + facet_wrap(~combined.dataset$MULTI_ID + combined.dataset$Exp, ncol = 6)

png("PJSC23_Human_clusters_mapped_res0.35.png", width = 12, height = 8, units = "cm", res = 300)
DimPlot(combined.dataset, group.by = "PJSC20_clusters", pt.size = 0.4, label.size = 5)
dev.off()

png("PJSC23_human_clusters_split.png", width = 34, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, group.by = "PJSC20_clusters", pt.size = 0.4, label.size = 5, split.by = "MULTI_ID") 
dev.off()

png("PJSC23_UMAP_Split_Time_No_label_2.png", width = 25, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", pt.size = 0.1, ncol = 3)
dev.off()

######################################

#From finding all markers, we found that infection (by looking at the gene expression pattern) goes from C3 > C4 > C3 > C1 > C0
#Here we are going to change the order of the cluster to correlate with infection timepoints
#The old order = C3 > C1 > C0 > C2 > 
#The new order = C0 > C1 > C2 > C3 >

#save old cluster system (RNA_snn_res.0.35) as a reference
Idents(combined.dataset) <- "integrated_snn_res.0.35"
combined.dataset[["Old.cluster"]] <- Idents(combined.dataset)
view(combined.dataset)

#Name new cluster system
combined.dataset <- RenameIdents(object = combined.dataset, '3' = "0", '4' = "1", '2' ="2", '1' = "3", "0" = '4')
combined.dataset[["New.cluster"]] <- Idents(combined.dataset)
DimPlot(combined.dataset, label = T, pt.size = 0.4, label.size = 5)

#Export UMAP Plot
png("PJSC23_NEW_UMAP_Colors_res0.35.png", width = 12, height = 8, units = "cm", res = 300)
DimPlot(combined.dataset, pt.size = 0.1)
dev.off()

PJSC23_KM_cellbarcodes_and_clusters <- combined.dataset@meta.data[,c(5,15)]
write.csv(PJSC23_KM_cellbarcodes_and_clusters, "PJSC23_KM_cellbarcodes_and_clusters.csv")

write.csv(combined.dataset@meta.data, "PJSC23_MetaData.csv")
#####################################

#Let's see in each cluster where are the cells from
png("PJSC23_New_quantification.png", width = 10, height = 8, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.35, Exp) %>% count() %>% group_by(integrated_snn_res.0.35) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.35, y=percent, fill = Exp)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_brewer(palette = "Set1")
dev.off()

#Now let's see according to the infection timepoints
png("PJSC23_NEW_timepoint_in_each_cluster.png", width = 10, height = 8, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.35, MULTI_ID) %>% count() %>% group_by(integrated_snn_res.0.35) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.35, y=percent, fill = MULTI_ID)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_brewer(palette = "Set1")
dev.off() 


#Let's see where are the infected cells
png("PJSC23_NEW_UMAP_percent.microsporidia.png", width = 12, height = 8, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = "percent.microsporidia", label = T)
dev.off()

png("PJSC23_Timepoint_Colors.png", width = 15, height = 12, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", pt.size = 0.1) + facet_wrap(~combined.dataset$MULTI_ID, ncol = 3)
dev.off()

#Switch to RNA slot to perform DGE analysis
DefaultAssay(combined.dataset) <- "RNA"
RNA_allmarkers <- FindAllMarkers(combined.dataset, min.pct = 0.25)
write.csv(RNA_allmarkers, "PJSC23_0.35_RNA_allmarkers_New.csv")

RNA_allmarkers.top10 <- RNA_allmarkers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(RNA_allmarkers.top10, "PJSC23_RNA_allmarkers_top10.csv")

RNA_allmarkers.top30 <- RNA_allmarkers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
write.csv(RNA_allmarkers.top30, "PJSC23_RNA_allmarkers_top30.csv")

my_level <- c(0,1,2,3,4)
combined.dataset@active.ident <- factor(x = combined.dataset@active.ident, levels = my_level)

png("PJSC23_DotPlot_RNA_top10.png", width = 35, height = 8, units = "cm", res = 300)
DotPlot(combined.dataset, features = unique(RNA_allmarkers.top10$gene)) + RotatedAxis() 
dev.off()

RNA_allmarkers.top5 <- RNA_allmarkers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

png("PJSC23_FeaturePlot_RNA_top10.png", width = 35, height = 100, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = RNA_allmarkers.top5$gene, label = T)
dev.off()

png("PJSC123_FeaturePlot_RNA_top5.png", width = 40, height = 50, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = RNA_allmarkers.top5$gene, label = T, ncol = 5)
dev.off()

P1_P2_P3 <- FindMarkers(combined.dataset, ident.1 = "1", ident.2 = "2", indent.3 = "3")
write.csv(P1_P2_P3, "PJSC23_P1_P2_P3.csv")

#Make heatmap
DefaultAssay(combined.dataset) <- "integrated"
png("PJSC23_Heatmap_top10_0.35_FIGURE_NEW_Font.png", width = 1850, height = 1850, res = 300)
DoHeatmap(combined.dataset, features = unique(RNA_allmarkers.top10$gene), group.bar = T, group.colors = c("white","white","white","white","white"), draw.lines = T, lines.width = 50) + scale_fill_gradientn(colours = rev(brewer.pal(n=9, name="RdBu")), na.value = "white") + theme(axis.text.y = element_text(size = 10))
dev.off()

#Try DGE between clusters
#1) C4 vs C3
#2) C4 vs C2
#3) C4 vs C1
#4) C4 vs C0
#5) C3 vs C2
#6) C3 vs C1
#7) C3 vs C0
#8) C2 vs C1
#9) C2 vs C0
#10) C1 vs C0

C4_vs_C5 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "5")
C0_vs_C5 <- FindMarkers(combined.dataset, ident.1 = "0", ident.2 = "5")
C1_vs_C5 <- FindMarkers(combined.dataset, ident.1 = "1", ident.2 = "5")

C5_vs_C0 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "0")
C5_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "1")
C5_vs_C4 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "4")
C3_vs_C0 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "0")
C2_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "2", ident.2 = "1")
C2_vs_C0 <- FindMarkers(combined.dataset, ident.1 = "2", ident.2 = "0")
C1_vs_C0 <- FindMarkers(combined.dataset, ident.1 = "1", ident.2 = "0")

write.csv(C4_vs_C3, "PJSC15_C4_vs_C3.csv")
write.csv(C4_vs_C2, "PJSC15_C4_vs_C2.csv")
write.csv(C4_vs_C1, "PJSC15_C4_vs_C1.csv")
write.csv(C4_vs_C0, "PJSC15_C4_vs_C0.csv")
write.csv(C3_vs_C2, "PJSC15_C3_vs_C2.csv")
write.csv(C3_vs_C1, "PJSC15_C3_vs_C1.csv")
write.csv(C3_vs_C0, "PJSC15_C3_vs_C0.csv")
write.csv(C2_vs_C1, "PJSC15_C2_vs_C1.csv")
write.csv(C2_vs_C0, "PJSC15_C2_vs_C0.csv")
write.csv(C1_vs_C0, "PJSC15_C1_vs_C0.csv")

Idents(combined.dataset) <- "integrated_snn_res.0.35"
DimPlot(combined.dataset)

#######################
#2022-06-04
png("PJSC16_VlmPlot_RNA_top10.png", width = 35, height = 100, units = "cm", res = 300)
VlnPlot(combined.dataset, features = RNA_allmarkers.top10$gene, pt.size = 0)
dev.off()

#2022-06-05
#We observed 5 different expression patterns of microsporidian genes
#1) genes that expressed in all clusters      Eint_030430, Eint_030710, Eint_100850
#2) genes that expressed in C1, C2, C3, C4    Eint_110360, Eint_020735, Eint_091100
#3) genes that expressed in C2, C3, C4        Eint_041120, Eint_060710
#4) genes that expressed in C3, C4            Eint_017020, Eint_011200
#5) genes that expressed in C4                Eint_060150, Eint_081760

png("PJSC16_UMAP_microsporidia_all_donors_split4.png", width = 24, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", split.by = "MULTI_ID", ncol = 3, cols = "Set1", pt.size = 0.2) 
dev.off()
#At 12hpi, we can see that microsporidia develop to C0 - C1 - C2

Idents(combined.dataset) <- "integrated_snn_res.0.35"
png("PJSC16_UMAP_microsporidia_all_donors_split5.png", width = 24, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", split.by = "MULTI_ID", ncol = 3, pt.size = 0.3, label = T, label.size = 5) 
dev.off()

DefaultAssay(combined.dataset) <- "RNA"
png("PJSC16_UMAP_signatures.png", width = 60, height = 8, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = c("Eint-030710","Eint-110360","Eint-041120","Eint-010720","Eint-060150"), ncol = 5, label = T, label.size = 5)
dev.off()

#Make heatmap
DefaultAssay(combined.dataset) <- "integrated"
png("PJSC16_Heatmap_top10_2.png", width = 2000, height = 1900, res = 300)
DoHeatmap(combined.dataset, features = unique(RNA_allmarkers.top10$gene), group.bar = T, group.colors = c("white","white","white","white","white"), draw.lines = T, lines.width = 50) + scale_fill_gradientn(colours = rev(brewer.pal(n=9, name="RdBu")), na.value = "white")
dev.off()


##############################Now Trying 0.35
combined.dataset <- RenameIdents(object = combined.dataset, '3' = "0", '4' = "1", '3' = "2", '1' = "3", '0' = "4")
DimPlot(combined.dataset, label = T, pt.size = 0.4, label.size = 5)

png("PJSC23_UMAP_res0.35.png", width = 12, height = 8, units = "cm", res = 300)
DimPlot(combined.dataset, pt.size = 0.1, label = T, label.size = 4)
dev.off()

png("PJSC23_Res0.35_quantification.png", width = 10, height = 8, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.35, Exp) %>% count() %>% group_by(integrated_snn_res.0.35) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.35, y=percent, fill = Exp)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_brewer(palette = "Set1")
dev.off()

png("PJSC23_Res0.35_quantification_timepoint_in_each_cluster.png", width = 10, height = 8, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.35, MULTI_ID) %>% count() %>% group_by(integrated_snn_res.0.35) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.35, y=percent, fill = MULTI_ID)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_brewer(palette = "Set1")
dev.off()

png("PJSC23_UMAP_percent.microsporidia_0.35.png", width = 12, height = 8, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = "percent.microsporidia", label = T)
dev.off()

DefaultAssay(combined.dataset) <- "RNA"
RNA_allmarkers <- FindAllMarkers(combined.dataset, min.pct = 0.25)

RNA_allmarkers.top5 <- RNA_allmarkers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

my_level <- c(0,1,2,3,4)
combined.dataset@active.ident <- factor(x = combined.dataset@active.ident, levels = my_level)

png("PJSC23_DotPlot_RNA_top10_0.35.png", width = 35, height = 8, units = "cm", res = 300)
DotPlot(combined.dataset, features = unique(RNA_allmarkers.top5$gene)) + RotatedAxis() 
dev.off()

png("PJSC23_FeaturePlot_RNA_top10_0.35_New.png", width = 35, height = 100, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = RNA_allmarkers.top10$gene, label = T)
dev.off()

DefaultAssay(combined.dataset) <- "integrated"
png("PJSC23_Heatmap_top10_0.1.png", width = 2000, height = 2000, res = 300)
DoHeatmap(combined.dataset, features = unique(RNA_allmarkers.top10$gene), group.bar = T, group.colors = c("white","white","white","white","white"), draw.lines = T, lines.width = 5) + scale_fill_gradientn(colours = rev(brewer.pal(n=9, name="RdBu")), na.value = "white")
dev.off()

combined.dataset$PJSC20_clusters = as.factor(combined.dataset$PJSC20_clusters)
png("PJSC23_human_clusters_quantified.png", width = 12, height = 10, units = "cm", res = 300)
percent_infected_human<- combined.dataset@meta.data %>% group_by(New.cluster, PJSC20_clusters) %>% count() %>% group_by(New.cluster) %>% mutate(percent=100*n/sum(n))
write.csv(percent_infected_human, "human_clusters_quanitifed_table_New.csv")

png("Human_Clusters_Correct.png", width = 12, height = 10, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(New.cluster, PJSC20_clusters) %>% count() %>% group_by(New.cluster) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=New.cluster, y=percent, fill = PJSC20_clusters)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 17, face = "plain"))
dev.off()

DimPlot(combined.dataset, group.by = "PJSC20_clusters", pt.size = 0.4, label.size = 5)

png("Damians_Proteins.png", width = 12, height = 10, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = c("Eint-010900", "Eint-011330", "Eint-061420", "Eint-070250", "Eint-030060","Eint-031080", "Eint-051130", "Eint-060290", "Eint-070470", "Eint-111480"), label = T, pt.size = 0.4, label.size = 5, ncol = 3)
dev.off()

VlnPlot(combined.dataset, features = c("Eint-010900", "Eint-011330", "Eint-061420", "Eint-070250", "Eint-030060","Eint-031080", "Eint-051130", "Eint-060290", "Eint-070470", "Eint-111480"), pt.size = 0.4, ncol = 3)

library(ggplot2)
png("Ricin_Domain_ggplot_New_Blue.png", width = 10, height = 12, units = "cm", res = 300)
ggplot(Ricin_Prediction_updated, aes(x = Cluster, y = Gene, fill= FC)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="royalblue3") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1.0))
dev.off()


FeaturePlot(combined.dataset, features = c("Eint-020020", "Eint-091560", "Eint-050120", "Eint-070340"), label = T, pt.size = 0.4, label.size = 5, ncol = 2)

Average_expression <- AverageExpression(object = combined.dataset,assays="RNA",slot="data")
Average_expression
write.csv(Average_expression, "Average_Expression.csv")


DefaultAssay(combined.dataset) <- "RNA"
C0_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "0", ident.2 = "1")
C0_vs_C2 <- FindMarkers(combined.dataset, ident.1 = "0", ident.2 = "2")
C0_vs_C3 <- FindMarkers(combined.dataset, ident.1 = "0", ident.2 = "3")
C0_vs_C4 <- FindMarkers(combined.dataset, ident.1 = "0", ident.2 = "4")

C1_vs_C2 <- FindMarkers(combined.dataset, ident.1 = "1", ident.2 = "2")
C1_vs_C3 <- FindMarkers(combined.dataset, ident.1 = "1", ident.2 = "3")
C1_vs_C4 <- FindMarkers(combined.dataset, ident.1 = "1", ident.2 = "4")

C2_vs_C3 <- FindMarkers(combined.dataset, ident.1 = "2", ident.2 = "3")
C2_vs_C4 <- FindMarkers(combined.dataset, ident.1 = "2", ident.2 = "4")

C3_vs_C4 <- FindMarkers(combined.dataset, ident.1 = "2", ident.2 = "4")


write.csv(C0_vs_C1, "PJSC23_C0_vs_C1.csv")
write.csv(C0_vs_C2, "PJSC23_C0_vs_C2.csv")
write.csv(C0_vs_C3, "PJSC23_C0_vs_C3.csv")
write.csv(C0_vs_C4, "PJSC23_C0_vs_C4.csv")
write.csv(C1_vs_C2, "PJSC23_C1_vs_C2.csv")
write.csv(C1_vs_C3, "PJSC23_C1_vs_C3.csv")
write.csv(C1_vs_C4, "PJSC23_C1_vs_C4.csv")
write.csv(C2_vs_C3, "PJSC23_C2_vs_C3.csv")
write.csv(C2_vs_C4, "PJSC23_C2_vs_C4.csv")
write.csv(C3_vs_C4, "PJSC23_C3_vs_C4.csv")

png("Ei_Expression_Patterns.png", width = 40, height = 70, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = c("Eint-081980","Eint-071030","Eint-081300", "Eint-100480", "Eint-100460", "Eint-100370","Eint-111430","Eint-050300","Eint-020670", "Eint-080330", "Eint-031370", "Eint-090460", "Eint-090450","Eint-071150", "Eint-030710", "Eint-080620", "Eint-010330", "Eint-041070"), label = T, pt.size = 0.4, label.size = 5, ncol = 3)
dev.off()

png("Ei_Expression_Patterns_Vlnplot.png", width = 70, height = 100, units = "cm", res = 300)
VlnPlot(combined.dataset, features = c("Eint-081980","Eint-071030","Eint-081300", "Eint-100480", "Eint-100460", "Eint-100370","Eint-111430","Eint-050300","Eint-020670", "Eint-080330", "Eint-031370", "Eint-090460", "Eint-090450","Eint-071150", "Eint-030710", "Eint-080620", "Eint-010330", "Eint-041070"), pt.size = 0.4, ncol = 3)
dev.off()

RNA_Markers_log2FC_1 <- RNA_allmarkers[RNA_allmarkers$avg_log2FC > 1,]

RNA_Markers_log2FC_1.3 <- RNA_allmarkers[RNA_allmarkers$avg_log2FC > 1.3,]

RNA_Markers_log2FC_1.5 <- RNA_allmarkers[RNA_allmarkers$avg_log2FC > 1.5,]

RNA_Markers_log2FC_0.5 <- RNA_allmarkers[RNA_allmarkers$avg_log2FC > 0.58,]

write_csv(RNA_Markers_log2FC_1.3, "RNA_Markers_log2FC_1.3.csv")

DefaultAssay(combined.dataset) <- "integrated"
png("PJSC23_Heatmap_FC_2.8.png", width = 4000, height = 6000, res = 300)
DoHeatmap(combined.dataset, features = unique(RNA_Markers_log2FC_1.5$gene), group.bar = T, group.colors = c("white","white","white","white","white"), draw.lines = T, lines.width = 5) + scale_fill_gradientn(colours = rev(brewer.pal(n=9, name="RdBu")), na.value = "white") + theme(axis.text.y = element_text(size = 15))
dev.off()

##Went with FC greater than 2.5
png("PJSC23_Heatmap_FC_2.5_large_updated.png", width = 4000, height = 6500, res = 300)
DoHeatmap(combined.dataset, features = unique(RNA_Markers_log2FC_1.3$gene), group.bar = T, group.colors = c("white","white","white","white","white"), draw.lines = T, lines.width = 50) + scale_fill_gradientn(colours = rev(brewer.pal(n=9, name="RdBu")), na.value = "white") + theme(axis.text.y = element_text(size = 15))
dev.off()


DefaultAssay(combined.dataset) <- "integrated"
png("PJSC23_Heatmap_FC_2_Large_Flip.png", width = 10000, height = 5000, res = 300)
DoHeatmap(combined.dataset, features = unique(RNA_Markers_log2FC_1$gene), group.bar = F, group.colors = c("white","white","white","white","white"), draw.lines = T, lines.width = 50) + scale_fill_gradientn(colours = rev(brewer.pal(n=9, name="RdBu")), na.value = "white") + theme(axis.text.y = element_text(size = 25)) + RotatedAxis() + coord_flip()
dev.off()

DefaultAssay(combined.dataset) <- "integrated"
png("PJSC23_Heatmap_top10_0.1.png", width = 2000, height = 4000, res = 300)
DoHeatmap(combined.dataset, features = unique(RNA_Markers_log2FC_0.5$gene), group.bar = T, group.colors = c("white","white","white","white","white"), draw.lines = T, lines.width = 5) + scale_fill_gradientn(colours = rev(brewer.pal(n=9, name="RdBu")), na.value = "white")
dev.off()


library(ggplot2)
png("Top10_FC.png", width = 10, height = 12, units = "cm", res = 300)
ggplot(RNA_Markers_Top10_ggplot_FC, aes(cluster, gene, fill= FC)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="Red") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1.0))


png("Ei_P0-P2_Expression_Patterns_nolabel_New.png", width = 60, height = 10, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = c("Eint-010330", "Eint-030430", "Eint-100850", "Eint-010720", "Eint-060150"), pt.size = 0.4, label.size = 5, ncol = 5, cols = c("whitesmoke", "royalblue4"))
dev.off()

FeaturePlot(combined.dataset, features = c("Eint-010720"), pt.size = 0.4, label.size = 5, ncol = 6)

FeaturePlot(combined.dataset, features = c("Eint-051210",
                                           "Eint-111650",
                                           "Eint-080350",
                                           "Eint-041320",
                                           "Eint-060720",
                                           "Eint-070810",
                                           "Eint-090740",
                                           "Eint-111310",
                                           "Eint-101570",
                                           "Eint-041220",
                                           "Eint-061030",
                                           "Eint-040450",
                                           "Eint-111230",
                                           "Eint-060980",
                                           "Eint-061210",
                                           "Eint-040120",
                                           "Eint-061450",
                                           "Eint-060190",
                                           "Eint-041450",
                                           "Eint-101680",
                                           "Eint-010990"), label = T, pt.size = 0.4, label.size = 5)


FeaturePlot(combined.dataset, features = c("Eint-051210",
                                           "Eint-111650",
                                           "Eint-080350",
                                           "Eint-041320",
                                           "Eint-060720",
                                           "Eint-070810",
                                           "Eint-090740",
                                           "Eint-111310",
                                           "Eint-101570",
                                           "Eint-041220",
                                           "Eint-061030",
                                           "Eint-040450",
                                           "Eint-111230",
                                           "Eint-060980",
                                           "Eint-061210",
                                           "Eint-040120",
                                           "Eint-061450",
                                           "Eint-060190",
                                           "Eint-041450",
                                           "Eint-101680",
                                           "Eint-010990",
                                           "Eint-111810",
                                           "Eint-080300",
                                           "Eint-090380",
                                           "Eint-091720",
                                           "Eint-091610",
                                           "Eint-081310",
                                           "Eint-090530",
                                           "Eint-041580",
                                           "Eint-101580",
                                           "Eint-081640",
                                           "Eint-091620",
                                           "Eint-080720",
                                           "Eint-081500",
                                           "Eint-050150"), label = T, pt.size = 0.4, label.size = 5)

write.csv(combined.dataset@meta.data,"PJSC23_MetaData.csv")

FeaturePlot(combined.dataset, reduction = "umap", split.by = "MULTI_ID", features = c("Eint-060140"), ncol = 3)

combined.dataset@meta.data %>% group_by(New.cluster, percent.microsporidia) %>% count() %>% group_by(New.cluster) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=New.cluster, y=percent, fill = percent.microsporidia)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 17, face = "plain")) + scale_fill_brewer(palette = "Greys")
dev.off()

library(ggplot2)
png("Avg_Expression_ggplot_Flip_Font_deeppink_Lines.png", width = 55, height = 10, units = "cm", res = 300)
ggplot(PJSC23_Top_50_Avg_Expression_Log_2, aes(x = Gene, y = Cluster, fill= Log)) + 
  geom_tile(color = "white",
            lwd = 0.1,
            linetype = 1) +
  scale_fill_gradient(low="white", high="deeppink3") +
  RotatedAxis() +
  coord_fixed(2.0) +
  theme(axis.text.x = element_text(size = 12, color = "black", family = "helvetica")) +
  theme(axis.text.y = element_text(size = 12, color = "black", family = "helvetica")) +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1.0))
dev.off()

write.csv(combined.dataset@meta.data, "PJSC23_Meta.csv")

FeaturePlot(combined.dataset, features = c("Eint-051620"))

DefaultAssay(combined.dataset) <- "RNA"
png("Violin_Fig4D_NoDots.png", width = 35, height = 10, units = "cm", res = 300)
VlnPlot(combined.dataset, features = c("Eint-111310","Eint-030430","Eint-100850","Eint-010720","Eint-060150"), ncol = 5, pt.size = 0) & theme(axis.title.x = element_blank())
VlnPlot(combined.dataset, features = c("Eint-111310","Eint-030430","Eint-100850","Eint-010720","Eint-060150"), ncol = 5)
dev.off()

setwd("/gpfs/data/bhabhaekiertlabs/home/km5313/scrnaseq/RStudio/PJSC23_KM/PJSC23_KM_3/PJSC23_Counts")
png("Counts_Donor_1.png", width = 33, height = 10, units = "cm", res = 300)
VlnPlot(Donor1.microsporidian, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.microsporidia"), ncol = 4)
dev.off()

png("Counts_Donor_2.png", width = 33, height = 10, units = "cm", res = 300)
VlnPlot(Donor2.microsporidian, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.microsporidia"), ncol = 4)
dev.off()

png("Counts_Donor_3.png", width = 33, height = 10, units = "cm", res = 300)
VlnPlot(Donor3.microsporidian, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.microsporidia"), ncol = 4)
dev.off()

png("Counts_Donor_4.png", width = 33, height = 10, units = "cm", res = 300)
VlnPlot(Donor4.microsporidian, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.microsporidia"), ncol = 4)
dev.off()
