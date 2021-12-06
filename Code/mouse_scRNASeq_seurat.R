# Load necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

##load data
load('Mat_RNA_list_seurat') 
## see the time points
names(Mat_RNA_list_seurat)
## and get the individual seurat objects and sample into metadata
E11 = Mat_RNA_list_seurat$E11
E11$sample <- plyr::mapvalues(x = E11$orig.ident, from = c('RNA'), to = c('E11'))
E12 = Mat_RNA_list_seurat$E12
E12$sample = E12$orig.ident
E14 = Mat_RNA_list_seurat$E14
E14$sample = E14$orig.ident
E16 = Mat_RNA_list_seurat$E16
E16$sample <- plyr::mapvalues(x = E16$orig.ident, from = c('RNA'), to = c('E16'))
E18 = Mat_RNA_list_seurat$E18
E18$sample = E18$orig.ident
P0 = Mat_RNA_list_seurat$P0
P0$sample <- plyr::mapvalues(x = P0$orig.ident, from = c('RNA'), to = c('P0'))
P2 = Mat_RNA_list_seurat$P2
P2$sample = P2$orig.ident
P5 = Mat_RNA_list_seurat$P5
head(x = P5[[]])
P5$sample <- plyr::mapvalues(x = P5$orig.ident, from = c('RNA'), to = c('P5'))
P8 = Mat_RNA_list_seurat$P8
P8$sample = P8$orig.ident

# Merge into one single Seurat object
mouse=merge(E11, y=c(E12,E14,E16,E18,P0,P2,P5,P8))

##check labels
table(mouse$sample)

#Store mitochondrial percentage in the Seurat object metadata
mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^mt-")

#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1, group.by = 'sample')
dev.copy2pdf(file="mouse_scrnaseq_all_0621.qc.pdf", width=20)   

#Filter the data
mouse <- subset(mouse, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
mouse=NormalizeData(mouse)
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
mouse=ScaleData(mouse)
mouse = RunPCA(mouse)
DimPlot(mouse, reduction='pca', group.by='sample')
dev.copy2pdf(file="mouse_scrnaseq_all_0620.pca.pdf", width=20)
ElbowPlot(mouse, ndims=30)
dev.copy2pdf(file="mouse_scrnaseq_all_0620.elbow_plot.pdf", width=20)


#Run SCTransform, set var.to.regress to percent.mt --- ended up using harmony
mouse_merged <- SCTransform(mouse, vars.to.regress = "percent.mt", verbose = FALSE)
##resolution 0.3 and 30 dims
#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
mouse_merged <- RunPCA(mouse_merged, verbose = FALSE)
mouse_merged <- RunUMAP(mouse_merged, dims = 1:30, verbose = FALSE)
mouse_merged <- FindNeighbors(mouse_merged, dims = 1:30, verbose=FALSE)
mouse_merged <- FindClusters(mouse_merged, resolution = 0.3)
#Plot UMAP
DimPlot(mouse_merged, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="mouse_scrnaseq.merged.UMAP.res0.3.sample_split.pdf", width = 20)
DimPlot(mouse_merged, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="mouse_scrnaseq.merged.UMAP.res0.3.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(mouse_merged$seurat_clusters, mouse_merged$sample)
write.csv(counts_cluster_sample, file='mouse_scrnaseq.merged.counts_cluster_sample.res0.3.csv')

#Harmony needs RunPCA to have been performed, so make sure previous analysis steps have bene performed.
mouse_harmony <- RunHarmony(object = mouse, group.by.vars = 'sample', kmeans_init_nstart=20, kmeans_init_iter_max=100)
DimPlot(mouse_harmony, reduction = 'harmony', group.by = 'sample')
dev.copy2pdf(file="mouse_scrnaseq.harmony.harmony_plot.pdf", width = 20)

#Run standard Seurat workflow steps, but set reduction to "harmony" for UMAP and Neighbors
mouse_harmony <- RunUMAP(mouse_harmony, dims = 1:30, reduction = 'harmony')
mouse_harmony <- FindNeighbors(mouse_harmony, reduction = 'harmony', dims = 1:30)

##change resolutions - 0.3
mouse_harmony <- FindClusters(mouse_harmony, resolution = 0.3)
#Plot UMAP
DimPlot(mouse_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="mouse_scrnaseq.harmony.UMAP.res0.3.sample_split.pdf", width = 20)
DimPlot(mouse_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="mouse_scrnaseq.harmony.UMAP.res0.3.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(mouse_harmony$seurat_clusters, mouse_harmony$sample)
write.csv(counts_cluster_sample, file='mouse_scrnaseq.harmony.counts_cluster_sample.res0.3.csv')

#Find all markers that define each cluster
mouse_harmony.markers <- FindAllMarkers(mouse_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mouse_harmony.markers, file='mouse_scrnaseq.harmony.markers.res0.3.csv')

##save
saveRDS(mouse_harmony, file = "mouse_harmony.rds")



