#Seurat-based analysis pipeline of enhancer wt vs ko organoid snRNA-seq data

#R version 3.6.3

#Set working directory to active RSS to access files (this is specific to our workflow)
setwd('/active/cherry_t')

#Load necessary libraries
library(dplyr) #version 1.0.4
library(Seurat) #version 3.1.5
library(uwot) #version 0.1.9
library(patchwork) #version 1.1.0
library(ggplot2) #version 3.3.2

#Adjust maximum memory usage to fit memory allotted to session (200 GB)
options(future.globals.maxSize = 200 * 1024^3)

#Import 10x data and convert each dataset to Seurat object (take from individual sample outputs, not cellranger's aggr output). Define sample with project= "samplename"
#These file paths are specific to our environment; direct the Read10X function to a directory containing the barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz files
wt5_1.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/5_week_wt_vs_ko/Wt2-1/outs/filtered_feature_bc_matrix/')
wt5_1= CreateSeuratObject(counts = wt5_1.data, project = "wt5-1", min.cells = 3, min.features = 200)

wt5_2.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/5_week_wt_vs_ko/Wt3-1/outs/filtered_feature_bc_matrix/')
wt5_2= CreateSeuratObject(counts = wt5_2.data, project = "wt5-2", min.cells = 3, min.features = 200)

ko5_1.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/5_week_wt_vs_ko/Mut2-6/outs/filtered_feature_bc_matrix/')
ko5_1= CreateSeuratObject(counts = ko5_1.data, project = "ko5-1", min.cells = 3, min.features = 200)

ko5_2.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/5_week_wt_vs_ko/Mut3-8/outs/filtered_feature_bc_matrix/')
ko5_2= CreateSeuratObject(counts = ko5_2.data, project = "ko5-2", min.cells = 3, min.features = 200)

wt12_1.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/12_week_wt_vs_ko/WT1/outs/filtered_feature_bc_matrix/')
wt12_1= CreateSeuratObject(counts = wt12_1.data, project = "wt12-1", min.cells = 3, min.features = 200)

wt12_2.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/12_week_wt_vs_ko/WT2/outs/filtered_feature_bc_matrix/')
wt12_2= CreateSeuratObject(counts = wt12_2.data, project = "wt12-2", min.cells = 3, min.features = 200)

wt12_3.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/12_week_wt_vs_ko/WT3/outs/filtered_feature_bc_matrix/')
wt12_3= CreateSeuratObject(counts = wt12_3.data, project = "wt12-3", min.cells = 3, min.features = 200)

ko12_1.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/12_week_wt_vs_ko/KO1/outs/filtered_feature_bc_matrix/')
ko12_1= CreateSeuratObject(counts = ko12_1.data, project = "ko12-1", min.cells = 3, min.features = 200)

ko12_2.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/12_week_wt_vs_ko/KO2/outs/filtered_feature_bc_matrix/')
ko12_2= CreateSeuratObject(counts = ko12_2.data, project = "ko12-2", min.cells = 3, min.features = 200)

ko12_3.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/12_week_wt_vs_ko/KO3/outs/filtered_feature_bc_matrix/')
ko12_3= CreateSeuratObject(counts = ko12_3.data, project = "ko12-3", min.cells = 3, min.features = 200)

wt20_1.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/20_week_wt_vs_ko/WT1/outs/filtered_feature_bc_matrix/')
wt20_1= CreateSeuratObject(counts = wt20_1.data, project = "wt20-1", min.cells = 3, min.features = 200)

wt20_2.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/20_week_wt_vs_ko/WT2/outs/filtered_feature_bc_matrix/')
wt20_2= CreateSeuratObject(counts = wt20_2.data, project = "wt20-2", min.cells = 3, min.features = 200)

wt20_3.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/20_week_wt_vs_ko/WT3/outs/filtered_feature_bc_matrix/')
wt20_3= CreateSeuratObject(counts = wt20_3.data, project = "wt20-3", min.cells = 3, min.features = 200)

ko20_1.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/20_week_wt_vs_ko/KO1/outs/filtered_feature_bc_matrix/')
ko20_1= CreateSeuratObject(counts = ko20_1.data, project = "ko20-1", min.cells = 3, min.features = 200)

ko20_2.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/20_week_wt_vs_ko/KO2/outs/filtered_feature_bc_matrix/')
ko20_2= CreateSeuratObject(counts = ko20_2.data, project = "ko20-2", min.cells = 3, min.features = 200)

ko20_3.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/20_week_wt_vs_ko/KO3/outs/filtered_feature_bc_matrix/')
ko20_3= CreateSeuratObject(counts = ko20_3.data, project = "ko20-3", min.cells = 3, min.features = 200)

wt28_1.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/28_week_wt_vs_ko/WT1/outs/filtered_feature_bc_matrix/')
wt28_1= CreateSeuratObject(counts = wt28_1.data, project = "wt28-1", min.cells = 3, min.features = 200)

wt28_2.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/28_week_wt_vs_ko/WT2/outs/filtered_feature_bc_matrix/')
wt28_2= CreateSeuratObject(counts = wt28_2.data, project = "wt28-2", min.cells = 3, min.features = 200)

wt28_3.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/28_week_wt_vs_ko/WT3/outs/filtered_feature_bc_matrix/')
wt28_3= CreateSeuratObject(counts = wt28_3.data, project = "wt28-3", min.cells = 3, min.features = 200)

ko28_1.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/28_week_wt_vs_ko/KO1/outs/filtered_feature_bc_matrix/')
ko28_1= CreateSeuratObject(counts = ko28_1.data, project = "ko28-1", min.cells = 3, min.features = 200)

ko28_2.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/28_week_wt_vs_ko/KO2/outs/filtered_feature_bc_matrix/')
ko28_2= CreateSeuratObject(counts = ko28_2.data, project = "ko28-2", min.cells = 3, min.features = 200)

ko28_3.data = Read10X(data.dir = './OrgManuscript_SingleCell_Data/org_scRNA/28_week_wt_vs_ko/KO3/outs/filtered_feature_bc_matrix/')
ko28_3= CreateSeuratObject(counts = ko28_3.data, project = "ko28-3", min.cells = 3, min.features = 200)

# Merge into one single Seurat object
org = merge(wt5_1, y = c(wt5_2, ko5_1, ko5_2, wt12_1, wt12_2, wt12_3, ko12_1, ko12_2, ko12_3, wt20_1, wt20_2, wt20_3, ko20_1, ko20_2, ko20_3, wt28_1, wt28_2, wt28_3, ko28_1, ko28_2, ko28_3), add.cell.ids = c('wt5-1', 'wt5-2', 'ko5-1', 'ko5-2', 'wt12-1', 'wt12-2', 'wt12-3', 'ko12-1', 'ko12-2', 'ko12-3', 'wt20-1', 'wt20-2', 'wt20-3', 'ko20-1', 'ko20-2', 'ko20-3', 'wt28-1', 'wt28-2', 'wt28-3', 'ko28-1', 'ko28-2', 'ko28-3'))

#Store mitochondrial percentage in the Seurat object metadata
org[["percent.mt"]] <- PercentageFeatureSet(org, pattern = "^MT-")

#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(org, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

#Add sample and condition information to metadata

org$sample <- plyr::mapvalues(
  x = org$orig.ident, 
  from = c('wt5-1', 'wt5-2', 'ko5-1', 'ko5-2', 'wt12-1', 'wt12-2', 'wt12-3', 'ko12-1', 'ko12-2', 'ko12-3', 'wt20-1', 'wt20-2', 'wt20-3', 'ko20-1', 'ko20-2', 'ko20-3', 'wt28-1', 'wt28-2', 'wt28-3', 'ko28-1', 'ko28-2', 'ko28-3'), 
  to = c('wt5-1', 'wt5-2', 'ko5-1', 'ko5-2', 'wt12-1', 'wt12-2', 'wt12-3', 'ko12-1', 'ko12-2', 'ko12-3', 'wt20-1', 'wt20-2', 'wt20-3', 'ko20-1', 'ko20-2', 'ko20-3', 'wt28-1', 'wt28-2', 'wt28-3', 'ko28-1', 'ko28-2', 'ko28-3')
)

org$timepoint <- plyr::mapvalues(
  x = org$orig.ident, 
  from = c('wt5-1', 'wt5-2', 'ko5-1', 'ko5-2', 'wt12-1', 'wt12-2', 'wt12-3', 'ko12-1', 'ko12-2', 'ko12-3', 'wt20-1', 'wt20-2', 'wt20-3', 'ko20-1', 'ko20-2', 'ko20-3', 'wt28-1', 'wt28-2', 'wt28-3', 'ko28-1', 'ko28-2', 'ko28-3'), 
  to = c('wk5', 'wk5', 'wk5', 'wk5', 'wk12', 'wk12', 'wk12', 'wk12', 'wk12', 'wk12', 'wk20', 'wk20', 'wk20', 'wk20', 'wk20', 'wk20', 'wk28', 'wk28', 'wk28', 'wk28', 'wk28', 'wk28')
)

org$genotype <- plyr::mapvalues(
  x = org$orig.ident, 
  from = c('wt5-1', 'wt5-2', 'ko5-1', 'ko5-2', 'wt12-1', 'wt12-2', 'wt12-3', 'ko12-1', 'ko12-2', 'ko12-3', 'wt20-1', 'wt20-2', 'wt20-3', 'ko20-1', 'ko20-2', 'ko20-3', 'wt28-1', 'wt28-2', 'wt28-3', 'ko28-1', 'ko28-2', 'ko28-3'), 
  to = c('WT', 'WT', 'KO', 'KO', 'WT', 'WT', 'WT', 'KO', 'KO', 'KO', 'WT', 'WT', 'WT', 'KO', 'KO', 'KO', 'WT', 'WT', 'WT', 'KO', 'KO', 'KO')
)

#Subset the data based on QC metric cutoffs
org <- subset(org, subset = nFeature_RNA > 400 & nFeature_RNA < 5000 & percent.mt < 5)

#Save merged/subset Seurat object 
saveRDS(org, file = "/filepath/org.rds") #file path will obviously be user-dependent

#Run SCTransform and PCA, then check samples in PCA space
org <- SCTransform(org, vars.to.regress = "percent.mt", verbose = FALSE)
org <- RunPCA(org, verbose = FALSE)
DimPlot(org, reduction = 'pca')
DimPlot(org, reduction = 'pca', split.by = 'timepoint')
DimPlot(org, reduction = 'pca', split.by = 'genotype')
DimPlot(org, reduction = 'pca', split.by = 'sample', ncol = 5)
ElbowPlot(org, ndims=30)

#Samples align in PC space; no need to run Harmony. Generate/visualize UMAP
org <- RunUMAP(org, dims = 1:30, verbose = FALSE)
DimPlot(org, group.by = 'genotype')
DimPlot(org, group.by = 'timepoint')
DimPlot(org, group.by = 'genotype', split.by = 'genotype')
DimPlot(org, group.by = 'timepoint', split.by = 'timepoint')

#Finish standard Seurat workflow (FindNeighbors, FindClusters)
org <- FindNeighbors(org, dims = 1:30, verbose=FALSE)
org <- FindClusters(org, resolution = 0.4) 
DimPlot(org, reduction='umap', label=TRUE)+NoLegend()
DimPlot(org, reduction='umap', label=TRUE, split.by = 'genotype')+NoLegend()
DimPlot(org, reduction='umap', label=TRUE, split.by = 'timepoint')+NoLegend()

#Set default assay to RNA, then re-normalize (SCT assay is good for data visualization but not quantitative comparisons, according to Satija lab)
DefaultAssay(org) = 'RNA'
org = NormalizeData(org)

#Save Seurat object
saveRDS(org, file = "/filepath/org.rds")

#Re-examine QC metrics on a per-cluster basis
VlnPlot(org, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), pt.size=0)+NoLegend()

#Find all markers that define each cluster
markers <- FindAllMarkers(org, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file= '/filepath/cluster markers.csv')

#Examine number of cells per cluster per sample
table(org$seurat_clusters, org$sample)

#Subset out clusters 10 (low nCount_RNA spread AND enrichment of mitochondrial genes in marker list) and 18 (99% in one week12 sample - ko12-2)
Idents(org) = 'seurat_clusters'
org = subset(org, idents = c('10', '18'), invert=TRUE)

#Re-cluster cells after subsetting
DefaultAssay(org) = 'SCT' #UMAP is generated in SCT assay, so need to reset here otherwise clustering will fail
org <- FindNeighbors(org, dims = 1:30, verbose=FALSE)
org <- FindClusters(org, resolution = 0.4)
DimPlot(org, reduction='umap', label=TRUE)+NoLegend()

#Find all markers that define new clusters after re-clustering
DefaultAssay(org) = 'RNA'
markers <- FindAllMarkers(org, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file= '/filepath/cluster markers_recluster.csv')

#Examine expression of key retinal cell class marker genes
rods = c('PROM1', 'CRX', 'RCVRN', 'OTX2', 'RHO', 'NR2E3', 'NRL')
mg = c('SLC1A3', 'SLN', 'RLBP1', 'LINC00461', 'SOX2', 'NFIA', 'CRYM', 'CLU')
cones = c('ARR3', 'PROM1', 'CRX', 'GNAT2', 'OPN1SW', 'PDE6H')
hc = c('ONECUT1', 'ONECUT2', 'ONECUT3', 'TFAP2B', 'LHX1')
bc = c('VSX1', 'VSX2', 'OTX2', 'GRM6', 'PRKCA', 'LHX4', 'PROX1', 'PCP2', 'PCP4', 'TRPM1')
ac = c('SLC6A9', 'GAD1', 'SLC32A1', 'TFAP2B', 'GAD2', 'SLC18A3', 'LHX9', 'TFAP2A', 'TFAP2C')
rgc = c('POU4F2', 'RBPMS', 'NEFM', 'GAP43', 'POU4F1', 'ELAVL4')
rpe = c('BEST1', 'RPE65', 'TIMP3', 'TRPM1')
prog = c('VIM', 'SOX2', 'SFRP2', 'MKI67', 'UBE2C', 'FGF19', 'CCND1', 'ID3')
prbcpre = c('C8ORF46', 'RXRG', 'ATOH7', 'PRDM1', 'GADD45G', 'DCT', 'LHX4', 'OTX2')

DotPlot(org, features = rods)+RotatedAxis()
DotPlot(org, features = cones)+RotatedAxis()
DotPlot(org, features = mg)+RotatedAxis()
DotPlot(org, features = hc)+RotatedAxis()
DotPlot(org, features = bc)+RotatedAxis()
DotPlot(org, features = ac)+RotatedAxis()
DotPlot(org, features = rgc)+RotatedAxis()
DotPlot(org, features = rpe)+RotatedAxis()
DotPlot(org, features = prog)+RotatedAxis()
DotPlot(org, features = prbcpre)+RotatedAxis()

#Based on gene expression patterns across clusters, assign cell class info to clusters and add to metadata
org$celltype <- plyr::mapvalues(
  x = org$seurat_clusters, 
  from = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'), 
  to = c('Progenitors', 'Rods', 'Muller Glia', 'Cones', 'Neurogenic Progenitors', 'Progenitors', 'PC/BC Precursors', 'RGCs', 'Horizontal Cells', 'Amacrine Cells', 'Glial Precursors', 'Primary Progenitors', 'Bipolar Cells', 'RPE Precursors', 'RPE')
)

#Examine UMAP with cell class info
#Create list of new colors for UMAPS:
new.cols = c('Rods' = 'purple', 'Cones' = 'red', 'Muller Glia' = 'lightcoral', 'Bipolar Cells' = 'magenta', 'Horizontal Cells' = 'navy', 'Amacrine Cells' = 'darkgreen', 'RGCs' = 'gold3', 'Progenitors' = 'skyblue1', 'Primary Progenitors' = 'turquoise1', 'Neurogenic Progenitors' = 'green2', 'Glial Precursors' = 'khaki', 'RPE' = 'grey56', 'RPE Precursors' = 'gray81', 'PC/BC Precursors' = 'darkorange')

DimPlot(org, group.by = 'celltype', cols = new.cols, label=T)+NoLegend()
DimPlot(org, group.by = 'celltype', cols = new.cols,label=T, split.by = 'genotype')+NoLegend()
DimPlot(org, group.by = 'celltype', cols = new.cols,label=T, split.by = 'timepoint')+NoLegend()

#Calculate average expression of genes within/near enhancer TAD per sample
tad = c('LINC00461', 'TMEM161B', 'MEF2C', 'TMEM161B-AS1', 'MEF2C-AS2', 'RASA1', 'CCNH', 'AC091826.2')
Idents(org) = 'sample'
avg = AverageExpression(org, features = tad, assays = 'RNA', group.by = 'sample')
write.csv(avg$RNA, file = '/filepath/avg expression_global.csv')

#Examine gene expression as a feature plot
FeaturePlot(org, features = 'LINC00461', split.by = 'genotype')
FeaturePlot(org, features = 'TMEM161B-AS1', split.by = 'genotype')

#Subset Seurat object by cell class to examine gene expression in a cell class-specific manner
Idents(org) = 'celltype'
mg = subset(org, idents = 'Muller Glia')
pp = subset(org, idents = 'Primary Progenitors')
np = subset(org, idents = 'Neurogenic Progenitors')
prog = subset(org, idents = 'Progenitors')

#Calculate average expression of LINC00461 per sample in each cell class
Idents(mg) = 'sample'
Idents(pp) = 'sample'
Idents(np) = 'sample'
Idents(prog) = 'sample'
avg.mg = AverageExpression(mg, features = 'LINC00461', assays = 'RNA', group.by = 'sample')
avg.pp = AverageExpression(pp, features = 'LINC00461', assays = 'RNA', group.by = 'sample')
avg.np = AverageExpression(np, features = 'LINC00461', assays = 'RNA', group.by = 'sample')
avg.prog = AverageExpression(prog, features = 'LINC00461', assays = 'RNA', group.by = 'sample')

write.csv(avg.mg$RNA, file = '/filepath/avg linc expression_mullers.csv')
write.csv(avg.pp$RNA, file = '/filepath/avg linc expression_primary prog.csv')
write.csv(avg.np$RNA, file = '/filepath/avg linc expression_neurogenic prog.csv')
write.csv(avg.prog$RNA, file = '/filepath/avg linc expression_early prog.csv')

#Determine differentially expressed genes between wt and ko in Muller Glia and wk12 progenitors (Primary and Neurogenic)
Idents(mg) = 'genotype'
Idents(pp) = 'genotype'
Idents(np) = 'genotype'

de.mg = FindAllMarkers(mg, logfc.threshold = 0)
de.mg$log2FC = log2(exp(de.mg$avg_logFC)) #LogFC values are reported as natural logs, so convert to log2fc
write.csv(de.mg, file = '/filepath/de_mullers.csv')
de.pp = FindAllMarkers(pp, logfc.threshold = 0)
de.pp$log2FC = log2(exp(de.pp$avg_logFC))
write.csv(de.pp, file = '/filepath/de_primary prog.csv')
de.np = FindAllMarkers(np, logfc.threshold = 0)
de.np$log2FC = log2(exp(de.np$avg_logFC))
write.csv(de.np, file = '/filepath/de_neurogenic prog.csv')

