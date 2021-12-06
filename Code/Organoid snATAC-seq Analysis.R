#ArchR Analysis of organoid snATAC-seq data

#R version 3.6.3

#Load necessary libraries
library(ArchR) #version 0.9.5
library(pheatmap) #version 1.0.1
library(ggplot2) #version 3.3.2
library(GenomicRanges) #version 1.38.0
library(chromVARmotifs) #version 0.2.0
library(ComplexHeatmap) #version 2.2.0
library(circlize) #version 0.4.1
library(harmony) #version 1.0

#Adjust maximum memory usage to fit memory allotted to session (200 GB)
options(future.globals.maxSize = 200 * 1024 ^ 3)

#Set default number of threads at 0.75x the number of cores used on the session (8 cores in this case)
addArchRThreads(threads = 6) 
#Set genome
addArchRGenome("hg38")

#To create arrow file, make character vector of all fragment files
# These file paths are specific to our environment; set file path to fragments file from cell ranger atac output

#Split into 2 groups, as only the first set of samples have an R of 0.9 or higher, and are thus good for doublet scoring

input1 = c('12wk1' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC2/wk12_wt_vs_ko/WT1/outs/fragments.tsv.gz",
           '12wk3' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC2/wk12_wt_vs_ko/WT3/outs/fragments.tsv.gz",
           '20wk1' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC/20wk/outs/fragments.tsv.gz",
           '20wk2' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC/20wk_c5_1/outs/fragments.tsv.gz",
           '28wk1' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC/28-1/outs/fragments.tsv.gz",
           '28wk2' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC/28-2/outs/fragments.tsv.gz")

input2 = c('5wk1' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC/5wk/outs/fragments.tsv.gz",
           '5wk2' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC/5wk_c5_1/outs/fragments.tsv.gz",
           '12wk2' = "/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC2/wk12_wt_vs_ko/WT2/outs/fragments.tsv.gz")

ArrowFiles1 <- createArrowFiles(
  inputFiles = input1,
  sampleNames = names(input1),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles2 <- createArrowFiles(
  inputFiles = input2,
  sampleNames = names(input2),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#Add doublet scores only to input1 (the samples that are good for doublet scoring)
doubScores <- addDoubletScores(
  input = ArrowFiles1,
  k = 10, 
  knnMethod = "UMAP",
  LSIMethod = 1
)

#Combine into one set of arrow files
ArrowFiles3 = c(ArrowFiles1, ArrowFiles2)

#Create an ArchR Project comprising all arrow files; set an output directory to which all files will be saved
org <- ArchRProject(
  ArrowFiles = ArrowFiles3, 
  outputDirectory = "/outputdirectory/",
  copyArrows = TRUE #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
)

#Plot QC metrics - log10(unique fragments) vs tss enrichment
df <- getCellColData(org, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = org, addDOC = FALSE)

#Make ridge plots for each sample for TSS enrichment scores
p1 = plotGroups(
  ArchRProj = org, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

#Make violin plots for each sample for TSS enrichment scores
p2 = plotGroups(
  ArchRProj = org, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

#Make ridge plots per sample for log10(unique nuclear fragments)
p3 = plotGroups(
  ArchRProj = org, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

#Make violin plots per sample for log10(unique nuclear fragments)
p4 = plotGroups(
  ArchRProj = org, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = org, addDOC = FALSE, width = 4, height = 4)


#Plot fragment size distributions and TSS enrichment per sample
p1 = plotFragmentSizes(ArchRProj = org)
p2 = plotTSSEnrichment(ArchRProj = org)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = org, addDOC = FALSE, width = 5, height = 5)

#Add timepoint info to ArchR object
cM <- confusionMatrix(org$Sample, org$cellNames)
labelOld <- rownames(cM)
labelOld
labelNew <- c('wk20', 'wk28', 'wk28', 'wk20', 'wk12', 'wk12', 'wk5', 'wk5', 'wk12')
org$Timepoint <- mapLabels(org$Sample, newLabels = labelNew, oldLabels = labelOld)

#Save ArchR Project
saveArchRProject(ArchRProj = org, outputDirectory = "/outputdirectory/", load = TRUE)

#Save counts per sample as .csv
write.csv(table(org$Sample), file = '/outputdirectory/counts per sample pre-filter.csv')

#Filter doublets
org = filterDoublets(org)

#Save counts per sample as .csv
write.csv(table(org$Sample), file = '/outputdirectory/counts per sample post-doublet-filter.csv')

#Run dimensionality reduction and harmony batch correction
org <- addIterativeLSI(
  ArchRProj = org,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1, 0.2, 0.4), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

library(harmony)
org <- addHarmony(
  ArchRProj = org,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

#Determine clusters
org <- addClusters(
  input = org,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.4
)
org <- addClusters(
  input = org,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "ClustersHarmony",
  resolution = 0.4
)

#Tabulate and plot number of cells per cluster
write.csv(table(org$Clusters, org$Sample), file='/outputdirectory/cells per cluster per sample.csv')
write.csv(table(org$ClustersHarmony, org$Sample), file='/outputdirectory/cells per cluster per sample harmony.csv')

#Visualize cells per cluster as heatmap
cM <- confusionMatrix(paste0(org$Clusters), paste0(org$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), 
                        border_color = "black")
plotPDF(plotList = p, name = "cluster_sample_heatmap.pdf", 
        ArchRProj = org, addDOC = FALSE, width = 5, height = 5)
cMh <- confusionMatrix(paste0(org$ClustersHarmony), paste0(org$Sample))
cMh <- cMh / Matrix::rowSums(cMh)
p <- pheatmap::pheatmap(mat = as.matrix(cMh), color = paletteContinuous("whiteBlue"), 
                        border_color = "black")
plotPDF(plotList = p, name = "cluster_sample_heatmap_harmony.pdf", 
        ArchRProj = org, addDOC = FALSE, width = 5, height = 5)

#Generate UMAPs
org <- addUMAP(
  ArchRProj = org, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

org <- addUMAP(
  ArchRProj = org, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 = plotEmbedding(ArchRProj = org, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 = plotEmbedding(ArchRProj = org, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 = plotEmbedding(ArchRProj = org, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 = plotEmbedding(ArchRProj = org, colorBy = "cellColData", name = "ClustersHarmony", embedding = "UMAPHarmony")
p5 = plotEmbedding(ArchRProj = org, colorBy = "cellColData", name = "Timepoint", embedding = "UMAPHarmony")
p6 = plotEmbedding(ArchRProj = org, colorBy = "cellColData", name = "ClustersHarmony", embedding = "UMAPHarmony")
plotPDF(p1,p2, name = "organoid UMAPs.pdf", ArchRProj = org, addDOC = FALSE, width = 5, height = 5)
plotPDF(p3,p4, name = "organoid UMAPs_harmony.pdf", ArchRProj = org, addDOC = FALSE, width = 5, height = 5)
plotPDF(p5,p6, name = "organoid UMAPs_harmony by timepoint.pdf", ArchRProj = org, addDOC = FALSE, width = 5, height = 5)

#Move forward with the harmony clusters/UMAP from now on
#Save ArchR Project
saveArchRProject(ArchRProj = org, outputDirectory = "/outputdirectory/", load = TRUE)

#Identify marker genes
markers <- getMarkerFeatures(
  ArchRProj = org, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "ClustersHarmony",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markers, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.csv(markerList, file = '/outputdirectory/marker_list.csv')


#Impute and visualize marker genes
org <- addImputeWeights(org)
cell.markers = read.csv(file = '/filepath/Table S2.csv') # read in marker genes from Table S2
rods = cell.markers$Rods
mullers = cell.markers$Muller.Glia
mullers = na.omit(mullers)
cones = cell.markers$Cones
cones = na.omit(cones)
hc = cell.markers$Horizontal.Cells
hc = na.omit(hc)
bc = cell.markers$Bipolar.Cells
bc = na.omit(bc)
ac = cell.markers$Amacrine.Cells
ac = na.omit(ac)
rgc = cell.markers$RGCs
rgc = na.omit(rgc)
rpe = cell.markers$RPE
rpe = na.omit(rpe)
rpc = cell.markers$Early.RPC
rpc = na.omit(rpc)
prbcpre = cell.markers$PR.BC.Precursors
prbcpre = na.omit(prbcpre)
eprpc = cell.markers$Early.Primary.RPC
eprpc = na.omit(eprpc)
lprpc = cell.markers$Late.Primary.RPC
lprpc = na.omit(lprpc)
enrpc = cell.markers$Early.Neurogenic.RPC
enrpc = na.omit(enrpc)
lnrpc = cell.markers$Late.Neurogenic.RPC
lnrpc = na.omit(lnrpc)

#Plot gene scores in UMAP space for celltype-specific markers
p.rods <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = rods, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.rods, 
        name = "Plot-UMAP-Rod-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.cones <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = cones, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.cones, 
        name = "Plot-UMAP-Cone-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.mullers <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = mullers, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.mullers, 
        name = "Plot-UMAP-MullerGlia-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.hc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = hc, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.hc, 
        name = "Plot-UMAP-Horizontal-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.ac <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = ac, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.ac, 
        name = "Plot-UMAP-Amacrine-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.bc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = bc, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.bc, 
        name = "Plot-UMAP-Bipolar-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.rgc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = rgc, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.rgc, 
        name = "Plot-UMAP-RGC-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.rpe <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = rpe, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.rpe, 
        name = "Plot-UMAP-RPE-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.rpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = rpc, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.rpc, 
        name = "Plot-UMAP-Early-RPC-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.eprpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = eprpc, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.eprpc, 
        name = "Plot-UMAP-Early-Primary-RPC-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.lprpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = lprpc, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.lprpc, 
        name = "Plot-UMAP-Late-Primary-RPC-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.enrpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = enrpc, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.enrpc, 
        name = "Plot-UMAP-Early-Neurogenic-RPC-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.lnrpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = lnrpc, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.lnrpc, 
        name = "Plot-UMAP-Late-Neurogenic-RPC-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.prbcpre <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneScoreMatrix", 
  name = prbcpre, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.prbcpre, 
        name = "Plot-UMAP-PR-BC-Precursor-Genes-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)


#Integrate scRNAseq data from Seurat analysis
seRNA =  readRDS('/filepath/org.wt.rds')
seRNA

#Unconstrained integration
org <- addGeneIntegrationMatrix(
  ArchRProj = org, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "celltype",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
pal <- paletteDiscrete(values = seRNA$celltype)
p1 <- plotEmbedding(org, colorBy = "cellColData", 
                    name = "predictedGroup_Un", embedding = 'UMAPHarmony', pal = pal)
plotPDF(p1, name = "UMAP_clusters_scRNA overlap.pdf", ArchRProj = org, addDOC = FALSE, width = 5, height = 5)

#Overlap looks fairly good, go ahead and re-run integration and add RNA info to arrow files
#Set threads to 1 for this part
addArchRThreads(threads = 1)
org <- addGeneIntegrationMatrix(
  ArchRProj = org, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seRNA,
  addToArrow = TRUE,
  force = TRUE,
  groupRNA = "celltype",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
addArchRThreads(threads = 12)
saveArchRProject(ArchRProj = org, outputDirectory = "/outputdirectory/", load = TRUE)

#Plot gene integration in UMAP space for cell class-specific markers for harmony dataset
p.rods <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = rods,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.rods, 
        name = "Plot-UMAP-Rod-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.cones <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = cones,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.cones, 
        name = "Plot-UMAP-Cone-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.mullers <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = mullers,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.mullers, 
        name = "Plot-UMAP-MullerGlia-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.hc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = hc,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.hc, 
        name = "Plot-UMAP-Horizontal-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.ac <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = ac,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.ac, 
        name = "Plot-UMAP-Amacrine-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.bc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = bc,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.bc, 
        name = "Plot-UMAP-Bipolar-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.rgc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = rgc,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.rgc, 
        name = "Plot-UMAP-RGC-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.rpe <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = rpe,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.rpe, 
        name = "Plot-UMAP-RPE-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)
p.rpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = rpc,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.rpc, 
        name = "Plot-UMAP-Early-RPC-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.eprpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = eprpc,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.eprpc, 
        name = "Plot-UMAP-Early-Primary-RPC-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.lprpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = lprpc,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.lprpc, 
        name = "Plot-UMAP-Late-Primary-RPC-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.enrpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = enrpc,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.enrpc, 
        name = "Plot-UMAP-Early-Neurogenic-RPC-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.lnrpc <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = lnrpc,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.lnrpc, 
        name = "Plot-UMAP-Late-Neurogenic-RPC-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

p.prbcpre <- plotEmbedding(
  ArchRProj = org, 
  colorBy = "GeneIntegrationMatrix", 
  name = prbcpre,
  log2Norm = T,
  continuousSet = "solarExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(org)
)
plotPDF(plotList = p.prbcpre, 
        name = "Plot-UMAP-PR-BC-Precursor-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

#Add scRNAseq cluster info to ArchR object
cM <- confusionMatrix(org$ClustersHarmony, org$predictedGroup_Un)
cM
labelOld <- rownames(cM)
labelOld
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
org$CellTypes_RNA <- mapLabels(org$ClustersHarmony, newLabels = labelNew, oldLabels = labelOld)
p1 <- plotEmbedding(org, colorBy = "cellColData", name = "CellTypes_RNA", embedding = 'UMAPHarmony')
plotPDF(p1, 
        name = "Plot-UMAP-Celltypes-scRNA-Integrated.pdf", 
        ArchRProj = org, 
        addDOC = FALSE, width = 5, height = 5)

#Based on expression patterns of key retinal genes (and info from snRNA-seq integration), add cell class id's to ArchR object
cM <- confusionMatrix(org$ClustersHarmony, org$cellNames)
labelOld <- rownames(cM)
labelOld
pal_org = c('#154360', '#145A32', '#AD1457', '#B71C1C', '#F4D03F', '#5DADE2',
            '#4CAF50', '#EC7063', '#E67E22', '#B7950B', '#7D3C98', '#BDC3C7')
labelNew <- c('Cones', 'Muller Glia', 'Rods', 'Bipolar Cells', 'Amacrine/Horizontal Cells', 'Late RPCs', 'Early RPCs', 'Late RPCs', 'AC/HC Precursors', 'PR/BC Precursors', 'Early RPCs', 'Developing RGCs', 'RGCs', 'RPE', 'PR/BC Precursors')
org$CellTypes <- mapLabels(org$ClustersHarmony, newLabels = labelNew, oldLabels = labelOld)
p1 = plotEmbedding(ArchRProj = org, colorBy = "cellColData", name = "CellTypes", embedding = "UMAPHarmony", pal = pal_org, rastr = F)
plotPDF(p1, name = "Plot-UMAP-Celltypes.pdf", ArchRProj = org, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = org, outputDirectory = "/outputdirectory/", load = TRUE)

#get cell ids in order to split bam files for correlations
ids = rownames(getCellColData(org))
cell_class = getCellColData(org)$CellTypes
cell.info = data.frame(ids,cell_class)
write.csv(cell.info, file= '/outputdirectory/org_all.cell_class_info.101221.csv')

#Make pseudo-bulk replicates
org <- addGroupCoverages(ArchRProj = org, groupBy = "CellTypes")

#Call peaks with macs2
pathToMacs2 <- findMacs2()
org <- addReproduciblePeakSet(
  ArchRProj = org, 
  groupBy = "CellTypes", 
  pathToMacs2 = pathToMacs2,
  cutOff = 0.000001
)
getPeakSet(org)

#Add peak matrix
org = addPeakMatrix(org)

##Identifying Marker peaks between cell classes
#get peaks specific to cell class
markersPeaks <- getMarkerFeatures(
  ArchRProj = org, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellTypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

##graph marker peaks
heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks, 
                                  cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "heatmap_cell_class_marker_peaks.macs2_q0.000001.pdf", width = 8, height = 6, ArchRProj = org, addDOC = FALSE)
##instead of drawing plot make a matrix
heatmap.matrix <- plotMarkerHeatmap(seMarker = markersPeaks, 
                                    cutOff = "FDR <= 0.01 & Log2FC >= 1", returnMat = TRUE)
##get csv file from matrix
write.csv(heatmap.matrix, file='/outputdirectory/heatmap_cell_class_marker_peaks.macs2_q0.000001.matrix.csv')
##get all peaks 
heatmap.matrix <- plotMarkerHeatmap(seMarker = markersPeaks, 
                                    cutOff = "FDR <= 1000 & Log2FC >= -1000", returnMat = TRUE)
##get csv file from matrix
write.csv(heatmap.matrix, file='/outputdirectory/heatmap_cell_class_all_peaks.macs2_q0.000001.matrix.csv')

saveArchRProject(ArchRProj = org, outputDirectory = "/outputdirectory/", load = TRUE)

#Motif enrichment
org <- addMotifAnnotations(ArchRProj = org, motifSet = "cisbp", name = "Motif")

#Motif enrichment in marker peaks
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = org,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = org, addDOC = FALSE)

#ChromVAR deviations enrichment
org <- addBgdPeaks(org)
org <- addDeviationsMatrix(
  ArchRProj = org, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(org, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = org, addDOC = FALSE)

#Get peak to gene data and atac/rna index info - repeat and save in current directory
org <- addPeak2GeneLinks(
  ArchRProj = org,
  reducedDims = "Harmony",
  maxDist = 300000,
)

p2g <- getPeak2GeneLinks(
  ArchRProj = org,
  corCutOff = 0.4,
  resolution = 10000,
  returnLoops = F
)

write.table(p2g, '/outputdirectory/organoid_all.min_cell20.p2g.m2_q0.000001.cor4.txt', row.names = F, sep="\t", quote = FALSE)
rnaidx.info = metadata(p2g)[[2]]
rnaidx.df = data.frame(rnaidx.info)
write.table(rnaidx.df, '/outputdirectory/organoid_all.min_cell20.m2_q0.000001.rnaseq_info.txt', row.names = T, sep="\t", quote = FALSE)
peaks.gr <- getPeakSet(org)
df.peaks.gr = data.frame(peaks.gr)
write.table(df.peaks.gr, '/outputdirectory/organoid_all.min_cell20.peak_info.m2_q0.000001.txt', row.names = T, sep="\t", quote = FALSE)

saveArchRProject(ArchRProj = org, outputDirectory = "/outputdirectory/", load = TRUE)
