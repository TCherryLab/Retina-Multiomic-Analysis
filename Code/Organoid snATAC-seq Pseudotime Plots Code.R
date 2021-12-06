#Generate pseudotime plots for accessibility/expression of key retinal cell class genes in organoid snATAC-seq dataset

#R version 3.6.3

#Load necessary libraries
library(ArchR) #version 0.9.5

#Load organoid snATAC-seq ArchR project
proj1 = loadArchRProject(path = '/filepath/')

## Rod genes:

#Establish order of trajectory
rodtrajectory = c('Early RPCs', 'Late RPCs', 'PR/BC Precursors', 'Rods')
#Add trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "rodU", 
                       groupBy = "CellTypes",trajectory = rodtrajectory, 
                       reducedDims = "Harmony", embedding = "UMAPHarmony", force = TRUE)
#Plot pseudotime trajectory for relevant genes for both GeneScore and GeneIntegration
p1 <- plotTrajectory(proj1, trajectory = "rodU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "NRL", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "NRL", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "rodU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "NOTCH1", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "rodU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "NOTCH1", continuousSet = "blueYellow",  log2norm = TRUE)
#Save all plots in one giant pdf
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]],
        name = "trajectory_rod_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

## Cone genes
conetrajectory = c('Early RPCs', 'PR/BC Precursors', 'Cones')
proj1 <- addTrajectory(ArchRProj = proj1, name = "coneU", 
                       groupBy = "CellTypes",trajectory = conetrajectory, 
                       reducedDims = "Harmony", embedding = "UMAPHarmony", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "coneU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "ARR3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "ARR3", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "coneU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "NRL", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "coneU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "NRL", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]],
        name = "trajectory_cone_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

## Horizontal/amacrine genes
achctrajectory = c('Early RPCs', 'AC/HC Precursors', 'Amacrine/Horizontal Cells')
proj1 <- addTrajectory(ArchRProj = proj1, name = "achcU", 
                       groupBy = "CellTypes",trajectory = achctrajectory, 
                       reducedDims = "Harmony", embedding = "UMAPHarmony", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "achcU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "TFAP2B", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "achcU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "TFAP2B", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "achcU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "ONECUT2", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "achcU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "ONECUT2", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "achcU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "FOXN4", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "achcU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "FOXN4", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]],
        name = "trajectory_amacrine_horizontal_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


## Bipolar genes
bptrajectory = c('Early RPCs', 'Late RPCs', 'PR/BC Precursors', 'Bipolar Cells')
proj1 <- addTrajectory(ArchRProj = proj1, name = "cone_bipolarU", 
                       groupBy = "CellTypes",trajectory = bptrajectory, 
                       reducedDims = "Harmony", embedding = "UMAPHarmony", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "VSX1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "VSX1", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "PRDM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "PRDM1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]],
        name = "trajectory_bipolar_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)



## Ganglion genes
rgctrajectory = c('Early RPCs', 'Developing RGCs', 'RGCs')
proj1 <- addTrajectory(ArchRProj = proj1, name = "ganglion2U", 
                       groupBy = "CellTypes",trajectory = rgctrajectory, 
                       reducedDims = "Harmony", embedding = "UMAPHarmony", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "POU4F1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "POU4F1", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "ganglion2U", embedding = 'UMAPHarmony', colorBy = "GeneScoreMatrix", name = "ATOH7", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "ganglion2U", embedding = 'UMAPHarmony', colorBy = "GeneIntegrationMatrix", name = "ATOH7", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]],
        name = "trajectory_ganglion_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


## Muller glia genes
mgtrajectory = c('Early RPCs', 'Late RPCs', 'Muller Glia')
proj1 <- addTrajectory(ArchRProj = proj1, name = "mullerU", 
                       groupBy = "CellTypes",trajectory = mgtrajectory, 
                       reducedDims = "Harmony", embedding = "UMAPHarmony", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", embedding = "UMAPHarmony", colorBy = "GeneScoreMatrix", name = "HES1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", embedding = "UMAPHarmony", colorBy = "GeneIntegrationMatrix", name = "HES1", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "mullerU", embedding = "UMAPHarmony", colorBy = "GeneScoreMatrix", name = "SOX9", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "mullerU", embedding = "UMAPHarmony", colorBy = "GeneIntegrationMatrix", name = "SOX9", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]],
        name = "trajectory_muller_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

