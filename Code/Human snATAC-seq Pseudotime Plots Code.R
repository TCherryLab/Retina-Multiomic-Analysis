#Generate pseudotime plots for accessibility/expression of key retinal cell class genes in human snATAC-seq dataset

#R version 3.6.3

#Load necessary libraries
library(ArchR) #version 0.9.5

#Load human snATAC-seq ArchR project
proj1 = loadArchRProject(path = '/filepath/')

## Rod genes:

#Establish order of trajectory
rod_trajectory <- c("C20", "C6", "C2", "C1")
#Add trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "rodU", 
                       groupBy = "Clusters",trajectory = rod_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
#Plot pseudotime trajectory for relevant genes for both GeneScore and GeneIntegration
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "RHO", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "RHO", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "MEF2D", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "MEF2D", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "MEF2C", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "MEF2C", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "NRL", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "NRL", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "NR2E3", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "NR2E3", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "SAG", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "SAG", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "HES1", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "HES1", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "NOTCH1", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "NOTCH1", continuousSet = "blueYellow",  log2norm = TRUE)
p21 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p22 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p23 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p24 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p25 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p26 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p27 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "PRDM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p28 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "PRDM1", continuousSet = "blueYellow",  log2norm = TRUE)
#Save all plots in one giant pdf
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], p10[[2]], p11[[2]], p12[[2]], p13[[2]], p14[[2]],
        p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p20[[2]], p21[[2]], p22[[2]], p23[[2]], p25[[2]], p26[[2]], p27[[2]], p28[[2]],
        name = "trajectory_rod_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

## Cone genes
cone_trajectory <- c("C20", "C6", "C7", "C8")
proj1 <- addTrajectory(ArchRProj = proj1, name = "coneU", 
                       groupBy = "Clusters",trajectory = cone_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ARR3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ARR3", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "NRL", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "NRL", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "OPN1SW", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "OPN1SW", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "OPN1LW", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "OPN1LW", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ONECUT1", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ONECUT1", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ONECUT2", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ONECUT2", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ONECUT3", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ONECUT3", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p21 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p22 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p23 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p24 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p25 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "PRDM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p26 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "PRDM1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], p10[[2]], p11[[2]], p12[[2]], p13[[2]], p14[[2]],
        p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p20[[2]], p21[[2]], p23[[2]], p24[[2]], p25[[2]], p26[[2]],
        name = "trajectory_cone_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

## Horizontal genes
horizontal_trajectory <- c("C20", "C5", "C15", "C13")
proj1 <- addTrajectory(ArchRProj = proj1, name = "horizontalU", 
                       groupBy = "Clusters",trajectory = horizontal_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "LHX1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "LHX1", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "ONECUT1", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "ONECUT1", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "ONECUT2", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "ONECUT2", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "ONECUT3", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "ONECUT3", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "POU4F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "POU4F2", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "CALB1", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "CALB1", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "PTF1A", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "PTF1A", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "TFAP2B", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "TFAP2B", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p21 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p22 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p23 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "PROX1", continuousSet = "horizonExtra",  log2norm = TRUE)
p24 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "PROX1", continuousSet = "blueYellow",  log2norm = TRUE)
p25 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "FOXN4", continuousSet = "horizonExtra",  log2norm = TRUE)
p26 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "FOXN4", continuousSet = "blueYellow",  log2norm = TRUE)

plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], p10[[2]], p11[[2]], p12[[2]], p13[[2]], p14[[2]],
        p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p21[[2]], p22[[2]], p23[[2]], p24[[2]], p25[[2]], p26[[2]],
        name = "trajectory_horizontal_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


## Cone bipolar genes
cone_bipolar_trajectory <- c("C20", "C6", "C11", "C10")
proj1 <- addTrajectory(ArchRProj = proj1, name = "cone_bipolarU", 
                       groupBy = "Clusters",trajectory = cone_bipolar_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "PRDM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "PRDM1", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "NEUROD4", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NEUROD4", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "GRM6", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "GRM6", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "TRPM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "TRPM1", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIA", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIA", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIB", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIB", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIX", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIX", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p21 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p22 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p23 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p24 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p25 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "VSX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p26 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "VSX2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], p10[[2]], p11[[2]], p12[[2]], p13[[2]], p14[[2]],
        p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p20[[2]], p21[[2]], p23[[2]], p24[[2]], p25[[2]], p26[[2]],
        name = "trajectory_cone_bipolar_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

## Rod bipolar genes
rod_bipolar_trajectory <- c("C20", "C6", "C11", "C9")
proj1 <- addTrajectory(ArchRProj = proj1, name = "rod_bipolarU", 
                       groupBy = "Clusters",trajectory = rod_bipolar_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "PRDM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "PRDM1", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "NEUROD4", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NEUROD4", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "GRM6", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "GRM6", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "TRPM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "TRPM1", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIA", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIA", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIB", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIB", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIX", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIX", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p21 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p22 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p23 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p24 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p25 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "VSX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p26 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "VSX2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], p10[[2]], p11[[2]], p12[[2]], p13[[2]], p14[[2]],
        p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p20[[2]], p21[[2]], p23[[2]], p24[[2]], p25[[2]], p26[[2]],
        name = "trajectory_rod_bipolar_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)



##looking at the amacrines during development
amacrine_trajectory <- c("C20", "C5", "C16", "C14")
amacrine_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "amacrineU", 
                       groupBy = "Clusters",trajectory = amacrine_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "ONECUT1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "ONECUT1", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "ONECUT2", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "ONECUT2", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "ONECUT3", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "ONECUT3", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "PTF1A", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "PTF1A", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "FOXN4", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "FOXN4", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "PAX6", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "PAX6", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "NEUROD1", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "NEUROD1", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "MEIS1", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "MEIS1", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "MEIS2", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "MEIS2", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "MEIS3", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "MEIS3", continuousSet = "blueYellow",  log2norm = TRUE)
p21 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "GAD1", continuousSet = "horizonExtra",  log2norm = TRUE)
p22 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "GAD1", continuousSet = "blueYellow",  log2norm = TRUE)
p23 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "GAD2", continuousSet = "horizonExtra",  log2norm = TRUE)
p24 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "GAD2", continuousSet = "blueYellow",  log2norm = TRUE)
p25 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "SLC6A9", continuousSet = "horizonExtra",  log2norm = TRUE)
p26 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "SLC6A9", continuousSet = "blueYellow",  log2norm = TRUE)
p27 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p28 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p29 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p30 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p31 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p32 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p33 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "TFAP2A", continuousSet = "horizonExtra",  log2norm = TRUE)
p34 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "TFAP2A", continuousSet = "blueYellow",  log2norm = TRUE)
p35 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "TFAP2B", continuousSet = "horizonExtra",  log2norm = TRUE)
p36 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "TFAP2B", continuousSet = "blueYellow",  log2norm = TRUE)
p37 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "NFIA", continuousSet = "horizonExtra",  log2norm = TRUE)
p38 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "NFIA", continuousSet = "blueYellow",  log2norm = TRUE)

plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], p10[[2]], p11[[2]], p12[[2]], p13[[2]], p14[[2]],
        p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p20[[2]], p21[[2]], p22[[2]], p23[[2]], p24[[2]], p25[[2]], p26[[2]],
        p27[[2]], p28[[2]], p29[[2]], p31[[2]], p32[[2]], p33[[2]], p34[[2]], p35[[2]], p36[[2]], p37[[2]], p38[[2]],
        name = "trajectory_amacrine_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

## Ganglion genes
ganglion_trajectory2 <- c("C19", "C4", "C3", "C17", "C18")
proj1 <- addTrajectory(ArchRProj = proj1, name = "ganglion2U", 
                       groupBy = "Clusters",trajectory = ganglion_trajectory2, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "POU4F1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "POU4F1", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "POU4F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "POU4F2", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "POU4F3", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "POU4F3", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "EBF1", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "EBF1", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "EBF2", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "EBF2", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "EBF3", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "EBF3", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "RBPMS", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "RBPMS", continuousSet = "blueYellow",  log2norm = TRUE)
p21 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "ATOH7", continuousSet = "horizonExtra",  log2norm = TRUE)
p22 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "ATOH7", continuousSet = "blueYellow",  log2norm = TRUE)

plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], 
        p11[[2]], p12[[2]], p13[[2]], p14[[2]], p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p20[[2]], p21[[2]], p22[[2]],
        name = "trajectory_ganglion_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


## Muller glia genes
muller_trajectory <- c("C19", "C20", "C12")
proj1 <- addTrajectory(ArchRProj = proj1, name = "mullerU", 
                       groupBy = "Clusters",trajectory = muller_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "MEIS1", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "MEIS1", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "MEIS2", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "MEIS2", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "MEIS3", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "MEIS3", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "HES1", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "HES1", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "NOTCH1", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "NOTCH1", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "SOX9", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "SOX9", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "RLBP1", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "RLBP1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], p10[[2]], 
        p11[[2]], p12[[2]], p13[[2]], p14[[2]], p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p20[[2]],
        name = "trajectory_muller_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

## Bipolar genes (not split into rod and cone bps)
bipolar_trajectory <- c('Late Progenitors', 'Photoreceptor/Bipolar Precursors', 'Developing Bipolars', 'Mature Bipolars')
proj1 <- addTrajectory(ArchRProj = proj1, name = "bipolarU", 
                       groupBy = "Clusters2",trajectory = bipolar_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
p3 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p4 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
p5 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "PRDM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p6 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "PRDM1", continuousSet = "blueYellow",  log2norm = TRUE)
p7 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "NEUROD4", continuousSet = "horizonExtra",  log2norm = TRUE)
p8 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "NEUROD4", continuousSet = "blueYellow",  log2norm = TRUE)
p9 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "GRM6", continuousSet = "horizonExtra",  log2norm = TRUE)
p10 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "GRM6", continuousSet = "blueYellow",  log2norm = TRUE)
p11 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "TRPM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p12 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "TRPM1", continuousSet = "blueYellow",  log2norm = TRUE)
p13 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "NFIA", continuousSet = "horizonExtra",  log2norm = TRUE)
p14 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIA", continuousSet = "blueYellow",  log2norm = TRUE)
p15 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "NFIB", continuousSet = "horizonExtra",  log2norm = TRUE)
p16 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIB", continuousSet = "blueYellow",  log2norm = TRUE)
p17 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "NFIX", continuousSet = "horizonExtra",  log2norm = TRUE)
p18 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIX", continuousSet = "blueYellow",  log2norm = TRUE)
p19 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p20 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
p21 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p22 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
p23 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p24 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
p25 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "VSX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p26 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "VSX2", continuousSet = "blueYellow",  log2norm = TRUE)
p27 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "VSX1", continuousSet = "horizonExtra",  log2norm = TRUE)
p28 <- plotTrajectory(proj1, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "VSX1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1[[2]], p2[[2]], p3[[2]], p4[[2]], p5[[2]], p6[[2]], p7[[2]], p8[[2]], p9[[2]], p10[[2]], p11[[2]], p12[[2]], p13[[2]], p14[[2]],
        p15[[2]], p16[[2]], p17[[2]], p18[[2]], p19[[2]], p20[[2]], p21[[2]], p23[[2]], p24[[2]], p25[[2]], p26[[2]], p27[[2]], p28[[2]],
        name = "trajectory_bipolar_genes.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
