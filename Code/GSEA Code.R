#Gene set enrichment analysis of differentially expressed genes in enhancer wt vs ko organoid snRNA-seq dataset

#R version 4.0.4

#Load necessary libraries
library(Seurat) #version 4.0.1
library(presto) #version 1.0.0
library(dplyr) #version 1.0.6
library(tidyverse) #version 1.3.1
library(ggplot2) #version 3.3.5
library(fgsea) #version 1.16.0

#Adjust maximum memory usage to fit memory allotted to session (200 GB) 
options(future.globals.maxSize = 200 * 1024 ^ 3)

#Load organoid Seurat object (saved from Seurat analysis pipeline)
org <- readRDS("/filepath/org.rds")

#Subset object by celltype, to pull out Muller glia and week 12 (primary, neurogenic) progenitors
Idents(org) = 'celltype'
mg = subset(org, idents = 'Muller Glia')
pp = subset(org, idents = 'Primary Progenitors')
np = subset(org, idents = 'Neurogenic Progenitors')

#Use presto to determine differentially expressed genes (faster computation, and returns an AUC statistic)
de.mg = wilcoxauc(mg, group_by = 'genotype', assay = 'data', seurat_assay = 'RNA')
de.pp = wilcoxauc(pp, group_by = 'genotype', assay = 'data', seurat_assay = 'RNA')
de.np = wilcoxauc(np, group_by = 'genotype', assay = 'data', seurat_assay = 'RNA')

#Subset differential expression output to limit to genes in the knockout
mg.ko <- de.mg %>%
  dplyr::filter(group == "KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pp.ko <- de.pp %>%
  dplyr::filter(group == "KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
np.ko <- de.np %>%
  dplyr::filter(group == "KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

#Generate ranked list of genes by AUC
ranks.mg = deframe(mg.ko)
ranks.pp = deframe(pp.ko)
ranks.np = deframe(np.ko)

#Load Hallmark annotated pathways (downloaded from MSigDB)
hallmark <- gmtPathways("/filepath/h.all.v7.4.symbols.gmt")

#Run GSEA
gsea.mg <- fgseaMultilevel(hallmark, stats = ranks.mg, eps = 0)
gsea.pp <- fgseaMultilevel(hallmark, stats = ranks.pp, eps = 0)
gsea.np <- fgseaMultilevel(hallmark, stats = ranks.np, eps = 0)

#Tidy the data and plot
mg.tidy <- gsea.mg %>%
  as_tibble() %>%
  arrange(desc(NES))
np.tidy <- gsea.np %>%
  as_tibble() %>%
  arrange(desc(NES))
pp.tidy <- gsea.pp %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(mg.tidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways - Muller Glia") + 
  theme_minimal()

ggplot(np.tidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways - Neurogenic RPCs") + 
  theme_minimal()

ggplot(pp.tidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways - Primary RPCs") + 
  theme_minimal()

