#Examine cell class-specific differential gene expression (from enhancer wt vs ko organoid snRNA-seq dataset) as volcano plots

#R version 3.6.3

# Load necessary libraries
library(dplyr) #version 1.0.4
library(uwot) #version 0.1.9
library(patchwork) #version 1.1.0
library(ggplot2) #version 3.3.2
library(ggrepel) #version 0.8.2
library(EnhancedVolcano) #version 1.4.0

#Adjust maximum memory usage to fit memory allotted to session (200 GB) 
options(future.globals.maxSize = 200 * 1024 ^ 3)

# Generate plot for differentially expressed genes in knockout Muller glia:

#Load input files
input = read.csv('/filepath/mullers_ko.csv') #output of Seurat's differential expression, subset by KO
#NOTE also that all genes in the differential expression list with an adjusted p-value of 0 were set to 2.23E-308, which is the smallest value R can calculate, in order to avoid undefined numbers due to log transformation of p-values for volcano plot
targets = read.csv(file = '/filepath/mactel_nearby_genes.csv') #table of genes located within presumptive TAD of enhancer, generated via UCSC genome browser
targets2 = read.csv(file = '/filepath/mir9_predicted.csv') #table of predicted miR-9-2 target genes, generated via TargetScan

#Manipulate dataframes for downstream analysis
input$mir9_target = ifelse(input$gene %in% targets2$predicted_target, 'yes', 'no')
input$enh_target = ifelse(input$gene %in% targets$gene, 'yes', 'no')
input$mir9_up = ifelse(input$mir9_target == 'yes' & input$log2FC >= 0.5, 'yes', 'no')
targs3 = c('FLT1', 'SEMA3A', 'COL18A1') #in order to label interesting vascular genes
input$endo = ifelse(input$gene %in% targs3, 'yes', 'no')

#Set color and label parameters for volcano plot
keyvals = ifelse(input$enh_target == 'yes', 'darkgreen', ifelse(input$mir9_up == 'yes', '#7D3C98',  ifelse(input$endo == 'yes', '#F4D03F', ifelse(input$p_val_adj > 0.00001, 'gray29', ifelse(input$log2FC <= 0.5 & input$log2FC >= -0.5, '#21618C', '#B71C1C'))))) 
names(keyvals)[keyvals == 'darkgreen'] <- 'Potential e5q14.3 Target'
names(keyvals)[keyvals == '#7D3C98'] <- 'Predicted miR-9-2 Target'
names(keyvals)[keyvals == '#B71C1C'] <- 'Sig. p-value & log2FC'
names(keyvals)[keyvals == '#21618C'] <- 'Sig. p-value only'
names(keyvals)[keyvals == '#F4D03F'] <- 'Gliovascular Signaling Genes '
names(keyvals)[keyvals == 'gray29'] <- 'NS'

#Define lists of genes to be labeled in plot
enh = input[input$enh_target == 'yes',]
mir9 = input[input$mir9_up == 'yes',]
targs = as.character(enh$gene)
targs2 = as.character(mir9$gene)
targs3 = c('FLT1', 'SEMA3A', 'COL18A1')
targs = c(targs, targs2, targs3)

#Generate labeled volcano plot, with genes of interest labeled, colored, and increased in size
EnhancedVolcano(input,
                lab = as.character(input$gene),
                x = 'log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = targs,
                FCcutoff = 0.5,
                colAlpha = 3/4,
                colCustom = keyvals,
                labSize = 2,
                pointSize = c(ifelse(input$enh_target == 'yes', 6, ifelse(input$endo == 'yes', 6, ifelse(input$mir9_up == 'yes', 6, 1)))),
                drawConnectors = T,
                legendPosition = 'none',
                legendLabSize = 12,
                legendIconSize = 4.0,
)

# Genrate volcano plot with no labels (so we can label ourselves in Illustrator/Inkscape)
EnhancedVolcano(input,
                lab = as.character(input$gene),
                x = 'log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = targs,
                FCcutoff = 0.5,
                colAlpha = 3/4,
                colCustom = keyvals,
                labSize = 0,
                pointSize = c(ifelse(input$enh_target == 'yes', 6, ifelse(input$endo == 'yes', 6, ifelse(input$mir9_up == 'yes', 6, 1)))),
                legendPosition = 'none',
                legendLabSize = 12,
                legendIconSize = 4.0)

#Make plots for primary and neurogenic progenitors:
#Primary progenitors:
input = read.csv('/filepath/primary_prog_ko.csv')
input$mir9_target = ifelse(input$gene %in% targets2$predicted_target, 'yes', 'no')
input$enh_target = ifelse(input$gene %in% targets$gene, 'yes', 'no')
input$mir9_up = ifelse(input$mir9_target == 'yes' & input$log2FC >= 0.5, 'yes', 'no')

keyvals = ifelse(input$enh_target == 'yes', 'darkgreen', ifelse(input$mir9_up == 'yes', '#7D3C98', ifelse(input$p_val_adj > 0.00001, 'gray29', ifelse(input$log2FC <= 0.5 & input$log2FC >= -0.5, '#21618C', '#B71C1C')))) 
names(keyvals)[keyvals == 'darkgreen'] <- 'Potential e5q14.3 Target'
names(keyvals)[keyvals == '#7D3C98'] <- 'Predicted miR-9-2 Target'
names(keyvals)[keyvals == '#B71C1C'] <- 'Sig. p-value & log2FC'
names(keyvals)[keyvals == '#21618C'] <- 'Sig. p-value only'
names(keyvals)[keyvals == 'gray29'] <- 'NS'

enh = input[input$enh_target == 'yes',]
mir9 = input[input$mir9_up == 'yes',]
targs = as.character(enh$gene)
targs2 = as.character(mir9$gene)
targs = c(targs, targs2)

#Labeled plot
EnhancedVolcano(input,
                lab = as.character(input$gene),
                x = 'log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = targs,
                FCcutoff = 0.5,
                colAlpha = 1/2,
                colCustom = keyvals,
                labSize = 2.0,
                pointSize = c(ifelse(input$enh_target == 'yes', 6, ifelse(input$mir9_up == 'yes', 6, 1))),
                drawConnectors = F,
                legendPosition = 'none',
                legendLabSize = 12,
                legendIconSize = 4.0)
#Unlabeled plot
EnhancedVolcano(input,
                lab = as.character(input$gene),
                x = 'log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = targs,
                FCcutoff = 0.5,
                colAlpha = 3/4,
                colCustom = keyvals,
                labSize = 0,
                pointSize = c(ifelse(input$enh_target == 'yes', 6, ifelse(input$mir9_up == 'yes', 6, 1))),
                legendPosition = 'none',
                legendLabSize = 12,
                legendIconSize = 4.0)

#Neurogenic progenitors:
input = read.csv('/filepath/neuro_prog_ko.csv')
input$mir9_target = ifelse(input$gene %in% targets2$predicted_target, 'yes', 'no')
input$enh_target = ifelse(input$gene %in% targets$gene, 'yes', 'no')
input$mir9_up = ifelse(input$mir9_target == 'yes' & input$log2FC >= 0.5, 'yes', 'no')

keyvals = ifelse(input$enh_target == 'yes', 'darkgreen', ifelse(input$mir9_up == 'yes', '#7D3C98', ifelse(input$p_val_adj > 0.00001, 'gray29', ifelse(input$log2FC <= 0.5 & input$log2FC >= -0.5, '#21618C', '#B71C1C')))) 
names(keyvals)[keyvals == 'darkgreen'] <- 'Potential e5q14.3 Target'
names(keyvals)[keyvals == '#7D3C98'] <- 'Predicted miR-9-2 Target'
names(keyvals)[keyvals == '#B71C1C'] <- 'Sig. p-value & log2FC'
names(keyvals)[keyvals == '#21618C'] <- 'Sig. p-value only'
names(keyvals)[keyvals == 'gray29'] <- 'NS'

enh = input[input$enh_target == 'yes',]
mir9 = input[input$mir9_up == 'yes',]
targs = as.character(enh$gene)
targs2 = as.character(mir9$gene)
targs = c(targs, targs2)

#Labeled plot
EnhancedVolcano(input,
                lab = as.character(input$gene),
                x = 'log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = targs,
                FCcutoff = 0.5,
                colAlpha = 1/2,
                colCustom = keyvals,
                labSize = 2.0,
                pointSize = c(ifelse(input$enh_target == 'yes', 6, ifelse(input$mir9_up == 'yes', 6, 1))),
                drawConnectors = F,
                legendPosition = 'none',
                legendLabSize = 12,
                legendIconSize = 4.0)
#Unlabeled plot
EnhancedVolcano(input,
                lab = as.character(input$gene),
                x = 'log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                selectLab = targs,
                FCcutoff = 0.5,
                colAlpha = 3/4,
                colCustom = keyvals,
                labSize = 0,
                pointSize = c(ifelse(input$enh_target == 'yes', 6, ifelse(input$mir9_up == 'yes', 6, 1))),
                legendPosition = 'none',
                legendLabSize = 12,
                legendIconSize = 4.0)
