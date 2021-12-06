##load libraries
library("pheatmap")
library("ggplot2")
library('dendsort')

##function for sort clusters
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

##graph human/mouse data -- removing hu_Mature Ganglions
hm_hm <- read.table('human_mouse.ret_genes.heatmap_all.int.txt', row.names=1, header=T, sep='\t')
hm_hm = subset(hm_hm, select=-c(hu_Mature.Ganglions))
pheatmap(hm_hm, show_rownames = F)
dev.copy2pdf(file='human_mouse.ret_genes.heatmap_all.no_RGC.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(hm_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(hm_hm)))
pheatmap(hm_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_mouse.ret_genes.heatmap_all.no_RGC.cluster2.pdf', width = 9, height = 6)

##graph human/organoid data -- removing hu_Mature Ganglions
ho_hm <- read.table('human_organoid.ret_genes.heatmap_all.int.txt', row.names=1, header=T, sep='\t')
ho_hm = subset(ho_hm, select=-c(hu_Mature.Ganglions))
pheatmap(ho_hm, show_rownames = F)
dev.copy2pdf(file='human_org.ret_genes.heatmap_all.no_RGC.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(ho_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(ho_hm)))
pheatmap(ho_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_org.ret_genes.heatmap_all.no_RGC.cluster2.pdf', width = 9, height = 6)


