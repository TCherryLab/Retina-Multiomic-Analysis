##using r4.0.3
library(UpSetR)

##files to graph
file_to_graph = 'diff_peaks_rep.bg.28wk_organoid_adult_human.all_peaks.balanced.bt_int.mature_cc.upset.txt'
file_to_graph = 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc1_p0.05.bt_int.mature_cc.upset.txt'
file_to_graph = 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc2_p0.001.bt_int.mature_cc.upset.txt'
file_to_graph = 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc3_p0.001.bt_int.mature_cc.upset.txt'
file_to_graph = 'diff_peaks_rep.bg.adult_human_28wk_organoid.all_peaks.balanced.bt_int.mature_cc.upset.txt'
file_to_graph = 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.bt_int.mature_cc.upset.txt'
file_to_graph = 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.bt_int.mature_cc.upset.txt'
file_to_graph = 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.bt_int.mature_cc.upset.txt'

##cmds
upset_pdf <- gsub(".txt",".pdf",file_to_graph)
wide_csv = gsub(".txt",".wide.txt",file_to_graph)
##read in data -- adding summary meg/mic to phenotypes
links <- read.table(file_to_graph, header=T, sep='\t')
links_wide = reshape(links, idvar = "peak", timevar = "cell_class", direction = "wide")
write.csv(links_wide, file = wide_csv)
##graph interactions
upset(links_wide, sets = c('count.human.Mature_Amacrines', 'count.human.Mature_Bipolars', 'count.human.Mature_Cones', 'count.human.Mature_Horizontals', 'count.human.Mature_Mullers', 'count.human.Mature_Rods', 'count.organoid.Cones', 'count.organoid.RGCs', 'count.organoid.Rods', 'count.organoid.Amacrine_Horizontal_Cells', 'count.organoid.Bipolar_Cells', 'count.organoid.Muller_Glia'), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on", keep.order = TRUE)
dev.copy2pdf(file=upset_pdf, width = 12, height = 6)



