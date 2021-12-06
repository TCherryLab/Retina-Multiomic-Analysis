# Generate scatterplots for pairwise correlations of human and organoid mature cell classes

#R version 3.6.3

#Load necessary libraries
library(ggplot2) #version 3.3.2
library(cowplot) #version 1.1.0

#Amacrine cells:


#Load in table of correlation values (generated in a previous analysis)
ac = read.delim(file = '/active/cherry_t/OrgManuscript_SingleCell_Data/pairwise_correlations_1021/amacrines_correlation.100d.nonesize.macs2_q0.000001_summits.txt')

#Plot as a scatter plot
p1 = ggplot(ac, aes(x=human.Mature_Amacrines, y = organoid.Amacrine_Horizontal_Cells))+geom_point(size=0.5, color = '#145A32')+
  geom_smooth(method=lm, color = 'black', se=F, fullrange =T)+theme_classic()+
  theme(axis.text.x= element_text(size=10))+theme(axis.text.y= element_text(size=10))+
  scale_x_continuous(limits = c(0,600))+scale_y_continuous(limits = c(0,600))+theme(aspect.ratio = 1)+
  xlab('Human mAC')+ylab('Organoid AC/HC')

#Determine r value
cor(ac$human.Mature_Amacrines, ac$organoid.Amacrine_Horizontal_Cells)

#Bipolar cells:
bc = read.delim(file = '/active/cherry_t/OrgManuscript_SingleCell_Data/pairwise_correlations_1021/bipolars_correlation.100d.nonesize.macs2_q0.000001_summits.txt')
p2 = ggplot(bc, aes(x=human.Mature_Bipolars, y = organoid.Bipolar_Cells))+geom_point(size=0.5, color = '#AD1457')+
  geom_smooth(method=lm, color = 'black', se=F, fullrange =T)+theme_classic()+
  theme(axis.text.x= element_text(size=10))+theme(axis.text.y= element_text(size=10))+
  scale_x_continuous(limits = c(0,600))+scale_y_continuous(limits = c(0,600))+theme(aspect.ratio = 1)+
  xlab('Human mBC')+ylab('Organoid BC')
cor(bc$human.Mature_Bipolars, bc$organoid.Bipolar_Cells)

#Cones:
cones = read.delim(file = '/active/cherry_t/OrgManuscript_SingleCell_Data/pairwise_correlations_1021/cones_correlation.100d.nonesize.macs2_q0.000001_summits.txt')
p3 = ggplot(cones, aes(x=human.Mature_Cones, y = organoid.Cones))+geom_point(size=0.5, color = '#B71C1C')+
  geom_smooth(method=lm, color = 'black', se=F, fullrange =T)+theme_classic()+
  theme(axis.text.x= element_text(size=10))+theme(axis.text.y= element_text(size=10))+
  scale_x_continuous(limits = c(0,1000))+scale_y_continuous(limits = c(0,1000))+theme(aspect.ratio = 1)+
  xlab('Human mCones')+ylab('Organoid Cones')
cor(cones$human.Mature_Cones, cones$organoid.Cones)

#Early Progenitors:
ep = read.delim(file = '/active/cherry_t/OrgManuscript_SingleCell_Data/pairwise_correlations_1021/early_progenitor_correlation.100d.nonesize.macs2_q0.000001_summits.txt')
p4 = ggplot(ep, aes(x=human.Early_Progenitors, y = organoid.Early_RPCs))+geom_point(size=0.5, color = '#5DADE2')+
  geom_smooth(method=lm, color = 'black', se=F, fullrange = T)+theme_classic()+
  theme(axis.text.x= element_text(size=10))+theme(axis.text.y= element_text(size=10))+
  scale_x_continuous(limits = c(0,900))+scale_y_continuous(limits = c(0,900))+theme(aspect.ratio = 1)+
  xlab('Human Early RPCs')+ylab('Organoid Early RPCs')
cor(ep$human.Early_Progenitors, ep$organoid.Early_RPCs)

#Horizontal cells:
hc = read.delim(file = '/active/cherry_t/OrgManuscript_SingleCell_Data/pairwise_correlations_1021/horizontals_correlation.100d.nonesize.macs2_q0.000001_summits.txt')
p5 = ggplot(hc, aes(x=human.Mature_Horizontals, y = organoid.Amacrine_Horizontal_Cells))+geom_point(size=0.5, color = '#154360')+
  geom_smooth(method=lm, color = 'black', se=F, fullrange =T)+theme_classic()+
  theme(axis.text.x= element_text(size=10))+theme(axis.text.y= element_text(size=10))+
  scale_x_continuous(limits = c(0,800))+scale_y_continuous(limits = c(0,800))+theme(aspect.ratio = 1)+
  xlab('Human mHC')+ylab('Organoid AC/HC')
cor(hc$human.Mature_Horizontals, hc$organoid.Amacrine_Horizontal_Cells)

#Late progenitors:
lp = read.delim(file = '/active/cherry_t/OrgManuscript_SingleCell_Data/pairwise_correlations_1021/late_progenitor_correlation.100d.nonesize.macs2_q0.000001_summits.txt')
p6 = ggplot(lp, aes(x=human.Late_Progenitors, y = organoid.Late_RPCs))+geom_point(size=0.5, color = '#4CAF50')+
  geom_smooth(method=lm, color = 'black', se=F, fullrange = T)+theme_classic()+
  theme(axis.text.x= element_text(size=10))+theme(axis.text.y= element_text(size=10))+
  scale_x_continuous(limits = c(0,600))+scale_y_continuous(limits = c(0,600))+theme(aspect.ratio = 1)+
  xlab('Human Late RPCs')+ylab('Organoid Late RPCs')
cor(lp$human.Late_Progenitors, lp$organoid.Late_RPCs)

#Muller glia:
mg = read.delim(file = '/active/cherry_t/OrgManuscript_SingleCell_Data/pairwise_correlations_1021/mullers_correlation.100d.nonesize.macs2_q0.000001_summits.txt')
p7 = ggplot(mg, aes(x=human.Mature_Mullers, y = organoid.Muller_Glia))+geom_point(size=0.5, color = '#EC7063')+
  geom_smooth(method=lm, color = 'black', se=F, fullrange =T)+theme_classic()+
  theme(axis.text.x= element_text(size=10))+theme(axis.text.y= element_text(size=10))+
  scale_x_continuous(limits = c(0,1200))+scale_y_continuous(limits = c(0,1200))+theme(aspect.ratio = 1)+
  xlab('Human MG')+ylab('Organoid MG')
cor(mg$human.Mature_Mullers, mg$organoid.Muller_Glia)

#Rods:
rods = read.delim(file = '/active/cherry_t/OrgManuscript_SingleCell_Data/pairwise_correlations_1021/rods_correlation.100d.nonesize.macs2_q0.000001_summits.txt')
p8 = ggplot(rods, aes(x=human.Mature_Rods, y = organoid.Rods))+geom_point(size=0.5, color = '#7D3C98')+
  geom_smooth(method=lm, color = 'black', se=F, fullrange = T)+theme_classic()+
  theme(axis.text.x= element_text(size=10))+theme(axis.text.y= element_text(size=10))+
  scale_x_continuous(limits = c(0,1500))+scale_y_continuous(limits = c(0,1500))+theme(aspect.ratio = 1)+
  xlab('Human mRods')+ylab('Organoid Rods')
cor(rods$human.Mature_Rods, rods$organoid.Rods)

#Plot all plots on a grid
plot_grid(p4, p6, p5, p1, p8, p3, p7, p2, ncol=4)
