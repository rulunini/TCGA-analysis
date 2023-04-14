
# 1_correlation.R

.libPaths('/home/sohee/myProject/Library')
# devtools::install_github('dviraran/xCell')
library(xCell)
library(dplyr)
library(reshape2)
library(ggplot2)
# xcellnalysis ( 한 환자별, cellytpe알 수 잇음)
# 
# HNSC: TCGA data
# clinical: 환자정보
# 
# neg, pos만 
# interminate (빼고 사용)

setwd('/home/data/HNSCC/mydata/tcga/')
load('TCGAdata.RData')
save_path = '/home/data/HNSCC/results/tcga/'

# (1) TCGA data -> cell type annotation
xCellres <- xCellAnalysis(HNSC)
xCellres[rownames(xCellres) %in% 'Adipocytes',]

xCell.neg<- xCellres[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']
xCell.pos<- xCellres[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='positive']

# (2) HPV- 환자만 추출
HNSC.neg <- HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']
HNSC.pos <- HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='positive']

# (3) gene ratio in TCGA ----------------------------------------------------------
gene <- c('THBS3','AGRN','FN1','SELPLG')
df.all = data.frame()
for (i in 1:length(gene)){
  neg.df <- data.frame(melt(HNSC[rownames(HNSC) == gene[i], clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']), id = 'HPV-', ligand = gene[i])
  pos.df <- data.frame(melt(HNSC[rownames(HNSC) == gene[i], clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='positive']), id = 'HPV+', ligand = gene[i])
  
  df <- cbind(mean(neg.df$value), mean(pos.df$value))
  colnames(df)<- c('HPV-', 'HPV+')
  rownames(df) <- gene[i]
  
  df <- prop.table(df)*100
  df.all <- rbind(df.all, df)
}
df.all

df.all$ligand <- gene
df.all$ligand <- factor(df.all$ligand)

# levels(df.all$ligand)[levels(df.all$ligand) == 'THBS3'] <- 'TCGA-THBS3'
# levels(df.all$ligand)[levels(df.all$ligand) == 'AGRN'] <- 'TCGA-AGRN'
# levels(df.all$ligand)[levels(df.all$ligand) == 'FN1'] <- 'TCGA-FN1'
# levels(df.all$ligand)[levels(df.all$ligand) == 'SELPLG'] <- 'TCGA-SELPLG'
df.all$ligand<- factor(df.all$ligand, levels = c('THBS3','AGRN','FN1','SELPLG'))

df.neg <- df.all[1:3,]
df.pos <- df.all[4,]
df.neg <- melt(df.neg)
df.pos <- melt(df.pos)


b <- ggplot(df.pos, aes(variable, value ,fill = variable))+
  geom_bar(stat = 'identity')+
  facet_wrap(.~ligand)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = '',
        strip.placement = '',
        # panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c('#008000','#A9D18E'))+
  xlab('')+
  ylab('Proportion (%)')+
  ggtitle('TCGA HNSCC')
b
tiff(filename=paste0(save_path,"barplot_neg.tiff"), width = 80, height = 50, unit = "mm", bg = "transparent", res = 300)
a
dev.off()

tiff(filename=paste0(save_path,"barplot_pos.tiff"), width = 35, height = 50, unit = "mm", bg = "transparent", res = 300)
b
dev.off()

# (4) correlation analysis: between 'CAF' and 'ligand expression'  -------------
rownames(xCellres)

gene <- c('THBS3','AGRN','FN1','SELPLG')
df.all = data.frame()
for (i in 1:length(gene)){
  neg.df <- data.frame(melt(HNSC[rownames(HNSC) == gene[i], clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']), id = 'HPV-', ligand = gene[i])
  pos.df <- data.frame(melt(HNSC[rownames(HNSC) == gene[i], clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='positive']), id = 'HPV+', ligand = gene[i])
  df <- rbind(neg.df, pos.df)
  
  df.all <- rbind(df.all, df)
}

colnames(df.all) <- c('barcode','value','HPV','ligand')
# y축: ligand 발현량
# x축: CAF비율

head(df.all)
lig.df <- data.frame(t(acast(df.all, ligand~barcode)))
lig.df$barcode <- rownames(lig.df)
lig.df <- lig.df[sort(lig.df$barcode),]

cell.df <- data.frame(xCellres[rownames(xCellres)=='DC',colnames(xCellres) %in% df.all$barcode])
cell.df$barcode <- rownames(cell.df)
cell.df <- cell.df[sort(cell.df$barcode),]

lig.df$cell <- cell.df$xCellres.rownames.xCellres......DC...colnames.xCellres...in..
rownames(lig.df) <- NULL

lig.neg <- lig.df[lig.df$barcode %in% neg.df$variable,]
lig.neg$HPV <- 'HPV-'
lig.pos <- lig.df[lig.df$barcode %in% pos.df$variable,]
lig.pos$HPV <- 'HPV+'

lig.df <- rbind(lig.neg, lig.pos)
rownames(lig.df) <- lig.df$barcode
lig.df$barcode <- NULL

head(lig.df)

library(PerformanceAnalytics)

test <- lig.df[,c('AGRN','cell','HPV')]
plot(test)
chart.Correlation(test, histogram=T, pch = 19)


tiff(filename=paste0(save_path,"scatter_agrn.tiff"), width = 100, height = 60, unit = "mm", bg = "transparent", res = 300)

ggplot(lig.df, aes(cell, AGRN)) +
  geom_point() + geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size= 10),
        strip.background = element_blank())+
  ylab('Expression level')+
  xlab('Fibroblast expression')+
  ggtitle('AGRN')+
  facet_wrap(.~HPV, ncol = 1)

dev.off()

tiff(filename=paste0(save_path,"scatter_selplg.tiff"), width = 100, height = 60, unit = "mm", bg = "transparent", res = 300)

ggplot(lig.df, aes(cell, SELPLG)) +
  geom_point() + geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())+
  ylab('Expression level')+
  xlab('Dendritic expression')+
  ggtitle('SELPLG')+
  facet_wrap(.~HPV)
dev.off()
