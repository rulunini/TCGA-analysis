# 4_corr_SDC1_with_immune.R

.libPaths('/home/sohee/myProject/Library')
# devtools::install_github('dviraran/xCell')
library(xCell)
library(dplyr)
library(reshape2)
library(ggplot2)
# xcellnalysis ( 한 환자별, cellytpe알 수 잇음)

# HNSC: TCGA data
# clinical: 환자정보

# neg, pos만 
# interminate (빼고 사용)

setwd('/home/sohee/analysis/data/HNSCC/mydata/tcga/')
load('TCGAdata.RData')
save_path = '/home/sohee/analysis/data/HNSCC/results/tcga/'

# (1) TCGA data -> cell type annotation ----------
xCellres <- xCellAnalysis(HNSC)


# (2) preparing the dataframe  -------------------

xCellres <- as.matrix(xCellres)

my.df <- xCellres[!grepl('B-cells', rownames(xCellres)),]
my.df <- my.df[!grepl('CD4+', rownames(my.df)),]
my.df <- my.df[!grepl('CD8+', rownames(my.df)),]
my.df <- my.df[!grepl('Endothelial', rownames(my.df)),]
my.df <- my.df[!grepl('ImmuneScore', rownames(my.df)),]
my.df <- my.df[!grepl('StromaScore', rownames(my.df)),]
my.df <- my.df[!grepl('MicroenvironmentScore', rownames(my.df)),]

B <- colSums(xCellres[grepl('B-cells', rownames(xCellres)),])
CD4T <- colSums(xCellres[grepl('CD4+', rownames(xCellres)),])
CD8T <- colSums(xCellres[grepl('CD8+', rownames(xCellres)),])
Endothelial <- colSums(xCellres[grepl('Endothelial', rownames(xCellres)),])
my.df <- rbind(my.df, B, CD4T, CD8T, Endothelial)
my.df <- melt(my.df)
colnames(my.df) <- c('celltype','barcode','value')
my.df$HPV <- NA
head(my.df)

neg.barcode <- colnames(HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative'])
pos.barcode <- colnames(HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='positive'])
my.df[my.df$barcode %in% neg.barcode,'HPV'] = 'Negative'
my.df[my.df$barcode %in% pos.barcode,'HPV'] = 'Positive'


# (3) add SDC1 expression value to my dataframe  -------------------

SDC1.df <- t(HNSC[rownames(HNSC) %in% 'SDC1',])

my.ct <- c('B','Plasma cells','CD4T','CD8T','Tregs')
my.df <- my.df[order(my.df$barcode),]
my.df$SDC1 <- NA

for (i in 1:length(my.ct)){
  my.df[my.df$celltype %in% my.ct[i],]$SDC1 <- SDC1.df
}

cell.df <- my.df[my.df$celltype %in% my.ct[1:5],]
cell.df <- cell.df[is.na(cell.df$HPV) == FALSE,]

# (4) Visualizing correlation plot  -------------------

tiff(filename=paste0(save_path,"scatter_agrn.tiff"), width = 100, height = 60, unit = "mm", bg = "transparent", res = 300)

ggplot(cell.df, aes(value, SDC1)) +
  geom_point() + geom_smooth(mapping = aes(x = value, y = SDC1), method = "lm")+
  stat_cor(method = "pearson")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size= 10),
        strip.background = element_blank())+
  facet_wrap(celltype~HPV, scale = 'free')

dev.off()
