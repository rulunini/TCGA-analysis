# 6_chisq.R

.libPaths('/home/sohee/myProject/Library')
# devtools::install_github('dviraran/xCell')
library(xCell)
library(dplyr)
library(reshape)
library(ggplot2)
library(survival)
library(gmodels)
library(reshape2)
# xcellnalysis ( 한 환자별, cellytpe알 수 잇음)
# 
# HNSC: TCGA data
# clinical: 환자정보
# 
# neg, pos만 
# interminate (빼고 사용)

setwd('/home/sohee/analysis/data/HNSCC/mydata/tcga/')
load('TCGAdata.RData') 

# (1) TCGA data -> cell type annotation ----------
xCellres <- xCellAnalysis(HNSC)
xCellres[rownames(xCellres) %in% 'Adipocytes',]

xCell.neg<- xCellres[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']
xCell.pos<- xCellres[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='positive']


# (2) HPV- 환자만 추출 ----------
HNSC.neg <- HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']


# (3) my gene, stage 정보 포함한 dataframe 만들기 ----------
my.gene <- c('EFNA1','EFNA3','EPHA2','EPHA3')
my.gene <- c('TGFB1','ACVR1','TGFBR1','TGFBR2','ACVR1B')

i = 1
my.df <- cbind(t(HNSC.neg[rownames(HNSC.neg) %in% my.gene[i],]),
               clinical[clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative',][9])
my.df$barcode <- rownames(my.df)
rownames(my.df)<- NULL
colnames(my.df)[1] <- 'gene'
colnames(my.df)[2] <- 'stage'
my.df <- na.omit(my.df)
head(my.df)

my.df$mean <- ifelse(my.df$gene > mean(my.df$gene), 'high', 'low')
my.df$median <- ifelse(my.df$gene > median(my.df$gene), 'high', 'low')


# count table -------------------------

df.mean <- acast(my.df, stage~mean)
chisq.test(df.mean)
chisq.test(rbind(colSums(df.mean[1:3,]),df.mean[4:6,]))
chisq.test(rbind(df.mean[1:3,],colSums(df.mean[4:6,])))
chisq.test(rbind(colSums(df.mean[1:3,]),colSums(df.mean[4:6,])))


df.median <- acast(my.df, stage~median)

chisq.test(df.median)
chisq.test(rbind(colSums(df.median[1:3,]),df.median[4:6,]))
chisq.test(rbind(df.median[1:3,],colSums(df.median[4:6,])))
chisq.test(rbind(colSums(df.median[1:3,]),colSums(df.median[4:6,])))


# proportion table -----------------

d1 <- melt(prop.table(df.mean, margin = 2)*100)

ggplot(d1, aes(Var1, value, group = Var2, color = Var2))+
  geom_line()+
  geom_point()+
  theme_bw()+
  ylab('count')+
  xlab('')+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = 'bottom')+
  ggtitle(my.gene[i])
# facet_wrap(.~gene)

ggplot(d1, aes(Var1, value, fill = Var2))+
  geom_bar(stat = 'identity', position = 'dodge')+
  # facet_wrap(.~gene)
  theme_bw()+
  ylab('proportion')+
  xlab('')+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = 'bottom')+
  ggtitle(my.gene[i])















# (4) Visualization ----------
library(PerformanceAnalytics)

ggplot(my.df2, aes(stage, value))+
  geom_point()#+geom_smooth(method = "lm")+
facet_wrap(.~variable) # , ncol = 1

# stat_cor(method = "pearson")+
theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size= 10),
        strip.background = element_blank())
ylab('Expression level')+
  xlab('Fibroblast expression')+
  ggtitle('AGRN')+
  