# Goal: Stage와 특정 gene의 연관성.

# We need the two dataset from TCGA.
# First, clinical tables that include the patients information.
# Second, the tables gene expression values (rows) per patients (columns).

# From this two dataset, we will use both cancer stages information and gene expression.

# However, we will not use just gene expression value, but frequencies of cancer stages by high or low gene expression.


# < Load Library > ----------

library(reshape2)
library(dplyr)


# < Load your TCGA data > ----------

my.path <- '/home/sohee/analysis/data/HNSCC/mydata/tcga/'
load(paste0(my.path, 'TCGAdata.RData'))


# < Create my dataframe > ---------
my.df <- HNSC
HNSC.neg <- HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']

my.gene <- c('write_gene_symbols_that_you_are_interesting')
# for example 
# my.gene <- c('EFNA1','EFNA3','EPHA2','EPHA3')
# my.gene <- c('TGFB1','ACVR1','TGFBR1','TGFBR2','ACVR1B')

i = 1

a <- t(my.df[rownames(my.df) %in% my.gene[i],clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative'])
b <- clinical[clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative',][9]
cbind(a, b)
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


# < Convert to count table > ----------
# We will use count data by decided groups in advance.

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


# < Proportion visualization > ----------

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