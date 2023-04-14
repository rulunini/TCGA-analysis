# 5_survival.R

.libPaths('/home/sohee/myProject/Library')
# devtools::install_github('dviraran/xCell')
library(xCell)
library(dplyr)
library(reshape)
library(ggplot2)
library(survival)
library(survminer)
library(gmodels)
# xcellnalysis ( �� ȯ�ں�, cellytpe�� �� ����)
# 
# HNSC: TCGA data
# clinical: ȯ������
# 
# neg, pos�� 
# interminate (���� ���)

# https://www.cancer.gov/about-cancer/diagnosis-staging/staging -> TNM �����ؼ�
# https://www.ncbi.nlm.nih.gov/books/NBK65719/table/CDR0000062913__892/ -> TNM �����ؼ� (iva, ivb, ivc)
# https://bioinformaticsandme.tistory.com/224 -> ����� Ʃ�丮��
# https://mustlearning.tistory.com/18 -> ����� �׷���

setwd('/home/sohee/analysis/data/HNSCC/mydata/tcga/')
load('TCGAdata.RData') 

save_path = '/home/data/HNSCC/results/tcga/'

# (1) ������ ��ó�� -----------

# df <- clinical[clinical$patient.hpv_test_results.hpv_test_result.hpv_status == 'negative',]
# HNSC.neg <- HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']

df <- clinical[clinical$patient.hpv_test_results.hpv_test_result.hpv_status == 'positive',]
HNSC.neg <- HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='positive']

# my.gene <- c('EFNA1','EFNA3','EPHA2','EPHA3')
# my.gene <- c('TGFB1','ACVR1','TGFBR1','TGFBR2','ACVR1B')
my.gene <- c('COL4A4','SDC4','SDC1','CD44','ITGA11','ITGA1','ITGA2','ITGA3','ITGAV','ITGB1','ITGB8')
i = 1
df$gene <- t(HNSC.neg[rownames(HNSC.neg) %in% my.gene[i],])


# (2) �׷쳪���� to high expression, low expression -----------

df$quantile <- NA
df$quantile[df$gene > quantile(df$gene)[4]] <- 'high'
df$quantile[df$gene < quantile(df$gene)[2]] <- 'low'


df$mean <- ifelse(df$gene > mean(df$gene), 'high', 'low')
df$median <- ifelse(df$gene > median(df$gene), 'high', 'low')
df$patient.stage_event.clinical_stage <- matrix(ifelse(df$patient.stage_event.clinical_stage %in% 
                                                         c('stage i', 'stage ii', 'stage iii'),
                                                       'non-metastasis','metastasis'))

# (3) �����Ⱓ ----------

df$patient.months_to_death <- as.numeric(df$patient.days_to_death)/30
df <- df[!is.na(df$patient.months_to_death),]
head(df)
df <- subset(df, patient.months_to_death < 36)
# (4) censored ���� ----------

df$patient.vital_status <- factor(df$patient.vital_status,
                                  levels= c('alive','dead'),
                                  labels = c(0, 1))
df$patient.vital_status<-as.numeric(df$patient.vital_status)

# (5) survival obj ����� ----------

surv_obj <- Surv(time = df$patient.months_to_death, event = df$patient.vital_status)

# survival fit�� ���� kaplan-meier � fitting ----------

p <- survfit(formula = surv_obj ~ quantile, data = df)


# (6) �׷��� ----------

ggsurvplot(p, pval = TRUE,
           conf.int = TRUE, # �ŷڱ��� ǥ�� ����
           # risk.table = TRUE, # ���̺� ǥ�� ����
           # risk.table.height = 0.4, # ���̺� ���� ����
           ggtheme = theme_bw(), # ������ �׸� ����
           palette = c("#2E9FDF", "#FC4E07"),
           legend = 'bottom', xlim = c(0,25),
           title = my.gene[i])+
  xlab('month')


#          font.x = c(10), # x�� ���� ũ�� ����
#          font.y = c(10), # y�� ���� ũ�� ����
#          font.tickslab = c(10), # �� �� ũ�� ����
#          # surv.median.line = "hv", # 50% �������� ǥ��
#          break.time.by = 5
#         )
