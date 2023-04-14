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
# xcellnalysis ( 한 환자별, cellytpe알 수 잇음)
# 
# HNSC: TCGA data
# clinical: 환자정보
# 
# neg, pos만 
# interminate (빼고 사용)

# https://www.cancer.gov/about-cancer/diagnosis-staging/staging -> TNM 관련해서
# https://www.ncbi.nlm.nih.gov/books/NBK65719/table/CDR0000062913__892/ -> TNM 관련해서 (iva, ivb, ivc)
# https://bioinformaticsandme.tistory.com/224 -> 생존곡선 튜토리얼
# https://mustlearning.tistory.com/18 -> 생존곡선 그래프

setwd('/home/sohee/analysis/data/HNSCC/mydata/tcga/')
load('TCGAdata.RData') 

save_path = '/home/data/HNSCC/results/tcga/'

# (1) 데이터 전처리 -----------

# df <- clinical[clinical$patient.hpv_test_results.hpv_test_result.hpv_status == 'negative',]
# HNSC.neg <- HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='negative']

df <- clinical[clinical$patient.hpv_test_results.hpv_test_result.hpv_status == 'positive',]
HNSC.neg <- HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status=='positive']

# my.gene <- c('EFNA1','EFNA3','EPHA2','EPHA3')
# my.gene <- c('TGFB1','ACVR1','TGFBR1','TGFBR2','ACVR1B')
my.gene <- c('COL4A4','SDC4','SDC1','CD44','ITGA11','ITGA1','ITGA2','ITGA3','ITGAV','ITGB1','ITGB8')
i = 1
df$gene <- t(HNSC.neg[rownames(HNSC.neg) %in% my.gene[i],])


# (2) 그룹나누기 to high expression, low expression -----------

df$quantile <- NA
df$quantile[df$gene > quantile(df$gene)[4]] <- 'high'
df$quantile[df$gene < quantile(df$gene)[2]] <- 'low'


df$mean <- ifelse(df$gene > mean(df$gene), 'high', 'low')
df$median <- ifelse(df$gene > median(df$gene), 'high', 'low')
df$patient.stage_event.clinical_stage <- matrix(ifelse(df$patient.stage_event.clinical_stage %in% 
                                                         c('stage i', 'stage ii', 'stage iii'),
                                                       'non-metastasis','metastasis'))

# (3) 생존기간 ----------

df$patient.months_to_death <- as.numeric(df$patient.days_to_death)/30
df <- df[!is.na(df$patient.months_to_death),]
head(df)
df <- subset(df, patient.months_to_death < 36)
# (4) censored 여부 ----------

df$patient.vital_status <- factor(df$patient.vital_status,
                                  levels= c('alive','dead'),
                                  labels = c(0, 1))
df$patient.vital_status<-as.numeric(df$patient.vital_status)

# (5) survival obj 만들기 ----------

surv_obj <- Surv(time = df$patient.months_to_death, event = df$patient.vital_status)

# survival fit을 통해 kaplan-meier 곡선 fitting ----------

p <- survfit(formula = surv_obj ~ quantile, data = df)


# (6) 그래프 ----------

ggsurvplot(p, pval = TRUE,
           conf.int = TRUE, # 신뢰구간 표현 여부
           # risk.table = TRUE, # 테이블 표시 여부
           # risk.table.height = 0.4, # 테이블 높이 설정
           ggtheme = theme_bw(), # 데이터 테마 설정
           palette = c("#2E9FDF", "#FC4E07"),
           legend = 'bottom', xlim = c(0,25),
           title = my.gene[i])+
  xlab('month')


#          font.x = c(10), # x축 제목 크기 설정
#          font.y = c(10), # y축 제목 크기 설정
#          font.tickslab = c(10), # 축 값 크기 설정
#          # surv.median.line = "hv", # 50% 생존지점 표시
#          break.time.by = 5
#         )

