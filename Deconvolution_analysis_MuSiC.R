# Tutorial: https://xuranw.github.io/MuSiC/articles/MuSiC.html
# (1) 클러스터와 비슷한 세포타입을 그룹화함.
# (2) 클러스터 비율을 추정.
# (3) 각 클러스터에서 이 과정을 유사한 방식(recursively) 으로 반복

# < Pre-determinination of cell type > ----------
# scRNAseq 데이터에서 recursively(재귀적으로) TCGA의 세포타입을 추정함  
# (1) 높은 레벨의 클러스터의 비율 추정
# (2) 각 클러스터 내의 세포 유형 추정


# < load packages > ----------
.libPaths('/home/sohee/myProject/Library/')
library(MuSiC) 
library(Seurat)
library(SingleCellExperiment)
library(Biobase)
library(reshape)
library(ggplot2)

# < Data preparation > ----------
# Bulk data
# Download Mouse bulk dataset from Github
load('/home/sohee/analysis/data/HNSCC/mydata/tcga/TCGAdata.RData')
Mouse.bulk.eset$sampleID
head(Mouse.bulk.eset@phenoData@data$sampleID)
rownames(Mouse.bulk.eset@featureData@data)


my.assay <- (HNSC[,clinical$patient.hpv_test_results.hpv_test_result.hpv_status %in% c('negative','positive')])
my.assay <- as.matrix(my.assay)
dim(my.assay) 

clinical <- clinical[clinical$patient.hpv_test_results.hpv_test_result.hpv_status %in% c('negative','positive'),]
my.pheno <- clinical[,'patient.hpv_test_results.hpv_test_result.hpv_status']
my.pheno <- data.frame(my.pheno)
rownames(my.pheno) <- colnames(my.assay)
colnames(my.pheno) <- 'HPV'

# class(my.pheno)
my.pheno <- new("AnnotatedDataFrame", data = my.pheno)
my.bulk.eset <- ExpressionSet(assayData = my.assay, phenoData = my.pheno)
my.bulk.eset

# Download Mouse single cell dataset from Github
load('/home/sohee/analysis/data/HNSCC/mydata/seurat/v3/HNSCC_subcelltype.RData')
my.sce <- SingleCellExperiment(assays = list(counts=HNSCC@assays$RNA@counts),
                               colData = HNSCC@meta.data)
my.sce
my.sce$cellType <- my.sce$sub.celltype
levels(my.sce$cellType)


# < Estimation of cell type proportions > ----------
# Produce the first step information
my.basis = music_basis(my.sce, 
                       clusters = 'cellType', 
                       samples = 'HPV') 
                       # select.ct = c('Endo', 'Podo', 'PT', 'LOH', 'DCT', 'CD-PC', 'CD-IC', 'Fib','Macro', 'Neutro','B lymph', 'T lymph', 'NK'))

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(my.basis$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(my.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')

# We manually specify the cluster and annotated single cell data with cluster information.
clusters.type = list(
  C1 = c('Pericyte'),
  C2 = c('Mast'),
  C3 = c('Dendritic','Plasma','Macrophage','ImmatureT','B','NK','CD8T','Treg','CD4T'),
  C4 = c('Tumor','Normal'),
  C5 = c('Endothelial','Myofibroblast','NAF','CAF')
)

cl.type = as.character(my.sce$cellType)

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
my.sce$clusterType = factor(cl.type, levels = c(names(clusters.type))) #, 'CD-Trans', 'Novel1', 'Novel2'))

# 13 selected cell types
s.my = unlist(clusters.type)
s.my



load('/home/sohee/analysis/data/HNSCC/mydata/tcga/IEmarkers.RData')
# This RData file provides two vectors of gene names Epith.marker and Immune.marker

# We now construct the list of group marker
C4 <- c('EPCAM', 'KRT14')
C5 <- c('COL1A1','COL1A2')
IEmarkers = list(NULL, NULL, NULL, C4, C5) # Epith.marker, Immune.marker)
names(IEmarkers) = c('C1', 'C2', 'C3', 'C4', 'C5')
# The name of group markers should be the same as the cluster names

Est.my.bulk = music_prop.cluster(bulk.mtx = exprs(my.bulk.eset),
                                    sc.sce = my.sce,
                                    group.markers = IEmarkers, 
                                    clusters = 'cellType',
                                    groups = 'clusterType', 
                                    samples = 'HPV', 
                                    clusters.type = clusters.type)
Est.y.bulk






# < Deconvolution with TCGA > ----------

# This RData file provides two vectors of gene names Epith.marker and Immune.marker


clusters.type = list(C1 = 'Neutro', C2 = 'Podo', C3 = c('Endo', 'CD-PC', 'LOH', 'CD-IC', 'DCT', 'PT'), C4 = c('Macro', 'Fib', 'B lymph', 'NK', 'T lymph'))

cl.type = as.character(Mousesub.sce$cellType)

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
Mousesub.sce$clusterType = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))

# 13 selected cell types
s.mouse = unlist(clusters.type)
s.mouse


# We now construct the list of group marker
IEmarkers = list(NULL, NULL, NULL, NULL, NULL)
names(IEmarkers) = c('C1', 'C2', 'C3', 'C4', 'C5')




# < Deconvolution to TCGA > ----------
Est.my.bulk = music_prop.cluster(bulk.mtx = exprs(my.bulk), 
                                    sc.sce = my.sce, 
                                    group.markers = IEmarkers, 
                                    clusters = 'clusterType', 
                                    groups = 'sub.celltype',
                                    samples = 'HPV', 
                                    clusters.type = clusters.type)

