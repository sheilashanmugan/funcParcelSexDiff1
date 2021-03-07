
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)
library(parallel)
source("/cbica/projects/funcParcelSexDiff/scripts/rFunctions/perm.sphere.p.R")

#brain input is from SVM run 100x
# read in LH parc
lh_schaefer1000 <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/fsaverage5/lh_Schaefer1000_17net_fs5_L.txt")
lh_schaefer1000 <- lh_schaefer1000[1:10242,]
# read in color lookup table
lh_schaefer1000_ct <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/fsaverage5/lh_schaefer1000_17net_fs5_ct.txt")
lh_schaefer1000_ct_500 <- as.data.frame(lh_schaefer1000_ct[-1,])
colnames(lh_schaefer1000_ct_500) <- "V1"
# read in brain img data
data_brain1 <- readMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum.mat")
data_brain <-data_brain1$w.Brain.Sex.Abs.sum.lh
# parcellate data
data_vec <- as.vector(data_brain)
data_parc_lh <- cbind(data_vec, lh_schaefer1000) %>% as_tibble() %>% group_by(lh_schaefer1000) %>% summarise(wt_mean = mean(data_vec))

data_parc_lh_sort <- data_parc_lh[order(data_parc_lh$lh_schaefer1000),]
data_parc_lh_sort_500 <- subset(data_parc_lh_sort, lh_schaefer1000 != "65793", select=-lh_schaefer1000)





###read in gene data
geneinfo <- readMat("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/AHBAprocessed/ROIxGene_Schaefer1000_INT.mat")
geneinfo1 <- lapply(geneinfo, unlist, use.names=FALSE)
geneMat <- as.data.frame(geneinfo1$parcelExpression)
geneMat1 <- subset(geneMat, select = -V1)
gene <- read.csv("/Users/sheilash/Desktop/projects/pfn_sex_diff/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/schaefer1000/GeneSymbol.csv",header=TRUE)
colnames(gene) <- "V1"
colnames(geneMat1) <- gene$V1
geneMat2 <- cbind(lh_schaefer1000_ct_500, geneMat1)
geneMat2_sort <- geneMat2[order(geneMat2$V1),]
geneMat2_sort_500 <- subset(geneMat2_sort, select=-V1)



#read in permutations of schaefer left hemisphere (from Jakob)
schaefer_perm <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/permutations/Schaefer1000_perm.csv", header=TRUE)
#schaefer_perm1k <- as.data.frame(schaefer_perm[,1:1000])

#Select genes for chromosome enrichement 
chrom <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_500 <-geneMat2_sort_500[,match(geneRefined,colnames(geneMat2_sort_500),nomatch = 0)] #target gene first. where in #2 is #1

#make vector for perm,sphere.p
data_parc_lh_sort_500_vec <- as.vector(data_parc_lh_sort_500$wt_mean)

set.seed(1)
img_gene_permSphereP <- mclapply(geneMatChrom_500, FUN=perm.sphere.p, y=data_parc_lh_sort_500_vec, perm.id= schaefer_perm, corr.type='pearson', mc.set.seed =  TRUE)
img_gene_permSphereP <- mclapply(geneMatChrom_500, FUN=perm.sphere.p, y=data_parc_lh_sort_500_vec, perm.id= schaefer_perm, corr.type='pearson', mc.set.seed =  TRUE, mc.cores=getOption("mc.cores", 20))


write.csv(img_gene_permSphereP, "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img100_schaefer10K_gene_permSphereP_df_20200727.csv")

