library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)
source("/cbica/projects/funcParcelSexDiff/scripts/rFunctions/perm.sphere.p.R")


# read in LH parc
lh_schaefer400 <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/lh.schaefer400.txt")
lh_schaefer400 <- lh_schaefer400[1:10242,]
# read in color lookup table
lh_schaefer400_ct <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/lh.schaefer400.ct.txt")
lh_schaefer400_ct_200 <- as.data.frame(lh_schaefer400_ct[-1,])
colnames(lh_schaefer400_ct_200) <- "V1"
# read in brain img data
data_brain1 <- readMat("/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov/w_Brain_Sex_Abs_sum.mat")
data_brain <-data_brain1$w.Brain.Sex.Abs.sum.lh
# parcellate data
data_vec <- as.vector(data_brain)
data_parc_lh <- cbind(data_vec, lh_schaefer400) %>% as_tibble() %>% group_by(lh_schaefer400) %>% summarise(wt_mean = mean(data_vec))

data_parc_lh_sort <- data_parc_lh[order(data_parc_lh$lh_schaefer400),]
data_parc_lh_sort_200 <- subset(data_parc_lh_sort, lh_schaefer400 != "65793", select=-lh_schaefer400)



###read in gene data
gene <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/genes_20647_schaefer_V2.csv",header=FALSE)
geneMat <- readMat("/cbica/projects/funcParcelSexDiff/inputData/genetics/gene_regional_expression_zscored_Schaefer_V2.mat")
geneMat1<- as.data.frame(geneMat$gene.regional.expression)
colnames(geneMat1) <- gene$V1
geneMat2 <- cbind(lh_schaefer400_ct_200, geneMat1)
geneMat2_sort <- geneMat2[order(geneMat2$V1),]
geneMat2_sort_200 <- subset(geneMat2_sort, select=-V1)

#read in permutations of schaefer left hemisphere (from Jakob)
schaefer_perm <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/permutations/Schaefer400_perm.csv", header=TRUE)

#Select genes for chromosome enrichement 
chrom <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_200 <-geneMat2_sort_200[,match(geneRefined,colnames(geneMat2_sort_200),nomatch = 0)] #target gene first. where in #2 is #1


#make vector for perm,sphere.p
data_parc_lh_sort_200_vec <- as.vector(data_parc_lh_sort_200$wt_mean)


geneMatChrom_200_sub <- as.data.frame(geneMatChrom_200[,1:4])
img_gene_permSphereP <- lapply(geneMatChrom_200_sub, FUN=perm.sphere.p,  y=data_parc_lh_sort_200_vec, perm.id= schaefer_perm,  corr.type='pearson')


#create outpt df
img_gene_permSphereP_df <-data.frame(matrix(NA,nrow=20647,ncol=4)) 
namelist <- c("gene_name", "r_value", "p_unadj", "p_fdr")
colnames(img_gene_permSphereP_df)<- namelist

#correlate sorted img and gene data (medial wall value removed in both)
for (i in 1:15036) 
{ 
  img_gene_permSphereP_df[i,1] <- colnames(geneMatChrom_200[i])
  geneMat2_col <- as.vector(geneMatChrom_200[,i])
  
  gene_r <- cor(data_parc_lh_sort_200_vec,geneMat2_col)
  img_gene_permSphereP_df[i,2] <- gene_r
  
  perm_pvalue <- perm.sphere.p(data_parc_lh_sort_200_vec, geneMat2_col, schaefer_perm, corr.type='pearson')
  img_gene_permSphereP_df[i,3] <- perm_pvalue
  p_fdr <- p.adjust(perm_pvalue, method="fdr", n=15036)
  img_gene_permSphereP_df[i,4] <- p_fdr
}


write.csv(img_gene_permSphereP_df, "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img_gene_permSphereP_df_20200614.csv")
