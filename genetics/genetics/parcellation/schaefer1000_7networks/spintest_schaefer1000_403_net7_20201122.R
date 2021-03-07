#brain input is from SVM run 100x
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)
library(parallel)
source("/cbica/projects/funcParcelSexDiff/scripts/rFunctions/perm.sphere.p.R")


geneMatChrom_500_comp_fin <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/geneMatChrom_500_comp_fin_7net.csv", header=TRUE)
data_parc_lh_sort_500_comp <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/data_parc_lh_sort_500_comp_7net.csv", header=TRUE)
img_geneMatChrom_roi_500_comp <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img_geneMatChrom_roi_500_comp_7net.csv", header=TRUE)

geneMatChrom_500_comp_fin <- geneMatChrom_500_comp_fin[,-1]
img_geneMatChrom_roi_500_comp <- img_geneMatChrom_roi_500_comp[,-1]

schaefer1000_perm_lh1 <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/permutations/Schaefer1000_perm_403_500_7net.csv", header = TRUE)
schaefer1000_perm_lh <- schaefer1000_perm_lh1[,-1]

#make vector for perm,sphere.p
data_parc_lh_sort_500_comp_vec <- as.vector(data_parc_lh_sort_500_comp$wt_mean)


set.seed(1)
img_gene_permSphereP <- mclapply(geneMatChrom_500_comp_fin, FUN=perm.sphere.p, y=data_parc_lh_sort_500_comp_vec, perm.id= schaefer1000_perm_lh, corr.type='pearson', mc.set.seed =  TRUE, mc.cores=getOption("mc.cores", 30))


write.csv(img_gene_permSphereP, "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img100_schaefer1000__7net_perm10K_SortParcNum_gene_permSphereP_df_20201122.csv") 
