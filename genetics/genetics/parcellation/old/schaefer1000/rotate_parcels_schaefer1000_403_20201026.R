

library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(data.table)
library(parallel)
source("/cbica/projects/funcParcelSexDiff/scripts/rFunctions/perm.sphere.p.R")
library(matrixStats)
source("/cbica/projects/funcParcelSexDiff/scripts/rFunctions/rotate.parcellation.R")


geneMatChrom_500_comp_fin <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/geneMatChrom_500_comp_fin.csv", header=TRUE)
data_parc_lh_sort_500_comp <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/data_parc_lh_sort_500_comp.csv", header=TRUE)
img_geneMatChrom_roi_500_comp <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img_geneMatChrom_roi_500_comp.csv", header=TRUE)

geneMatChrom_500_comp_fin <- geneMatChrom_500_comp_fin[,-1]
data_parc_lh_sort_500_comp <- data_parc_lh_sort_500_comp[,-1]
img_geneMatChrom_roi_500_comp <- img_geneMatChrom_roi_500_comp[,-1]

coord.l <- data.matrix(img_geneMatChrom_roi_500_comp[1:length(img_geneMatChrom_roi_500_comp$R),3:5])
rois <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mmCentroid_RAS.csv")
coord.r <- data.matrix(rois[501:1000,3:5])

schaefer1000_perm <- rotate.parcellation(coord.l,coord.r, nrot=10000)
schaefer1000_perm_lh <- schaefer1000_perm[1:length(img_geneMatChrom_roi_500_comp$R),]
write.csv(schaefer1000_perm_lh, "/cbica/projects/funcParcelSexDiff/inputData/genetics/permutations/Schaefer1000_perm_403_500.csv")

#make vector for perm,sphere.p
data_parc_lh_sort_500_comp_vec <- as.vector(data_parc_lh_sort_500_comp$wt_mean)

set.seed(1)
img_gene_permSphereP <- mclapply(geneMatChrom_500_comp_fin, FUN=perm.sphere.p, y=data_parc_lh_sort_500_comp_vec, perm.id= schaefer1000_perm_lh, corr.type='pearson', mc.set.seed =  TRUE, mc.cores=getOption("mc.cores", 20))


write.csv(img_gene_permSphereP, "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img100_schaefer1000_perm10K_SortParcNum_gene_permSphereP_df_20201023.csv")
