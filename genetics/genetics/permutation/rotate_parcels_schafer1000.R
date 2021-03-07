
library(matrixStats)
source("/cbica/projects/funcParcelSexDiff/scripts/rFunctions/rotate.parcellation.R")

rois <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mmCentroid_RAS.csv")
coord.l <- data.matrix(rois[1:500,3:5])
coord.r <- data.matrix(rois[501:1000,3:5])

schaefer1000_perm <- rotate.parcellation(coord.l,coord.r, nrot=10000)

schaefer1000_perm_lh <- schaefer1000_perm[1:500,]

write.csv(schaefer1000_perm_lh, "/cbica/projects/funcParcelSexDiff/inputData/genetics/permutations/Schaefer1000_perm.csv")
