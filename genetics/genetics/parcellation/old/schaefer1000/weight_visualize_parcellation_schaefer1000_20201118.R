
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)

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
#data_parc_lh_sort_500 <- subset(data_parc_lh_sort, lh_schaefer1000 != "65793", select=-lh_schaefer1000)
data_parc_lh_sort_500 <- subset(data_parc_lh_sort, lh_schaefer1000 != "65793")


parc_vert<-matrix(NA,nrow=10242,ncol=1)
for (i in 1:500){
  for (j in 1:10242){
    if (data_parc_lh_sort_500$lh_schaefer1000[i] == lh_schaefer1000[j]){
      parc_vert[j] <- data_parc_lh_sort_500$wt_mean[i] }}}
matPath <- "/cbica/projects/funcParcelSexDiff/results/genetics/parc_mat/"
writeMat("/cbica/projects/funcParcelSexDiff/results/genetics/parc_mat/w_Brain_Sex_Abs_sum_schaefer1000_lh.mat", parc_vert = parc_vert)

