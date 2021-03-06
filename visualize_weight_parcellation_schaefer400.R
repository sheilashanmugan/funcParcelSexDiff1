
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)

# read in LH parc
lh_schaefer400 <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/lh.schaefer400.txt")
lh_schaefer400 <- lh_schaefer400[1:10242,]
# read in color lookup table
lh_schaefer400_ct <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/lh.schaefer400.ct.txt")
lh_schaefer400_ct_200 <- as.data.frame(lh_schaefer400_ct[-1,])
colnames(lh_schaefer400_ct_200) <- "V1"
# read in brain img data
data_brain1 <- readMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum.mat")
data_brain <-data_brain1$w.Brain.Sex.Abs.sum.lh
# parcellate data
data_vec <- as.vector(data_brain)
data_parc_lh <- cbind(data_vec, lh_schaefer400) %>% as_tibble() %>% group_by(lh_schaefer400) %>% summarise(wt_mean = mean(data_vec))

data_parc_lh_sort <- data_parc_lh[order(data_parc_lh$lh_schaefer400),]
#data_parc_lh_sort_200 <- subset(data_parc_lh_sort, lh_schaefer400 != "65793", select=-lh_schaefer400)
data_parc_lh_sort_200 <- subset(data_parc_lh_sort, lh_schaefer400 != "65793")


parc_vert<-matrix(NA,nrow=10242,ncol=1)
for (i in 1:200){
  for (j in 1:10242){
    if (data_parc_lh_sort_200$lh_schaefer400[i] == lh_schaefer400[j]){
      parc_vert[j] <- data_parc_lh_sort_200$wt_mean[i] }}}
matPath <- "/cbica/projects/funcParcelSexDiff/results/genetics/parc_mat/"
writeMat("/cbica/projects/funcParcelSexDiff/results/genetics/parc_mat/w_Brain_Sex_Abs_sum_schaefer400_lh.mat", parc_vert = parc_vert)

