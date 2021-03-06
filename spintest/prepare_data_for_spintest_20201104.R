
library(R.matlab)

data_brainSVM <- readMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum.mat")
data_brainSVM_lh <-data.frame(t(data_brainSVM$w.Brain.Sex.Abs.sum.lh))
data_brainSVM_rh <-data.frame(t(data_brainSVM$w.Brain.Sex.Abs.sum.rh))


data_brainGAM <- readMat("/cbica/projects/funcParcelSexDiff/results//GamAnalysis/AtlasLoading/SexEffects/UnthreshAbsSum/Gam_Sex_Abs_sum.mat")
data_brainGAM_lh <-data.frame(t(data_brainGAM$Gam.Sex.Abs.sum.lh))
data_brainGAM_rh <-data.frame(t(data_brainGAM$Gam.Sex.Abs.sum.rh))

write.table(data_brainSVM_lh, "/cbica/projects/funcParcelSexDiff/inputData/spintest/wBrainSexAbssum_lh.csv", row.names = FALSE, col.names =FALSE)
write.table(data_brainSVM_rh, "/cbica/projects/funcParcelSexDiff/inputData/spintest/wBrainSexAbssum_rh.csv", row.names = FALSE, col.names =FALSE)
write.table(data_brainGAM_lh, "/cbica/projects/funcParcelSexDiff/inputData/spintest/GamSexAbssum_lh.csv", row.names = FALSE, col.names =FALSE)
write.table(data_brainGAM_rh, "/cbica/projects/funcParcelSexDiff/inputData/spintest/GamSexAbssum_rh.csv", row.names = FALSE, col.names =FALSE)
