

library(R.matlab)
library(dplyr)

#PredictionFolder <- '/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/'
Behavior_data<- as.data.frame(readMat('/cbica/projects/funcParcelSexDiff/inputData/Behavior_693.mat'))
BBLID <- as.data.frame(Behavior_data$BBLID)
colnames(BBLID) <- "bblid"
Behavior_data$tertile <- ntile(Behavior_data$AgeYears, 3)
Behavior_data_young <- subset(Behavior_data, tertile == 1)
Behavior_data_middle <- subset(Behavior_data, tertile == 2)
Behavior_data_old <- subset(Behavior_data, tertile == 3)
writeMat("/cbica/projects/funcParcelSexDiff/inputData/behavior/Behavior_data_young.mat", Behavior_data_young = Behavior_data_young)
writeMat("/cbica/projects/funcParcelSexDiff/inputData/behavior/Behavior_data_middle.mat", Behavior_data_middle = Behavior_data_middle)
writeMat("/cbica/projects/funcParcelSexDiff/inputData/behavior/Behavior_data_old.mat", Behavior_data_old = Behavior_data_old)


library(Hmisc)
describe(Behavior_data_young$AgeYears)
describe(Behavior_data_middle$AgeYears)
describe(Behavior_data_old$AgeYears)

