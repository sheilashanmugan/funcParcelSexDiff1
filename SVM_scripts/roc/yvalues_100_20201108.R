
library(R.matlab)

PredictionFolder <- '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/'
Behavior_data<- as.data.frame(readMat('/cbica/projects/pncSingleFuncParcel/Replication/Revision/PredictionAnalysis/Behavior_693.mat'))
BBLID <- as.data.frame(Behavior_data$BBLID)
colnames(BBLID) <- "bblid"

data_Ymat_10 <-data.frame(matrix(NA,nrow=693,ncol=100)) 

for (i in 1:100)
{ 
  Ymat = paste0(PredictionFolder, 'Time_', i, '/Y.mat');
  data_Ymat = readMat(Ymat)
  data_Ymat0 <-data.frame(matrix(NA,nrow=length(data_Ymat$Y.group0),ncol=2))
  colnames(data_Ymat0)<- c("Num_ID", "Y_value")
  data_Ymat0[,1] <- data_Ymat$Y.group0.subj.index
  data_Ymat0[,2] <- t(as.data.frame(data_Ymat$Y.group0)) #make numeric to data frame and transpose
  data_Ymat1 <-data.frame(matrix(NA,nrow=length(data_Ymat$Y.group1),ncol=2))
  colnames(data_Ymat1)<- c("Num_ID", "Y_value")
  data_Ymat1[,1] <- data_Ymat$Y.group1.subj.index
  data_Ymat1[,2] <- t(as.data.frame(data_Ymat$Y.group1)) 
  data_Ymat_all <- rbind(data_Ymat0, data_Ymat1)
  data_Ymat_all_sort <- data_Ymat_all[order(data_Ymat_all$Num_ID),]
  data_Ymat_10[,i] <- data_Ymat_all_sort[,2]
}

data_yval_ids <- cbind(BBLID, data_Ymat_10)
data_yval_mean <- data.frame(bblid=data_yval_ids[,1], Y_means=rowMeans(data_yval_ids[,-1]))

Behavior_data$sex1m1 <- ifelse(Behavior_data$Sex==1, -1,1)
data_sex <- subset(Behavior_data, select= c("BBLID", "sex1m1"))

data_yval_sex_id <- merge(data_yval_mean, data_sex, by.x="bblid", by.y="BBLID")

mat_yval <- as.matrix(data_yval_sex_id[,2])
mat_sex <- as.matrix(data_yval_sex_id[,3])

writeMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/ROC/mat_yval.mat", mat_yval = mat_yval)
writeMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/ROC/mat_sex.mat", mat_sex = mat_sex)


write.csv(data_yval_mean, "/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/ROC/yvalues100.csv")
