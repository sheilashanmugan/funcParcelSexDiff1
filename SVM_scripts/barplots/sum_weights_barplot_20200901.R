
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)

#read in brain data
data_brain1 <- readMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Matrix.mat")
data_brain <-data_brain1$w.Brain.Sex.Matrix

#function to sum negative weights
sumneg <- function(x) sum(x[x<0])
sums_neg<- as.data.frame(apply(data_brain,1,sumneg))

#function to sum positive weights
sumpos <- function(x) sum(x[x>0])
sums_pos<- as.data.frame(apply(data_brain,1,sumpos))

#name columns, make row number the network label, and combine positive and negative sums into 1 df
colnames(sums_neg) <- "weights"
setDT(sums_neg, keep.rownames = "network")
colnames(sums_pos) <- "weights"
setDT(sums_pos, keep.rownames = "network")
sums_all <- rbind(sums_pos, sums_neg)


#function to get sum of absolute value of weights and put in df w corresponding network
sumabs <- function(x) sum(abs(x))
sums_abs<- as.data.frame(apply(data_brain,1,sumabs))
colnames(sums_abs) <- "weights"
setDT(sums_abs, keep.rownames = "network")


#create color labels for the networks
network <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
  13, 14, 15, 16, 17)

netColor <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
                 "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
                 "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
                 "#4E31A8", "#F5BA2E")


netColorDf <- cbind(network, netColor)

#combine df of sum of abs value of weights with color table
sums_abs_col <- merge(sums_abs, netColorDf, by = "network")

#make color label an ordered factor so ggplot can match the color to network
sums_abs_col$netColorF <- factor(sums_abs_col$netColor, ordered =T)

#plot of sum of abs val of weights
ggplot(sums_abs_col, aes(x = reorder(network, weights), y = weights))+
  geom_bar(stat = "identity", position = "identity", fill= sums_abs_col$netColorF) + xlab("Network") + ylab("Sum of Weights") + geom_hline(yintercept = 0)



#get network ranking based on absolute value so other plots can have networks in the same order
sums_abs_col$netrank <- rank(sums_abs_col$weights)
sum_cols <- subset(sums_abs_col, select = -weights)

#merge df of rankings with df of pos +neg sums
sums_all_col_rank <- merge(sums_all, sum_cols, by = "network")

#plot of pos and neg sums
ggplot(sums_all_col_rank, aes(x = reorder(network, netrank), y = weights))+
  geom_bar(stat = "identity", position = "identity", fill= sums_all_col_rank$netColorF) + xlab("Network") + ylab("Sum of Weights") + geom_hline(yintercept = 0)



