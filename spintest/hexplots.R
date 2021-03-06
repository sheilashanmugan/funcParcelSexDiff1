library(R.matlab)
library(ggplot2)
library(hexbin)
theme_set(theme_classic(base_size = 16))

data_brainSVM <- readMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum_NoMedialWall.mat")
data_brainSVM_lh <-data.frame(t(data_brainSVM$w.Brain.Sex.Abs.sum.lh.NoMedialWall))
data_brainSVM_rh <-data.frame(t(data_brainSVM$w.Brain.Sex.Abs.sum.rh.NoMedialWall))
data_brainSVM_All <-data.frame(t(data_brainSVM$w.Brain.Sex.Abs.sum.All.NoMedialWall))
colnames(data_brainSVM_lh) <- "weights"
colnames(data_brainSVM_rh) <- "weights"
colnames(data_brainSVM_All) <- "weights"

data_brainGAM <- readMat("/cbica/projects/funcParcelSexDiff/results//GamAnalysis/AtlasLoading/SexEffects/UnthreshAbsSum/Gam_Sex_Abs_sum_NoMedialWall.mat")
data_brainGAM_lh <-data.frame(t(data_brainGAM$Gam.Sex.Abs.sum.lh.NoMedialWall))
data_brainGAM_rh <-data.frame(t(data_brainGAM$Gam.Sex.Abs.sum.rh.NoMedialWall))
data_brainGAM_All <-data.frame(t(data_brainGAM$Gam.Sex.Abs.sum.All.NoMedialWall))

colnames(data_brainGAM_lh) <- "loadings"
colnames(data_brainGAM_rh) <- "loadings"
colnames(data_brainGAM_All) <- "loadings"


cor.test(data_brainSVM_rh$weights, data_brainGAM_rh$loadings, method = "pearson")


myPalette <- c("#333333", "#4C4C4C", "#666666", "#7F7F7F", "#999999", "#B2B2B2", "#CCCCCC");
# rh
Data_tmp1 = data.frame(data_brainGAM_rh_loadings = as.numeric(data_brainGAM_rh$loadings));
Data_tmp1$data_brainSVM_rh_weights = as.numeric(data_brainSVM_rh$weights);
cor.test(Data_tmp1$data_brainGAM_rh_loadings, Data_tmp1$data_brainSVM_rh_weights, method = "pearson")

hexinfo <- hexbin(Data_tmp1$data_brainSVM_rh_weights, Data_tmp1$data_brainGAM_rh_loadings, xbins = 30);
data_hex <- data.frame(hcell2xy(hexinfo), count = hexinfo@count);

ggplot() +
  geom_hex(data = subset(data_hex, count > 10), aes(x, y, fill = count), stat = "identity") +
  scale_fill_gradientn(colours = myPalette, breaks = c(50, 100, 150)) + 
  geom_smooth(data = Data_tmp1, aes(x = data_brainSVM_rh_weights, y = data_brainGAM_rh_loadings), method = lm, color = "#FFFFFF", linetype = "dashed") +
  labs(x = "SVM Weights", y = "GAM Loadings") + 
  theme(axis.text=element_text(family="Arial", size=16, color='black'), axis.title=element_text(family="Arial", size=16), aspect.ratio = 1) +
  theme(legend.text = element_text(family="Arial", size = 16), legend.title = element_text(family="Arial", size = 16)) +
  theme(legend.justification = c(1, 1), legend.position = c(1.02, 1)) 

theme(axis.text.y = element_text (family="Arial", size = 16, colour= "black", hjust=0.5)) +
  

# All
Data_tmp1 = data.frame(data_brainGAM_All_loadings = as.numeric(data_brainGAM_All$loadings));
Data_tmp1$data_brainSVM_All_weights = as.numeric(data_brainSVM_All$weights);
cor.test(Data_tmp1$data_brainGAM_All_loadings, Data_tmp1$data_brainSVM_All_weights, method = "pearson")

hexinfo <- hexbin(Data_tmp1$data_brainSVM_All_weights, Data_tmp1$data_brainGAM_All_loadings, xbins = 30);
data_hex <- data.frame(hcell2xy(hexinfo), count = hexinfo@count);

tiff("/Users/sheilash/Desktop/projects/pfn_sex_diff/paper/figures/hexplot_all.tiff", width = 3.5, height = 3.5, units = 'in', res = 300)
ggplot() +
  geom_hex(data = subset(data_hex, count > 10), aes(x, y, fill = count), stat = "identity") +
  scale_fill_gradientn(colours = myPalette, breaks = c(100, 200, 300), name = "Count") + 
  geom_smooth(data = Data_tmp1, aes(x = data_brainSVM_All_weights, y = data_brainGAM_All_loadings), method = lm, color = "#FFFFFF", linetype = "dashed") +
  theme_classic() + labs(x = "SVM Weights", y = "GAM Loadings") + 
  theme(axis.text=element_text(family="Arial", size=16, color='black'), axis.title=element_text(family="Arial", size=16), aspect.ratio = 1) +
  theme(legend.text = element_text(family="Arial", size = 10), legend.title = element_text(family="Arial", size = 16)) +
  theme(legend.justification = c(1, 1), legend.position = c(1.02, 1)) + theme(axis.text.x = element_text(size= 10), axis.text.y = element_text(size= 10))
dev.off()
# 
# +
#   scale_x_continuous(limits = c(0, 0.04), breaks = c( 0, 0.01, .02, .03, 0.04)) + 
#   scale_y_continuous(limits = c(0, 26.5), breaks = c(0, 5, 10, 15, 20, 25));


