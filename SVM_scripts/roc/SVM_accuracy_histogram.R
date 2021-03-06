
df <- read.csv("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM_perm_accuracy.csv")


library(ggplot2)
theme_set(theme_classic(base_size = 16))

tiff("/Users/sheilash/Desktop/projects/pfn_sex_diff/paper/figures/SVM_perm_accuracy.tiff", width = 2, height = 2.5, units = 'in', res = 300)

ggplot(df, aes(x=Accuracy)) + 
  geom_histogram(binwidth=0.005, fill ="black") +   geom_vline(aes(xintercept = as.numeric(0.83)), color = "red", linetype = "dashed") +
  scale_y_continuous(name="Count", expand = c(0, 0)) +  
  scale_x_continuous(breaks=c(0.4, 0.6, 0.8), labels=c("0.4", "0.6", "0.8")) +
  theme(axis.text.y = element_text (family="Arial", size = 10)) +
  theme(axis.title = element_text (family="Arial", size = 12, colour= "black")) +
  theme(axis.text.x = element_text (family="Arial", size = 10)) 
 
dev.off()