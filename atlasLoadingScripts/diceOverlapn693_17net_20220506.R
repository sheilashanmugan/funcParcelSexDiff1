
library(ggplot2)
library(R.matlab)
library(dplyr)
library(RNifti)
library(data.table)
theme_set(theme_classic(base_size = 16))


#get subject list
Behavior_data<- as.data.frame(readMat('/cbica/projects/funcParcelSexDiff/inputData/Behavior_693.mat'))
bblidList <- as.list(Behavior_data$BBLID)

#directories with single subject atlas hard labels
mask1dir <- "/cbica/projects/funcParcelSexDiff/data/Revision/SingleParcellation/SingleAtlas_Analysis/Atlas_Visualize/Group_AtlasLabel_Network_"
mask2dir <- "/cbica/projects/funcParcelSexDiff/data/Revision/SingleParcellation/SingleAtlas_Analysis/Atlas_Visualize"



#calculate dice for all networks for each subjects
for(j in 1:17){
  print(j);
  dicepath <- paste0("/cbica/projects/funcParcelSexDiff/Revision/inputData/diceResults/dice_results_", j, ".csv")
  if (!file.exists(dicepath)){
    dice.results<-data.frame(matrix(NA,nrow=length(bblidList),ncol=2))
    colnames(dice.results)<- c("BBLID", paste0("net",j))
    for(i in 1:length(bblidList)){
      mask1path <- paste0(mask1dir, j, ".dlabel.nii");
      vol1 <- readNifti(mask1path)
      mask2path <- paste0(mask2dir, '/', as.character(bblidList[i]), '/', as.character(bblidList[i]), "_AtlasLabel_17Colors_Network_", j, ".dlabel.nii");
      vol2 <- readNifti(mask2path)
      dicenet <- 2 * sum(vol1*vol2) / (sum(vol1) + sum(vol2))
      dice.results[i, 1] <- bblidList[i]
      dice.results[i, 2] <- dicenet
      
    }
    write.csv(dice.results, dicepath)
  }
}



#create df to read in dice results
netList <- c("1", "2", "3","4", "5", "6","7", "8", "9","10", "11", "12","13", "14", "15","16", "17")
dice.results<-data.frame(matrix(NA,nrow=length(bblidList),ncol=17))
colnames(dice.results)<- netList
rownames(dice.results)<- bblidList

#read in dice results
for(i in 1:17){
  print(i);
  dicepath <- paste0("/cbica/projects/funcParcelSexDiff/Revision/inputData/diceResults/dice_results_", i, ".csv")
  tmp1 <- read.csv(dicepath)
  dice.results[i] <- tmp1[3]
}

#calculate medians for each network
setDT(dice.results, keep.rownames = "BBLID")
dice_medians<- as.data.frame(apply(dice.results[,2:18],2,median))
colnames(dice_medians) <- "dice_median"
setDT(dice_medians, keep.rownames = "network")
dice_medians$netrank <- rank(dice_medians$dice_median)
dice_medians <- subset(dice_medians, select = -c(dice_median))




#create color labels for the networks
network <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
             13, 14, 15, 16, 17)

netColor <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
              "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
              "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
              "#4E31A8", "#F5BA2E")


netColorDf <- cbind(network, netColor)


#combine df of medians with color table
dice_medians_col <- merge(dice_medians, netColorDf, by = "network")

mdata <- melt(dice.results, id=c("BBLID"))
colnames(mdata) <- c("BBLID", "network", "dice")
mdata_col <- merge(mdata, dice_medians_col, by="network")


#make color label an ordered factor so ggplot can match the color to network
mdata_col$netColorF <- factor(mdata_col$netColor, ordered =T)

mdata_col <- mdata_col %>% 
  mutate(network = as.factor(network),netColorF=as.character(netColorF))%>%
  mutate(network = reorder(network,netrank))%>%
  arrange (netrank)
colormap=mdata_col %>% select(network,netColorF)%>%unique()

tiff("/Users/sheilash/Desktop/projects/pfn_sex_diff/paper/figures/revision/violin_plot_dice.tiff", width = 8, height = 2.5, units = 'in', res = 300)
ggplot(mdata_col, aes(x=network, y=dice, fill=network, color=network)) +
  geom_violin(trim = FALSE) +
  xlab("Network") + ylab("Dice Coefficient") +
  scale_fill_manual(values = colormap$netColorF) +
  scale_color_manual(values = colormap$netColorF) +
  theme(legend.position="none") + theme(axis.title = element_text(size = 14), axis.text.x = element_text(size= 12), axis.text.y = element_text(size= 12))
dev.off()

