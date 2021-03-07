
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)
library(ggplot2)
library(hexbin)
library(matrixStats)
library(plotrix)
library(scales)
theme_set(theme_classic(base_size = 16))

#brain input is from SVM run 100x
# read in LH parc
lh_schaefer1000 <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/lh_Schaefer1000_7net_L.csv")
lh_schaefer1000 <- lh_schaefer1000[1:10242,]
# read in color lookup table
lh_schaefer1000_ct <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/lh_Schaefer1000_7net_cttable.csv")
lh_schaefer1000_ct_500 <- as.data.frame(lh_schaefer1000_ct[-1,])
colnames(lh_schaefer1000_ct_500) <- "V1"
setDT(lh_schaefer1000_ct_500, keep.rownames = "parcNum")

# read in brain img data
data_brain1 <- readMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum.mat")
data_brain <-data_brain1$w.Brain.Sex.Abs.sum.lh
# parcellate data
data_vec <- as.vector(data_brain)
data_parc_lh <- cbind(data_vec, lh_schaefer1000) %>% as_tibble() %>% group_by(lh_schaefer1000) %>% summarise(wt_mean = mean(data_vec))

#order by parcel number for rotations
data_parc_lh <- subset(data_parc_lh, lh_schaefer1000 != "65793")
data_parc_lh_parcnum <- data_parc_lh
data_parc_lh_parcnum$parcNum <- "NA"

for (i in 1:500){
  for (j in 1:500){
    if (data_parc_lh_parcnum$lh_schaefer1000[i] == lh_schaefer1000_ct_500$V1[j]){
      data_parc_lh_parcnum$parcNum[i] <- lh_schaefer1000_ct_500$parcNum[j] }}}

data_parc_lh_parcnum$parcNum <- as.numeric(data_parc_lh_parcnum$parcNum)
data_parc_lh_parcnum_sort <- data_parc_lh_parcnum[order(data_parc_lh_parcnum$parcNum),]


###read in gene data
geneinfo <- readMat("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/AHBAprocessed/ROIxGene_Schaefer1000_INT.mat")
geneinfo1 <- lapply(geneinfo, unlist, use.names=FALSE)
geneMat <- as.data.frame(geneinfo1$parcelExpression)
geneMat1 <- subset(geneMat, select = -V1)
gene <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/schaefer1000/GeneSymbol.csv",header=TRUE)
colnames(gene) <- "V1"
colnames(geneMat1) <- gene$V1
geneMat2 <- cbind(lh_schaefer1000_ct_500, geneMat1)
geneMat2$parcNum <- as.numeric(geneMat2$parcNum)
geneMat2_sort <- geneMat2[order(geneMat2$parcNum),]
geneMat2_sort_500 <- subset(geneMat2_sort, select=-V1)


#merge img, gene, and centroid coord. remove parcels with missing gene data
img_geneMatChrom_500 <- merge(data_parc_lh_parcnum_sort,geneMat2_sort_500, by="parcNum")


rois <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mmCentroid_RAS.csv")
roislh <- rois[1:500,]

img_geneMatChrom_roi_500 <- merge(roislh, img_geneMatChrom_500, by.x = "RoiLabel", by.y= "parcNum",all = FALSE)
img_geneMatChrom_roi_500_comp <- subset(img_geneMatChrom_roi_500, !is.nan(img_geneMatChrom_roi_500[[10]])) ##if add more columns, change "10". This is the complete DF. Any subsets below should come from here/match this order.


data_parc_lh_sort_500_comp <- subset(img_geneMatChrom_roi_500_comp, select=wt_mean)


othercols <- c(colnames(rois), "lh_schaefer1000", "wt_mean")
geneMatChrom_500_comp <- img_geneMatChrom_roi_500_comp[(length(othercols)+1):(length(colnames(img_geneMatChrom_roi_500_comp)))]



# 
# #correlate sorted img and gene data (medial wall value removed in both) for other enrichement analyses
# rankedGeneList <- sapply(1:15745, function(x) cor(data_parc_lh_sort_500_comp,geneMatChrom_500_comp[,x]))
# rankedGeneList1<- cbind(gene, rankedGeneList)
# rankedGeneList2 <- rankedGeneList1[with(rankedGeneList1, order(-rankedGeneList)),]
# rankedGeneList3 <- subset(rankedGeneList2, select = V1)
# write.csv(rankedGeneList3, "/cbica/projects/funcParcelSexDiff/results/genetics/GO/RankedGeneList_wSex100_Cor_schaefer1000Net7_20201123.csv", row.names = FALSE,  quote=FALSE) 
# 


#Select genes for chromosome enrichement 
chrom <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_500_comp_fin <-geneMatChrom_500_comp[,match(geneRefined,colnames(geneMatChrom_500_comp),nomatch = 0)] #target gene first. where in #2 is #1




geneRefined_df <- as.data.frame(geneRefined)
chromRefined_df <- as.data.frame(chromRefined)
chromGeneRefined <- cbind(chromRefined_df, geneRefined_df)






pval1 <- read.csv( "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img100_schaefer1000__7net_perm10K_SortParcNum_gene_permSphereP_df_20201122.csv")
pval <- pval1[,-1]
df2 <- data.frame(t(pval))
colnames(df2) <- "pval"
df2$fdrp <- p.adjust(df2$pval, method="fdr", n=12986)
setDT(df2, keep.rownames = "gene")
chromGeneRefined$geneRefined <- gsub("-", "." , chromGeneRefined$geneRefined)

df3a <- merge(df2, chromGeneRefined, by.x="gene", by.y="geneRefined", all=TRUE)

t2 <- sapply(1:12986, function(x) cor(data_parc_lh_sort_500_comp,geneMatChrom_500_comp_fin[,x]))
colt2<- as.vector(chromGeneRefined$geneRefined)
dft2 <- as.data.frame(cbind(colt2, t2))
colnames(dft2)<- c("gene", "corr")

df3 <- merge(dft2, df3a, by="gene")



dfsig <- subset(df3, fdrp <= 0.05)
dfsig$corrSigRanked <- rank(dfsig$corr)
sigt3 <- sapply(levels(dfsig$chromRefined)[-1], function(x) median(dfsig$corrSigRanked[which(dfsig$chromRefined==x)],na.rm = TRUE)) # which(chromRefined==x) --  get an index of all the trues,,,   indexing out of t2Ranked where chromRefined==x.... find the median of t2Ranked for each of these levels
sigt4<- sigt3 - median(dfsig$corrSigRanked) # for the average correlation ranking for each chromosome, how far away it is from the median
sigt4
sigt5 <- as.data.frame(sigt4)
setDT(sigt5, keep.rownames = "chromosome")
sigt5$chromosome <- as.factor(sigt5$chromosome)

sigt5$sigt7 <- ifelse(is.na(sigt5$sigt4), "0", sigt5$sigt4)
sigt5$sigt7 <- as.numeric(sigt5$sigt7)
p<-ggplot(data=sigt5, aes(x=reorder(chromosome, sigt7), y=sigt7)) +
  geom_bar(stat="identity") + coord_flip() + xlab("chromosome") + ylab("enrichment")
p




sigX <- subset(dfsig, chromRefined == "X")
dfX <- subset(df3, chromRefined == "X")




perm.pval <- data.frame(matrix(NA,nrow=length(sigt5$chromosome),ncol=3))
colnames(perm.pval)<- c("chromosome", "NumGenes", "pval")

sigt5vec <- as.list(as.character(sigt5$chromosome))

for (j in 1:length(sigt5$chromosome)){
  chrom_tmp <- sigt5vec[j]
  sigAll_tmp <- subset(dfsig, chromRefined == chrom_tmp)
  genenum <- length(sigAll_tmp$gene)
  
  if (genenum > 0) {
    null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
    colnames(null.dist.medians)<- "median_rank"
    for (i in 1:1000){
      set.seed(i)
      genelisttmp <- sample(dfsig$gene, genenum, replace = FALSE, prob = NULL)
      dfgenestmp <- filter(dfsig, gene %in% genelisttmp)
      median_rank <- median(dfgenestmp$corrSigRanked)
      null.dist.medians[i, 1] <- median_rank
    }
  } 
  newdata <- subset(sigt5, chromosome == chrom_tmp)
  if (newdata$sigt7 > 0){
    pval<- (sum(null.dist.medians$median_rank > median(sigAll_tmp$corrSigRanked)))/1000
  }
  if (newdata$sigt7 < 0){
    pval<- (sum(null.dist.medians$median_rank < median(sigAll_tmp$corrSigRanked)))/1000
  }
  perm.pval[j, 1]  <- chrom_tmp
  perm.pval[j, 2] <- genenum
  perm.pval[j, 3] <- pval
}


df_enrich_pval <- merge(sigt5, perm.pval, by = "chromosome")
df_enrich_pval1 <- df_enrich_pval[with(df_enrich_pval, order(-sigt7)),]




df3$corrSigRanked <- rank(df3$corr) #rank 15036 genes based on correlation with SVM weights

#plot for classical method, using df3
t3 <- sapply(levels(df3$chromRefined)[-1], function(x) median(df3$corrSigRanked[which(df3$chromRefined==x)],na.rm = TRUE)) # which(chromRefined==x) --  get an index of all the trues,,,   indexing out of t2Ranked where chromRefined==x.... find the median of t2Ranked for each of these levels
t4<- t3 - median(df3$corrSigRanked) # for the average correlation ranking for each chromosome, how far away it is from the median
t4
t5 <- as.data.frame(t4)
setDT(t5, keep.rownames = "chromosome")
t5$chromosome <- as.factor(t5$chromosome)

t5$t7 <- ifelse(is.na(t5$t4), "0", t5$t4)
t5$t7 <- as.numeric(t5$t7)
p<-ggplot(data=t5, aes(x=reorder(chromosome, t7), y=t7)) +
  geom_bar(stat="identity") + coord_flip() + xlab("Chromosome") + ylab("Enrichment")
p



#plot for classical method, using df3
t3 <- sapply(levels(df3$chromRefined)[-1], function(x) median(df3$corrSigRanked[which(df3$chromRefined==x)],na.rm = TRUE)) # which(chromRefined==x) --  get an index of all the trues,,,   indexing out of t2Ranked where chromRefined==x.... find the median of t2Ranked for each of these levels
t4<- t3 - median(df3$corrSigRanked) # for the average correlation ranking for each chromosome, how far away it is from the median
t4
t5 <- as.data.frame(t4)
setDT(t5, keep.rownames = "chromosome")
t5$chromosome <- as.factor(t5$chromosome)


t3se <- sapply(levels(df3$chromRefined)[-1], function(x) std.error(df3$corrSigRanked[which(df3$chromRefined==x)],na.rm = TRUE)) # which(chromRefined==x) --  get an index of all the trues,,,   indexing out of t2Ranked where chromRefined==x.... find the median of t2Ranked for each of these levels
t3se
t5se <- as.data.frame(t3se)
setDT(t5se, keep.rownames = "chromosome")
t5se$chromosome <- as.factor(t5$chromosome)
t8 <- merge (t5, t5se, by = "chromosome")

colnames(t8)<- c("chromosome", "Rank", "se")



fillcolor <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", 
                       "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                       "#FFFFFF", "#FFFFFF", "#FFFFFF", "grey39", "#FFFFFF",
                       "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                       "#FFFFFF", "#FFFFFF", "#FF0000", "#FFFFFF");


BorderColor <- c("grey39", "grey39", "grey39", "grey39", "grey39", 
                 "grey39", "grey39", "grey39", "grey39", "grey39",
                 "grey39", "grey39", "grey39", "grey39", "grey39",
                 "grey39", "grey39", "grey39", "grey39", "grey39",
                 "grey39", "grey39", "#FF0000", "#0000FF");

LineType <- c("dashed",  "dashed", "dashed", "dashed", "dashed",
              "dashed", "dashed", "dashed", "dashed", "dashed",
              "dashed", "dashed", "dashed", "solid", "dashed",
              "dashed", "dashed", "dashed", "dashed", "dashed",
              "dashed", "dashed", "solid", "dashed");

t9 <- cbind(t8, fillcolor, BorderColor, LineType)




#make color label an ordered factor so ggplot can match the color to network
t9$fillcolor <- factor(t9$fillcolor, ordered =T)
t9$BorderColor <- factor(t9$BorderColor, ordered =T)
t9$LineType <- factor(t9$LineType, ordered =T)


Fig <- ggplot(data_Gam_WholeNetworkSum, aes(EffectRank, SexEffects_Z)) +
  geom_bar(stat = "identity", position="identity", fill=c("#FF0000", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", 
                                                          "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                                                          "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                                                          "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                                                          "#FFFFFF", "#FFFFFF", "grey39", "#FFFFFF"),
           colour = BorderColor, linetype = LineType, width = 0.8) +
  labs(x = "Networks", y = expression(paste("Sex Association (", italic("Z"), ")"))) + theme_classic() + 
  theme(axis.text.x = element_text(size= 27, color = ColorScheme_XAxis),
        axis.text.y = element_text(size= 33, color = "black"),
        axis.title=element_text(size = 33)) +
  
  
  
  
  t9 <- t9 %>% 
  mutate(chromosome = as.factor(chromosome),fillcolor=as.character(fillcolor))%>%
  mutate(chromosome = as.factor(chromosome),BorderColor=as.character(BorderColor))%>%
  mutate(chromosome = as.factor(chromosome),LineType=as.character(LineType))%>%
  mutate(chromosome = reorder(chromosome,Rank))%>%
  arrange (chromosome)
colormap=t9 %>% select(chromosome,fillcolor)%>%unique()
linemap=t9 %>% select(chromosome,LineType)%>%unique()
bordermap=t9 %>% select(chromosome,BorderColor)%>%unique()


overall_se <- std.error(df3$corrSigRanked)


p<-ggplot(data=t9, aes(x=chromosome, y=Rank)) +
  geom_bar(stat="identity", fill=colormap$fillcolor ,
           colour = bordermap$BorderColor, linetype = linemap$LineType, width = 0.8) + geom_errorbar(aes(ymin=Rank-se, ymax=Rank+se), width=0, linetype = linemap$LineType, position=position_dodge(.9)) +
  geom_point(data=t9, aes(y=Rank, x=chromosome), size = 0.8) +
  coord_flip() + xlab("Chromosome") + ylab("Enrichment") +
  geom_hline(aes(yintercept = as.numeric(0))) +
  ylim(-2799, 1010) +
  #scale_y_discrete(breaks = c(-2000, -1000, 1000)) +
  theme(legend.position="none") + theme(axis.text.x = element_text(size= 10), axis.text.y = element_text(size= 10, color = bordermap$BorderColor), axis.title=element_text(size = 12))
p



p<-ggplot(data=t9, aes(x=chromosome, y=Rank)) +
  geom_bar(stat="identity", fill=colormap$fillcolor ,
           colour = bordermap$BorderColor, linetype = linemap$LineType, width = 0.8) + geom_errorbar(aes(ymin=Rank-se, ymax=Rank+se), width=0, size= 0.3, position=position_dodge(.9)) +
  geom_point(data=t9, aes(y=Rank, x=chromosome), size = 0.8) +
  coord_flip() + xlab("Chromosome") + ylab("Enrichment") +
  geom_hline(aes(yintercept = as.numeric(0))) +
  ylim(-2799, 1010) +
  #scale_y_discrete(breaks = c(-2000, -1000, 1000)) +
  theme(legend.position="none") + theme(axis.text.x = element_text(size= 10), axis.text.y = element_text(size= 10, color = bordermap$BorderColor), axis.title=element_text(size = 12))
p




perm.pval <- data.frame(matrix(NA,nrow=length(t5$chromosome),ncol=3))
colnames(perm.pval)<- c("chromosome", "NumGenes", "pval")

t5vec <- as.list(as.character(t5$chromosome))

for (j in 1:length(t5$chromosome)){
  chromosome_tmp <- t5vec[j]
  All_tmp <- subset(df3, chromRefined == chromosome_tmp)
  genenum <- length(All_tmp$gene)
  
  if (genenum > 0) {
    null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
    colnames(null.dist.medians)<- "median_rank"
    for (i in 1:1000){
      set.seed(i)
      genelisttmp <- sample(df3$gene, genenum, replace = FALSE, prob = NULL)
      dfgenestmp <- filter(df3, gene %in% genelisttmp)
      median_rank <- median(dfgenestmp$corrSigRanked)
      null.dist.medians[i, 1] <- median_rank
    }
  } 
  newdata <- subset(t5, chromosome == chromosome_tmp)
  if (newdata$t7 > 0){
    pval<- (sum(null.dist.medians$median_rank > median(All_tmp$corrSigRanked)))/1000
  }
  if (newdata$t7 < 0){
    pval<- (sum(null.dist.medians$median_rank < median(All_tmp$corrSigRanked)))/1000
  }
  perm.pval[j, 1]  <- chromosome_tmp
  perm.pval[j, 2] <- genenum
  perm.pval[j, 3] <- pval
}


df_enrich_pval <- merge(t5, perm.pval, by = "chromosome")
df_enrich_pval1 <- df_enrich_pval[with(df_enrich_pval, order(-t7)),]








########## ttest of X vs Y genes for each region, then correlate tstat with map

chromGeneRefined1 <- cbind(chromRefined_df, geneRefined_df)
Xgenes <- subset(df3, chromRefined == "X", select= "gene")
Ygenes <- subset(df3, chromRefined == "Y", select= "gene")

geneMatChrom_X_500 <- geneMatChrom_500_comp_fin[colnames(geneMatChrom_500_comp_fin) %in% Xgenes$gene]
geneMatChrom_X_500_t <- as.data.frame(t(geneMatChrom_X_500))

geneMatChrom_Y_500 <- geneMatChrom_500_comp_fin[colnames(geneMatChrom_500_comp_fin) %in% Ygenes$gene]
geneMatChrom_Y_500_t <- as.data.frame(t(geneMatChrom_Y_500))


XY_ttest <- data.frame(matrix(NA,nrow=403,ncol=2))
colnames(XY_ttest)<- c("parcNum", "tstat")

for (i in 1:403){
  XY_ttest[i, 1] <- colnames(geneMatChrom_X_500_t[i])
  tmp <- t.test(geneMatChrom_X_500_t[i], geneMatChrom_Y_500_t[i])
  XY_ttest[i,2] <- tmp$statistic
}

cor.test(img_geneMatChrom_roi_500_comp$wt_mean, XY_ttest$tstat)

XMeans <- colMeans(geneMatChrom_X_500_t)
cor.test(img_geneMatChrom_roi_500_comp$wt_mean, XMeans)

YMeans <- colMeans(geneMatChrom_Y_500_t)
cor.test(img_geneMatChrom_roi_500_comp$wt_mean, YMeans)


parcNum_vert_key <-data.frame(matrix(NA,nrow=10242,ncol=2))
colnames(parcNum_vert_key)<- c("parcNum", "vert")

for (i in 1:10242){
  for (j in 1:500){
    if (lh_schaefer1000[i] == data_parc_lh_parcnum$lh_schaefer1000[j]){
      parcNum_vert_key$parcNum[i] <- data_parc_lh_parcnum$parcNum[j]
      parcNum_vert_key$vert[i] <- data_parc_lh_parcnum$lh_schaefer1000[j]}}}



parcNum_vert_500 <- subset(data_parc_lh_parcnum_sort, select = -wt_mean)
XY_ttest_parc_lhschaefer1000 <- merge(XY_ttest, parcNum_vert_500, by = "parcNum")
XY_ttest_vert <-data.frame(matrix(NA,nrow=10242, ncol=1))
colnames(XY_ttest_vert) <- "XY_ttest_schaefer1000"

for (i in 1:10242){
  for (j in 1:403){
    if (lh_schaefer1000[i] == XY_ttest_parc_lhschaefer1000$lh_schaefer1000[j]){
      XY_ttest_vert$XY_ttest_schaefer1000[i] <- XY_ttest$tstat[j] }}}

writeMat("/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/XY_ttest_schaefer1000_net7.mat", XY_ttest_vert = XY_ttest_vert)






#scatterplot of tstat vs svm weights
Data_tmp1 = as.data.frame(img_geneMatChrom_roi_500_comp$wt_mean)
colnames(Data_tmp1) <- "wt_mean"
Data_tmp1$XYttest = as.numeric(XY_ttest$tstat);
cor.test(Data_tmp1$wt_mean, Data_tmp1$XYttest, method = "pearson")


ggplot(Data_tmp1, aes(x=wt_mean, y=XYttest, ymax=max(XYttest)*1.05)) +
  geom_point(position=position_jitter(width = NULL), stat="identity") +
  ylab("XY tstatistic") + xlab("Weights") + 
  stat_smooth(method = "lm") 









###X genes
Xgenes <- subset(df3, chromRefined == "X", select= "gene")
geneMatChrom_X_500 <- geneMatChrom_500_comp_fin[colnames(geneMatChrom_500_comp_fin) %in% Xgenes$gene]
geneMatChrom_X_500_t <- as.data.frame(t(geneMatChrom_X_500))


all_X_genes_median <- as.data.frame(rowMedians(as.matrix(geneMatChrom_X_500)))

colnames(all_X_genes_median)<- c("Xmed")
cor.test(img_geneMatChrom_roi_500_comp$wt_mean, all_X_genes_median$Xmed)


all_X_genes_mean <- as.data.frame(rowMeans(geneMatChrom_X_500))
all_X_genes_mean$parcNum <- rownames(geneMatChrom_X_500)
colnames(all_X_genes_mean)<- c("Xmean", "parcNum")
Xmean_parc_lhschaefer1000 <- merge(all_X_genes_mean, parcNum_vert_500, by = "parcNum")


all_X_genes_mean_vert <-data.frame(matrix(NA,nrow=10242, ncol=1))
colnames(all_X_genes_mean_vert) <- "mean_X_schaefer1000"

for (i in 1:10242){
  for (j in 1:403){
    if (lh_schaefer1000[i] == Xmean_parc_lhschaefer1000$lh_schaefer1000[j]){
      all_X_genes_mean_vert$mean_X_schaefer1000[i] <- all_X_genes_mean$Xmean[j] }}}


writeMat("/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/Xgenes_mean_schaefer1000_net7.mat", all_X_genes_mean_vert = all_X_genes_mean_vert)



###Y genes
Ygenes <- subset(df3, chromRefined == "Y", select= "gene")
geneMatChrom_Y_500 <- geneMatChrom_500_comp_fin[colnames(geneMatChrom_500_comp_fin) %in% Ygenes$gene]
geneMatChrom_Y_500_t <- as.data.frame(t(geneMatChrom_Y_500))

all_Y_genes_mean <- as.data.frame(rowMeans(geneMatChrom_Y_500))
all_Y_genes_mean$parcNum <- rownames(geneMatChrom_Y_500)
colnames(all_Y_genes_mean)<- c("Ymean", "parcNum")
Ymean_parc_lhschaefer1000 <- merge(all_Y_genes_mean, parcNum_vert_500, by = "parcNum")


all_Y_genes_mean_vert <-data.frame(matrix(NA,nrow=10242, ncol=1))
colnames(all_Y_genes_mean_vert) <- "mean_Y_schaefer1000"

for (i in 1:10242){
  for (j in 1:403){
    if (lh_schaefer1000[i] == Ymean_parc_lhschaefer1000$lh_schaefer1000[j]){
      all_Y_genes_mean_vert$mean_Y_schaefer1000[i] <- all_Y_genes_mean$Ymean[j] }}}


writeMat("/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/Ygenes_mean_schaefer1000_net7.mat", all_Y_genes_mean_vert = all_Y_genes_mean_vert)
