
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)
library(ggplot2)
library(hexbin)
library(matrixStats)

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







###read in gene data
#check cell enrichement 
cell3a <- read.csv("~/cbica/projects/funcParcelSexDiff/inputData/genetics/cellTypes/Lake18_celltypes.csv")
cell3 <- subset(cell3a, grepl("^Ex", cell3a$Cluster) | grepl("^In", cell3a$Cluster) )
 genebrain <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/cellTypes/brain_genes_HPA.csv", header = FALSE)
cell2 <- merge(cell3, genebrain, by.x="Gene", by.y="V1")

cellRefined <- cell2$Cluster[cell2$Gene %in% gene$V1] #which values of cell2$Gene are in gene$V1 
geneRefined <- cell2$Gene[cell2$Gene %in% gene$V1] 
geneMatCell_500 <-geneMatChrom_500_comp[,match(geneRefined,colnames(geneMatChrom_500_comp),nomatch = 0)] #target gene first. where in #2 is #1


geneRefined_df <- as.data.frame(geneRefined)
cellRefined_df <- as.data.frame(cellRefined)
cellGeneRefined <- cbind(cellRefined_df, geneRefined_df)





t2 <- sapply(1:length(geneMatCell_500), function(x) cor(data_parc_lh_sort_500_comp,geneMatCell_500[,x]))
colt2<- as.vector(cellGeneRefined$geneRefined)
dft2 <- as.data.frame(cbind(colt2, t2))
colnames(dft2)<- c("gene", "corr")

df3a <- merge(dft2, cellGeneRefined, by.x="gene", by.y="geneRefined")
df3 <- unique(df3a) 



#plot for classical method, using df3
#Rank these 4558 observations based on correlation
df3$corrSigRanked <- rank(df3$corr)
#Calculate median ranking for each cell type
t3 <- sapply(levels(df3$cellRefined), function(x) median(df3$corrSigRanked[which(df3$cellRefined==x)],na.rm = TRUE)) # which(chromRefined==x) --  get an index of all the trues,,,   indexing out of t2Ranked where chromRefined==x.... find the median of t2Ranked for each of these levels
#subtract the overall median rank of the 4965 observations from the median rank for each celltype 
t4<- t3 - median(df3$corrSigRanked) # for the average correlation ranking for each cellType, how far away it is from the median
t4
t5a <- as.data.frame(t4)
setDT(t5a, keep.rownames = "cellType")
t5a$cellType <- as.factor(t5a$cellType)
#remove the non neuron labels
t5 <- subset(t5a, grepl("^Ex", t5a$cellType) | grepl("^In", t5a$cellType) )

t5$t7 <- ifelse(is.na(t5$t4), "0", t5$t4)
t5$t7 <- as.numeric(t5$t7)
p<-ggplot(data=t5, aes(x=reorder(cellType, t7), y=t7)) +
  geom_bar(stat="identity") + coord_flip() + xlab("cellType") + ylab("enrichment") +ggtitle("brain expressed")
p



perm.pval <- data.frame(matrix(NA,nrow=length(t5$cellType),ncol=3))
colnames(perm.pval)<- c("cellType", "NumGenes", "pval")

t5vec <- as.list(as.character(t5$cellType))
df3gene_unique<- unique(df3$gene)

for (j in 1:length(t5$cellType)){
  cellType_tmp <- t5vec[j]
  All_tmp <- subset(df3, cellRefined == cellType_tmp)
  genenum <- length(All_tmp$gene)
  
  if (genenum > 0) {
    null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
    colnames(null.dist.medians)<- "median_rank"
    for (i in 1:1000){
      set.seed(i)
      genelisttmp <- sample(df3gene_unique, genenum, replace = FALSE, prob = NULL)
      dfgenestmp <- filter(df3, gene %in% genelisttmp)
      median_rank <- median(dfgenestmp$corrSigRanked)
      null.dist.medians[i, 1] <- median_rank
    }
  } 
  newdata <- subset(t5, cellType == cellType_tmp)
  if (newdata$t7 > 0){
    pval<- (sum(null.dist.medians$median_rank > median(All_tmp$corrSigRanked)))/1000
  }
  if (newdata$t7 < 0){
    pval<- (sum(null.dist.medians$median_rank < median(All_tmp$corrSigRanked)))/1000
  }
  perm.pval[j, 1]  <- cellType_tmp
  perm.pval[j, 2] <- genenum
  perm.pval[j, 3] <- pval
}


df_enrich_pval <- merge(t5, perm.pval, by = "cellType")
df_enrich_pval1 <- df_enrich_pval[with(df_enrich_pval, order(-t7)),]



