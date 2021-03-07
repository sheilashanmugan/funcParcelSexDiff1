
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)



geneMatChrom_500_comp_fin <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/geneMatChrom_500_comp_fin.csv", header=TRUE)
data_parc_lh_sort_500_comp <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/data_parc_lh_sort_500_comp.csv", header=TRUE)
img_geneMatChrom_roi_500_comp <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img_geneMatChrom_roi_500_comp.csv", header=TRUE)

geneMatChrom_500_comp_fin <- geneMatChrom_500_comp_fin[,-1]
img_geneMatChrom_roi_500_comp <- img_geneMatChrom_roi_500_comp[,-1]
system("echo line19")



library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(data.table)
library(parallel)
source("/cbica/projects/funcParcelSexDiff/scripts/rFunctions/perm.sphere.p.R")
library(matrixStats)
source("/cbica/projects/funcParcelSexDiff/scripts/rFunctions/rotate.parcellation.R")

#brain input is from SVM run 100x
# read in LH parc
lh_schaefer1000 <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/fsaverage5/lh_Schaefer1000_17net_fs5_L.txt")
lh_schaefer1000 <- lh_schaefer1000[1:10242,]
# read in color lookup table
lh_schaefer1000_ct <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/fsaverage5/lh_schaefer1000_17net_fs5_ct.txt")
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

system("echo line 43")

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




system("echo line 58")
#merge img, gene, and centroid coord. remove parcels with missing gene data
img_geneMatChrom_500 <- merge(data_parc_lh_parcnum_sort,geneMat2_sort_500, by="parcNum")


system("echo line 62")
rois <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mmCentroid_RAS.csv")
roislh <- rois[1:500,]
system("echo line 65")

img_geneMatChrom_roi_500 <- merge(roislh, img_geneMatChrom_500, by.x = "RoiLabel", by.y= "parcNum",all = FALSE)
img_geneMatChrom_roi_500_comp <- subset(img_geneMatChrom_roi_500, !is.nan(img_geneMatChrom_roi_500[[10]])) ##if add more columns, change "10". This is the complete DF. Any subsets below should come from here/match this order.

system("echo line 71")

data_parc_lh_sort_500_comp <- subset(img_geneMatChrom_roi_500_comp, select=wt_mean)


othercols <- c(colnames(rois), "lh_schaefer1000", "wt_mean")
geneMatChrom_500_comp <- img_geneMatChrom_roi_500_comp[(length(othercols)+1):(length(colnames(img_geneMatChrom_roi_500_comp)))]

#Select genes for chromosome enrichement 
chrom <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_500_comp_fin <-geneMatChrom_500_comp[,match(geneRefined,colnames(geneMatChrom_500_comp),nomatch = 0)] #target gene first. where in #2 is #1




geneRefined_df <- as.data.frame(geneRefined)
chromRefined_df <- as.data.frame(chromRefined)
chromGeneRefined <- cbind(chromRefined_df, geneRefined_df)


t2 <- sapply(1:12986, function(x) cor(data_parc_lh_sort_500_comp,geneMatChrom_500_comp_fin[,x]))
colt2<- as.vector(chromGeneRefined$geneRefined)
df2 <- as.data.frame(cbind(colt2, t2))
colnames(df2)<- c("gene", "corr")
df3 <- merge(df2, chromGeneRefined, by.x="gene", by.y="geneRefined")
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
  geom_bar(stat="identity") + coord_flip() + xlab("chromosome") + ylab("enrichment")
p


null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
colnames(null.dist.medians)<- "median_rank"
for (i in 1:1000){
  set.seed(i)
  genelist <- sample(df3$gene, 449, replace = FALSE, prob = NULL) #create sample of 597 random genes from 15036 genes -> genelist
  dfgenes <- filter(df3, gene %in% genelist) #get a df that includes all the columns in "gene" where values correspond to the genes in genelist
  median_rank <- median(dfgenes$corrSigRanked) #get the median rank of this sample of genes
  null.dist.medians[i, 1] <- median_rank # put median rank in df for comparison to emperical rank below
}
dfgenesAllX <- subset(df3, chromRefined == "X") #get subset of 597 genes that are on X
(sum(null.dist.medians$median_rank > median(dfgenesAllX$corrSigRanked)))/1000 #number of times median rank in null distribution is greater than epirical rank/1000

