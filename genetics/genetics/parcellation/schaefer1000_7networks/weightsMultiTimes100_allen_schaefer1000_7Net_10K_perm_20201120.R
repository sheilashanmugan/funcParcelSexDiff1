
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
rois <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/github/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mmCentroid_RAS.csv")
roislh <- rois[1:500,]
system("echo line 65")

img_geneMatChrom_roi_500 <- merge(roislh, img_geneMatChrom_500, by.x = "RoiLabel", by.y= "parcNum",all = FALSE)
img_geneMatChrom_roi_500_comp <- subset(img_geneMatChrom_roi_500, !is.nan(img_geneMatChrom_roi_500[[10]])) ##if add more columns, change "10". This is the complete DF. Any subsets below should come from here/match this order.
write.csv(img_geneMatChrom_roi_500_comp, "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img_geneMatChrom_roi_500_comp_7net.csv")

system("echo line 71")

data_parc_lh_sort_500_comp <- subset(img_geneMatChrom_roi_500_comp, select=wt_mean)
write.csv(data_parc_lh_sort_500_comp, "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/data_parc_lh_sort_500_comp_7net.csv")


othercols <- c(colnames(rois), "lh_schaefer1000", "wt_mean")
geneMatChrom_500_comp <- img_geneMatChrom_roi_500_comp[(length(othercols)+1):(length(colnames(img_geneMatChrom_roi_500_comp)))]

#Select genes for chromosome enrichement 
chrom <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_500_comp_fin <-geneMatChrom_500_comp[,match(geneRefined,colnames(geneMatChrom_500_comp),nomatch = 0)] #target gene first. where in #2 is #1

write.csv(geneMatChrom_500_comp_fin, "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/geneMatChrom_500_comp_fin_7net.csv")
# 
# coord.l <- data.matrix(img_geneMatChrom_roi_500_comp[1:length(img_geneMatChrom_roi_500_comp$R),3:5])
# coord.r <- data.matrix(rois[501:1000,3:5])
# 
# schaefer1000_perm <- rotate.parcellation(coord.l,coord.r, nrot=10000)
# schaefer1000_perm_lh <- schaefer1000_perm[1:length(img_geneMatChrom_roi_500_comp$R),]
# write.csv(schaefer1000_perm_lh, "/cbica/projects/funcParcelSexDiff/inputData/genetics/permutations/Schaefer1000_perm_403_500.csv")
# 
# 
# #read in permutations of schaefer left hemisphere 
# #schaefer_perm1 <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/permutations/Schaefer1000_perm.csv", header=TRUE)
# #schaefer_perm <- subset(schaefer_perm1, select=-X)
# #schaefer1000_perm_lh <- schaefer1000_perm[1:length(img_geneMatChrom_roi_500_comp$R),]
# 
# #make vector for perm,sphere.p
# # data_parc_lh_sort_500_comp_vec <- as.vector(data_parc_lh_sort_500_comp$wt_mean)
# # 
# # set.seed(1)
# # img_gene_permSphereP <- mclapply(geneMatChrom_500_comp_fin, FUN=perm.sphere.p, y=data_parc_lh_sort_500_comp_vec, perm.id= schaefer1000_perm_lh, corr.type='pearson', mc.set.seed =  TRUE, mc.cores=getOption("mc.cores", 20))
# # 
# # 
# # write.csv(img_gene_permSphereP, "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img100_schaefer1000_perm10K_SortParcNum_gene_permSphereP_df_20201023.csv")
# # 
