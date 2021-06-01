
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)
library(ggplot2)
library(hexbin)
library(plotrix)
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
#Cell type
cell3 <- read.csv("~/cbica/projects/funcParcelSexDiff/inputData/genetics/cellTypes/Lake18_celltypes.csv")
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


#Chromosome
chrom <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_500_comp_fin <-geneMatChrom_500_comp[,match(geneRefined,colnames(geneMatChrom_500_comp),nomatch = 0)] #target gene first. where in #2 is #1


geneRefined_df <- as.data.frame(geneRefined)
chromRefined_df <- as.data.frame(chromRefined)
chromGeneRefined <- cbind(chromRefined_df, geneRefined_df)


chromGeneRefined$geneRefined <- gsub("-", "." , chromGeneRefined$geneRefined)
t2g <- sapply(1:12986, function(x) cor(data_parc_lh_sort_500_comp,geneMatChrom_500_comp_fin[,x]))
colt2g<- as.vector(chromGeneRefined$geneRefined)
dft2g <- as.data.frame(cbind(colt2g, t2g))
colnames(dft2g)<- c("gene", "corr")
df3g <- merge(dft2g, chromGeneRefined, by.x="gene", by.y="geneRefined", all=TRUE)


#merge chromosome and cell tpye 
dfall <- merge(df3, df3g, by="gene")
df_SigCelltype <- subset(dfall, cellRefined == "Ex5b" | cellRefined == "Ex1" | cellRefined ==  "Ex3e" | cellRefined == "Ex6b" | cellRefined == "Ex2" | cellRefined == "Ast")
df_SigCelltype_X <- subset(df_SigCelltype, chromRefined == "X", select=c("cellRefined", "gene"))
