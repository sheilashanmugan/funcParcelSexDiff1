library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)

# read in LH parc
lh_schaefer400 <- read.table("/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/inputData/genetics/lh.schaefer400.txt")
lh_schaefer400 <- lh_schaefer400[1:10242,]
# read in color lookup table
lh_schaefer400_ct <- read.table("/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/inputData/genetics/lh.schaefer400.ct.txt")
lh_schaefer400_ct_200 <- as.data.frame(lh_schaefer400_ct[-1,])
colnames(lh_schaefer400_ct_200) <- "V1"
# read in brain img data
data_brain1 <- readMat("/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov/w_Brain_Sex_Abs_sum.mat")
data_brain <-data_brain1$w.Brain.Sex.Abs.sum.lh
# parcellate data
data_vec <- as.vector(data_brain)
data_parc_lh <- cbind(data_vec, lh_schaefer400) %>% as_tibble() %>% group_by(lh_schaefer400) %>% summarise(wt_mean = mean(data_vec))


# 
# # read in RH parc
# rh_schaefer400 <- read.table("/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/inputData/genetics/rh.schaefer400.txt")
# rh_schaefer400 <- rh_schaefer400[1:10242,]
# # read in color lookup table
# rh_schaefer400_ct <- read.table("/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/inputData/genetics/rh.schaefer400.ct.txt")
# # subset brain data
# data_brain <-data_brain1$w.Brain.Sex.Abs.sum.rh
# # parcellate data
# data_vec <- as.vector(data_brain)
# data_parc_rh <- cbind(data_vec, rh_schaefer400) %>% as_tibble() %>% group_by(rh_schaefer400) %>% summarise(wt_mean = mean(data_vec))

data_parc_lh_sort <- data_parc_lh[order(data_parc_lh$lh_schaefer400),]
data_parc_lh_sort_200 <- subset(data_parc_lh_sort, lh_schaefer400 != "65793", select=-lh_schaefer400)
# data_parc_rh_sort <- data_parc_rh[order(data_parc_rh$rh_schaefer400),]
# data_parc_rh_sort_200 <- subset(data_parc_rh_sort, rh_schaefer400 != "65793", select=-rh_schaefer400)

#data_parc_rm <- rowMeans(cbind(data_parc_lh,data_parc_rh)) #cant do this? left and right not symmetric?
#data_parc_lMr <- data_parc_lh - data_parc_rh



###read in gene data
gene <- read.csv("/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/inputData/genetics/genes_20647_schaefer_V2.csv",header=FALSE)
geneMat <- readMat("/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/inputData/genetics/gene_regional_expression_zscored_Schaefer_V2.mat")
geneMat1<- as.data.frame(geneMat$gene.regional.expression)
colnames(geneMat1) <- gene$V1
geneMat2 <- cbind(lh_schaefer400_ct_200, geneMat1)
geneMat2_sort <- geneMat2[order(geneMat2$V1),]
#geneMat2_sort_200 <- subset(geneMat2_sort, V1 != "65793", select=-V1)
geneMat2_sort_200 <- subset(geneMat2_sort, select=-V1)

#correlate sorted img and gene data (medial wall value removed in both)
img_gene_corr_mat <- sapply(1:20647, function(x) cor(data_parc_lh_sort_200,geneMat2_sort_200[,x]))
write.csv(gene[order(img_gene_corr_mat,decreasing = TRUE),], "~/Desktop/gene_ranked.csv", row.names = FALSE,  quote=FALSE) #use ordering of img_gene_corr_mat based on img_gene_corr_mat'sg values, ordeing based on correlations

#chech chromosome enrichement 
chrom <- read.csv("~/Desktop/projects/pfn_sex_diff/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
#geneMatChrom_200 <- geneMat2_sort_200[,gene$V1 %in% chrom$gene_symbol] #take geneMat2_sort_200 where gene$V1 and chrom$gene_symbol match

#genelist <- gene$V1[gene$V1 %in% chrom$gene_symbol] # take gene$V1 where gene$V1 and chrom$gene_symbol match
#geneMatChrom_200<- geneMat2_sort_200[,genelist] #


chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
#geneMatChrom_200<- geneMat2_sort_200[,geneRefined] #
geneMatChrom_200 <-geneMat2_sort_200[,match(geneRefined,colnames(geneMat2_sort_200),nomatch = 0)] #target gene first. where in #2 is #1
#geneMatChrom_200[,match(geneRefined,colnames(geneMatChrom_200))]


#when right, > colnames(geneMatChrom_200)[which(chromRefined=="Y")], will give same answer as geneMatChrom_200<- geneMatChrom_200[,geneRefined]



t2 <- sapply(1:15036, function(x) cor(data_parc_lh_sort_200,geneMatChrom_200[,x]))
t2Ranked <- rank(t2)
t3 <- sapply(levels(chromRefined)[-1], function(x) median(t2Ranked[which(chromRefined==x)],na.rm = TRUE)) # which(chromRefined==x) --  get an index of all the trues,,,   indexing out of t2Ranked where chromRefined==x.... find the median of t2Ranked for each of these levels
t4<- t3 - median(t2Ranked) # for the average correlation ranking for each chromosome, how far away it is from the median
t4
t5 <- as.data.frame(t4)
setDT(t5, keep.rownames = "chromosome")
t5$chromosome <- as.factor(t5$chromosome)
#t5$chromosome <- relevel(t5$chromosome, levels="1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "X", "Y")



p<-ggplot(data=t5, aes(x=reorder(chromosome, t4), y=t4)) +
  geom_bar(stat="identity") + coord_flip() + xlab("chromosome") + ylab("enrichment")
p





#####