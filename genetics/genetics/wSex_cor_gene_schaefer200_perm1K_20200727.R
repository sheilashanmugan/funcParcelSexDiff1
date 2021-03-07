
library(R.matlab)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(data.table)

# read in LH parc
lh_schaefer400 <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/lh.schaefer400.txt")
lh_schaefer400 <- lh_schaefer400[1:10242,]
# read in color lookup table
lh_schaefer400_ct <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/lh.schaefer400.ct.txt")
lh_schaefer400_ct_200 <- as.data.frame(lh_schaefer400_ct[-1,])
colnames(lh_schaefer400_ct_200) <- "V1"
# read in brain img data
data_brain1 <- readMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum.mat")
data_brain <-data_brain1$w.Brain.Sex.Abs.sum.lh
# parcellate data
data_vec <- as.vector(data_brain)
data_parc_lh <- cbind(data_vec, lh_schaefer400) %>% as_tibble() %>% group_by(lh_schaefer400) %>% summarise(wt_mean = mean(data_vec))

data_parc_lh_sort <- data_parc_lh[order(data_parc_lh$lh_schaefer400),]
data_parc_lh_sort_200 <- subset(data_parc_lh_sort, lh_schaefer400 != "65793", select=-lh_schaefer400)







###read in gene data
gene <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/genes_20647_schaefer_V2.csv",header=FALSE)
geneMat <- readMat("/cbica/projects/funcParcelSexDiff/inputData/genetics/gene_regional_expression_zscored_Schaefer_V2.mat")
geneMat1<- as.data.frame(geneMat$gene.regional.expression)
colnames(geneMat1) <- gene$V1
geneMat2 <- cbind(lh_schaefer400_ct_200, geneMat1)
geneMat2_sort <- geneMat2[order(geneMat2$V1),]
geneMat2_sort_200 <- subset(geneMat2_sort, select=-V1)

##correlate sorted img and gene data (medial wall value removed in both) for other enrichement analyses
#img_gene_corr_mat <- sapply(1:20647, function(x) cor(data_parc_lh_sort_200,geneMat2_sort_200[,x]))
#write.csv(gene[order(img_gene_corr_mat,decreasing = TRUE),], "~/Desktop/projects/pfn_sex_diff/genetics/gene_ranked.csv", row.names = FALSE,  quote=FALSE) #use ordering of img_gene_corr_mat based on img_gene_corr_mat'sg values, ordeing based on correlations

#chech chromosome enrichement 
chrom <- read.csv("~/Desktop/projects/pfn_sex_diff/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_200 <-geneMat2_sort_200[,match(geneRefined,colnames(geneMat2_sort_200),nomatch = 0)] #target gene first. where in #2 is #1


geneRefined_df <- as.data.frame(geneRefined)
chromRefined_df <- as.data.frame(chromRefined)
chromGeneRefined <- cbind(chromRefined_df, geneRefined_df)


pval1 <- read.csv( "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img100_gene_permSphereP_df_20200721.csv")
pval <- pval1[,-1]
df2 <- data.frame(t(pval))
colnames(df2) <- "pval"
df2$fdrp <- p.adjust(df2$pval, method="fdr", n=15036)
setDT(df2, keep.rownames = "gene")
df2$gene <- as.character(df2$gene)

chromGeneRefined$geneRefined <- as.character(chromGeneRefined$geneRefined)
df2$gene <- gsub("\\s+", "" , df2$gene)
chromGeneRefined$geneRefined <- gsub("\\s+", "" , df2$gene)


df3 <- merge(df2, chromGeneRefined, by.x="gene", by.y="geneRefined", all=TRUE)



t2 <- sapply(1:15036, function(x) cor(data_parc_lh_sort_200,geneMatChrom_200[,x]))
df3$corr <- t2

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



