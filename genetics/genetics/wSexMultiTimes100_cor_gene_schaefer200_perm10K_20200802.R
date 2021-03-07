
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
 #write.csv(gene[order(img_gene_corr_mat,decreasing = TRUE),], "/cbica/projects/funcParcelSexDiff/results/genetics/GO/gene_ranked_wSexMultiTimes100CorGene_20200804.csv", row.names = FALSE,  quote=FALSE) #use ordering of img_gene_corr_mat based on img_gene_corr_mat'sg values, ordeing based on correlations

#chech chromosome enrichement 
chrom <- read.csv("~/cbica/projects/funcParcelSexDiff/inputData/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_200 <-geneMat2_sort_200[,match(geneRefined,colnames(geneMat2_sort_200),nomatch = 0)] #target gene first. where in #2 is #1


geneRefined_df <- as.data.frame(geneRefined)
chromRefined_df <- as.data.frame(chromRefined)
chromGeneRefined <- cbind(chromRefined_df, geneRefined_df)


pval1 <- read.csv( "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img100_schaefer10K_gene_permSphereP_df_20200727.csv")
pval <- pval1[,-1]
df2 <- data.frame(t(pval))
colnames(df2) <- "pval"
df2$fdrp <- p.adjust(df2$pval, method="fdr", n=15036)
setDT(df2, keep.rownames = "gene")
chromGeneRefined$geneRefined <- gsub("-", "." , chromGeneRefined$geneRefined)

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



sigX <- subset(dfsig, chromRefined == "X" & corr > 0)
#write.csv(sigX, "/cbica/projects/funcParcelSexDiff/results/genetics/sigX.csv", row.names = FALSE,  quote=FALSE)



#using sample of 11 genes (all genes on X that survived FDR), compare to 483 that survive FDR
null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
colnames(null.dist.medians)<- "median_rank"
for (i in 1:1000){
  set.seed(i)
  genelist <- sample(dfsig$gene, 11, replace = FALSE, prob = NULL)
  df11genes <- filter(dfsig, gene %in% genelist)
  median_rank <- median(df11genes$corrSigRanked)
  null.dist.medians[i, 1] <- median_rank
}

dfsigX <- subset(dfsig, chromRefined == "X")
(sum(null.dist.medians$median_rank > median(dfsigX$corrSigRanked)))/1000



#using sample of 9 genes (genes with positive corr that survive FDR), compare to 274 w pos corr that survive FDR
dfsig_pos <- subset(dfsig, corr > 0)
null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
colnames(null.dist.medians)<- "median_rank"
for (i in 1:1000){
set.seed(i)
genelist <- sample(dfsig_pos$gene, 9, replace = FALSE, prob = NULL)
df9genes <- filter(dfsig_pos, gene %in% genelist)
median_rank <- median(df9genes$corrSigRanked)
null.dist.medians[i, 1] <- median_rank
}
median_rank_X <- median(sigX$corrSigRanked)
(sum(null.dist.medians$median_rank > median(sigX$corrSigRanked)))/1000



#using sample of 11 genes (all genes on X that survived FDR), compare to 15036 set
df3$corrSigRanked <- rank(df3$corr)
#using sample of 11 genes 
null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
colnames(null.dist.medians)<- "median_rank"
for (i in 1:1000){
  set.seed(i)
  genelist <- sample(df3$gene, 11, replace = FALSE, prob = NULL)
  dfgenes <- filter(df3, gene %in% genelist)
  median_rank <- median(dfgenes$corrSigRanked)
  null.dist.medians[i, 1] <- median_rank
}

dfsigX <- subset(dfsig, chromRefined == "X") #11 genes
geneSigX11 <- dfsigX$gene
dfgenesSigX11 <- filter(df3, gene %in% geneSigX11) #df that includes all the columns in "gene" where values correspond to the genes in geneSigX11
(sum(null.dist.medians$median_rank > median(dfgenesSigX11$corrSigRanked)))/1000


#using sample of 597 genes on X, compare to 15036 set
df3$corrSigRanked <- rank(df3$corr) #rank 15036 genes based on correlation with SVM weights
null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
colnames(null.dist.medians)<- "median_rank"
for (i in 1:1000){
  set.seed(i)
  genelist <- sample(df3$gene, 597, replace = FALSE, prob = NULL) #create sample of 597 random genes from 15036 genes -> genelist
  dfgenes <- filter(df3, gene %in% genelist) #get a df that includes all the columns in "gene" where values correspond to the genes in genelist
  median_rank <- median(dfgenes$corrSigRanked) #get the median rank of this sample of genes
  null.dist.medians[i, 1] <- median_rank # put median rank in df for comparison to emperical rank below
}
dfgenesAllX <- subset(df3, chromRefined == "X") #get subset of 597 genes that are on X
(sum(null.dist.medians$median_rank > median(dfgenesAllX$corrSigRanked)))/1000 #number of times median rank in null distribution is greater than epirical rank/1000



#using sample of 25 genes on y, compare to 15036 set
dfgenesAllY <- subset(df3, chromRefined == "Y") #get subset of genes that are on Y
numgenes <- length(dfgenesAllY$gene)
null.dist.medians <-data.frame(matrix(NA,nrow=1000,ncol=1))
colnames(null.dist.medians)<- "median_rank"
for (i in 1:1000){
  set.seed(i)
  genelist <- sample(df3$gene, numgenes, replace = FALSE, prob = NULL) #create sample of n random genes from 15036 genes -> genelist
  dfgenes <- filter(df3, gene %in% genelist) #get a df that includes all the columns in "gene" where values correspond to the genes in genelist
  median_rank <- median(dfgenes$corrSigRanked) #get the median rank of this sample of genes
  null.dist.medians[i, 1] <- median_rank # put median rank in df for comparison to emperical rank below
}
(sum(null.dist.medians$median_rank < median(dfgenesAllY$corrSigRanked)))/1000 #number of times median rank in null distribution is greater than epirical rank/1000



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


#odds ratio is (number of genes on X in 483 set/non-X genes in 483 set) / (genes on X in 15036 set that are not in 483 set/non-X genes in 15036 set that are not in 483 set)

NumGenesXin483 <- length(geneSigX11) #number of genes on X in 483 set
NumNonXin483 <- (length(dfsig$gene)-NumGenesXin483) #number of non-X genes in 483 set
dfgenesXfullNotThresh <- filter(dfgenesAllX, (gene %in% geneSigX11 == FALSE)) #genes on X in 15036 set that are not in 483 set
NumGenesXfullNotThresh <- length(dfgenesXfullNotThresh$gene)  #number of genes on X in 15036 set that are not in 483 set

dfNotX <- subset(df3, chromRefined != "X") #df of subset of nonX genes from 15036
geneDfgenesNotX <- dfNotX$gene #nonX genes from 15036 set
dfsigNotX <- subset(dfsig, chromRefined != "X") #df of non X genes in 483 set
geneDfsigNotX <- dfsigNotX$gene #nonX genes in 483 set
dfNonXfullNotThresh <- filter(dfNotX, (geneDfgenesNotX %in% geneDfsigNotX == FALSE)) #df of non-X genes in 15036 set that are not in 483 set
NumGenesNotXfullNotThresh <- length(dfNonXfullNotThresh$gene) #number of non-X genes in 15036 set that are not in 483 set
(NumGenesXin483 / NumNonXin483) / (NumGenesXfullNotThresh / NumGenesNotXfullNotThresh) #odds ratio

ORtable<-matrix(c(NumGenesXin483, 

