
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
t2 <- sapply(1:15036, function(x) cor(data_parc_lh_sort_200,geneMatChrom_200[,x]))
df2$corr <- t2
chromGeneRefined$geneRefined <- gsub("-", "." , chromGeneRefined$geneRefined)



df3 <- merge(df2, chromGeneRefined, by.x="gene", by.y="geneRefined", all=TRUE)

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


sigX <- subset(dfsig, chromRefined == "X")




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
    perm.pval[j, 3] <- pval
  }
  if (newdata$sigt7 < 0){
    pval<- (sum(null.dist.medians$median_rank < median(sigAll_tmp$corrSigRanked)))/1000
    perm.pval[j, 3] <- pval
  }
  perm.pval[j, 1]  <- chrom_tmp
  perm.pval[j, 2] <- genenum
  
}

df_enrich_pval <- merge(sigt5, perm.pval, by = "chromosome")
df_enrich_pval1 <- df_enrich_pval[with(df_enrich_pval, order(-sigt7)),]




df3$corrSigRanked <- rank(df3$corr)
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




