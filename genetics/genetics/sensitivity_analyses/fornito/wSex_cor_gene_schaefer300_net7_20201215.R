
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
lh_schaefer300 <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/schaefer300/lh_Schaefer300_7net_L.csv")
lh_schaefer300 <- lh_schaefer300[1:10242,]
# read in color lookup table
lh_schaefer300_ct <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/schaefer300/lh_Schaefer300_7net_cttable.csv")
lh_schaefer300_ct_150 <- as.data.frame(lh_schaefer300_ct[-1,])
colnames(lh_schaefer300_ct_150) <- "V1"
setDT(lh_schaefer300_ct_150, keep.rownames = "parcNum")

# read in brain img data
data_brain1 <- readMat("/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum.mat")
data_brain <-data_brain1$w.Brain.Sex.Abs.sum.lh
# parcellate data
data_vec <- as.vector(data_brain)
data_parc_lh <- cbind(data_vec, lh_schaefer300) %>% as_tibble() %>% group_by(lh_schaefer300) %>% summarise(wt_mean = mean(data_vec))

#order by parcel number for rotations
data_parc_lh <- subset(data_parc_lh, lh_schaefer300 != "65793")
data_parc_lh_parcnum <- data_parc_lh
data_parc_lh_parcnum$parcNum <- "NA"

for (i in 1:150){
  for (j in 1:150){
    if (data_parc_lh_parcnum$lh_schaefer300[i] == lh_schaefer300_ct_150$V1[j]){
      data_parc_lh_parcnum$parcNum[i] <- lh_schaefer300_ct_150$parcNum[j] }}}

data_parc_lh_parcnum$parcNum <- as.numeric(data_parc_lh_parcnum$parcNum)
data_parc_lh_parcnum_sort <- data_parc_lh_parcnum[order(data_parc_lh_parcnum$parcNum),]


###read in gene data
geneinfo <- readMat("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/AHBAprocessed/ROIxGene_Schaefer300_INT.mat")
geneinfo1 <- lapply(geneinfo, unlist, use.names=FALSE)
geneMat <- as.data.frame(geneinfo1$parcelExpression)
geneMat1 <- subset(geneMat, select = -V1)
gene <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/sensitivity_analyses/parcellation/data/genes/parcellations/schaefer300/GeneSymbol.csv",header=TRUE)
colnames(gene) <- "V1"
colnames(geneMat1) <- gene$V1
geneMat2 <- cbind(lh_schaefer300_ct_150, geneMat1)
geneMat2$parcNum <- as.numeric(geneMat2$parcNum)
geneMat2_sort <- geneMat2[order(geneMat2$parcNum),]
geneMat2_sort_150 <- subset(geneMat2_sort, select=-V1)


#merge img, gene, and centroid coord. remove parcels with missing gene data
img_geneMatChrom_150 <- merge(data_parc_lh_parcnum_sort,geneMat2_sort_150, by="parcNum")

img_geneMatChrom_150_comp <- subset(img_geneMatChrom_150, !is.nan(img_geneMatChrom_150[[5]])) ##if add more columns, change "10". This is the complete DF. Any subsets below should come from here/match this order.
data_parc_lh_sort_150_comp <- subset(img_geneMatChrom_150_comp, select=wt_mean)


othercols <- colnames(data_parc_lh_parcnum_sort)
geneMatChrom_150_comp <- img_geneMatChrom_150_comp[(length(othercols)+1):(length(colnames(img_geneMatChrom_150_comp)))]



#Select genes for chromosome enrichement 
chrom <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/Richiardi_Data_File_S2_ChrAdded.csv")
chromRefined <- chrom$chromosome_name[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneRefined <- chrom$gene_symbol[chrom$gene_symbol %in% gene$V1] #16000, which of those are in the 20, 000? indexing the 16000
geneMatChrom_150_comp_fin <-geneMatChrom_150_comp[,match(geneRefined,colnames(geneMatChrom_150_comp),nomatch = 0)] #target gene first. where in #2 is #1




geneRefined_df <- as.data.frame(geneRefined)
chromRefined_df <- as.data.frame(chromRefined)
chromGeneRefined <- cbind(chromRefined_df, geneRefined_df)






t2 <- sapply(1:12986, function(x) cor(data_parc_lh_sort_150_comp,geneMatChrom_150_comp_fin[,x]))
colt2<- as.vector(chromGeneRefined$geneRefined)
dft2 <- as.data.frame(cbind(colt2, t2))
colnames(dft2)<- c("gene", "corr")

df3 <- merge(dft2, chromGeneRefined, by.x="gene", by.y="geneRefined", all=TRUE)





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





