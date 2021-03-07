
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
cell2 <- read.csv("~/cbica/projects/funcParcelSexDiff/inputData/genetics/cellTypes/celltypes_PSP.csv")
cellRefined <- cell2$class[cell2$gene %in% gene$V1] #which values of cell2$Gene are in gene$V1 
geneRefined <- cell2$gene[cell2$gene %in% gene$V1] 
geneMatCell_200 <-geneMat2_sort_200[,match(geneRefined,colnames(geneMat2_sort_200),nomatch = 0)] #target gene first. where in #2 is #1


geneRefined_df <- as.data.frame(geneRefined)
cellRefined_df <- as.data.frame(cellRefined)
cellGeneRefined <- cbind(cellRefined_df, geneRefined_df)


pval1 <- read.csv( "/cbica/projects/funcParcelSexDiff/results/genetics/permutations/img100_schaefer10K_gene_permSphereP_df_20200727.csv")
pval <- pval1[,-1]
df2 <- data.frame(t(pval))
colnames(df2) <- "pval"
df2$fdrp <- p.adjust(df2$pval, method="fdr", n=15036)
setDT(df2, keep.rownames = "gene")
cellGeneRefined$geneRefined <- gsub("-", "." , cellGeneRefined$geneRefined)


t2 <- sapply(1:length(geneMatCell_200), function(x) cor(data_parc_lh_sort_200,geneMatCell_200[,x]))
cellGeneRefined$corr <- t2


#df3 <- merge(df2, cellGeneRefined, by.x="gene", by.y="geneRefined", all=TRUE)
df3 <- merge(df2, cellGeneRefined, by.x="gene", by.y="geneRefined")

#subset for genes surviving spin test
dfsig <- subset(df3, fdrp <= 0.05)
#Ranked these 184 observations based on correlation
dfsig$corrSigRanked <- rank(dfsig$corr)
#Calculate median ranking for each cell type
sigt3 <- sapply(levels(dfsig$cellRefined), function(x) median(dfsig$corrSigRanked[which(dfsig$cellRefined==x)],na.rm = TRUE)) # which(chromRefined==x) --  get an index of all the trues,,,   indexing out of t2Ranked where chromRefined==x.... find the median of t2Ranked for each of these levels
#subtract the overall median rank of the 184 observations from the median rank for each celltype
sigt4<- sigt3 - median(dfsig$corrSigRanked) # for the average correlation ranking for each chromosome, how far away it is from the median
sigt4
sigt5a <- as.data.frame(sigt4)
setDT(sigt5a, keep.rownames = "cellType")
sigt5a$cellType <- as.factor(sigt5a$cellType)

sigt5 <- sigt5a

sigt5$sigt7 <- ifelse(is.na(sigt5$sigt4), "0", sigt5$sigt4)
sigt5$sigt7 <- as.numeric(sigt5$sigt7)
p<-ggplot(data=sigt5, aes(x=reorder(cellType, sigt7), y=sigt7)) +
  geom_bar(stat="identity") + coord_flip() + xlab("cellType") + ylab("enrichment") +ggtitle("Subset surviving spin test")
p






perm.pval <- data.frame(matrix(NA,nrow=length(sigt5$cellType),ncol=3))
colnames(perm.pval)<- c("cellType", "NumGenes", "pval")

sigt5vec <- as.list(as.character(sigt5$cellType))

for (j in 1:length(sigt5$cellType)){
     cellType_tmp <- sigt5vec[j]
     sigAll_tmp <- subset(dfsig, cellRefined == cellType_tmp)
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
     newdata <- subset(sigt5, cellType == cellType_tmp)
      if (newdata$sigt7 > 0){
        pval<- (sum(null.dist.medians$median_rank > median(sigAll_tmp$corrSigRanked)))/1000
      }
     if (newdata$sigt7 < 0){
       pval<- (sum(null.dist.medians$median_rank < median(sigAll_tmp$corrSigRanked)))/1000
     }
    perm.pval[j, 1]  <- cellType_tmp
    perm.pval[j, 2] <- genenum
    perm.pval[j, 3] <- pval
}


df_enrich_pval <- merge(sigt5, perm.pval, by = "cellType")
df_enrich_pval1 <- df_enrich_pval[with(df_enrich_pval, order(-sigt7)),]






#plot for classical method, using df3
#Rank these 4558 observations based on correlation
df3$corrSigRanked <- rank(df3$corr)
#Calculate median ranking for each cell type
t3 <- sapply(levels(df3$cellRefined), function(x) median(df3$corrSigRanked[which(df3$cellRefined==x)],na.rm = TRUE)) # which(chromRefined==x) --  get an index of all the trues,,,   indexing out of t2Ranked where chromRefined==x.... find the median of t2Ranked for each of these levels
#subtract the overall median rank of the 4558 observations from the median rank for each celltype 
t4<- t3 - median(df3$corrSigRanked) # for the average correlation ranking for each cellType, how far away it is from the median
t4
t5a <- as.data.frame(t4)
setDT(t5a, keep.rownames = "cellType")
t5a$cellType <- as.factor(t5a$cellType)
t5 <- t5a

t5$t7 <- ifelse(is.na(t5$t4), "0", t5$t4)
t5$t7 <- as.numeric(t5$t7)
p<-ggplot(data=t5, aes(x=reorder(cellType, t7), y=t7)) +
  geom_bar(stat="identity") + coord_flip() + xlab("cellType") + ylab("enrichment") +ggtitle("classical")
p



perm.pval <- data.frame(matrix(NA,nrow=length(t5$cellType),ncol=3))
colnames(perm.pval)<- c("cellType", "NumGenes", "pval")

t5vec <- as.list(as.character(t5$cellType))

for (j in 1:length(t5$cellType)){
  cellType_tmp <- t5vec[j]
  All_tmp <- subset(df3, cellRefined == cellType_tmp)
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

