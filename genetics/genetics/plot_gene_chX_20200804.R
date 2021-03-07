


# read in LH parc
lh_schaefer400 <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/lh.schaefer400.txt")
lh_schaefer400 <- lh_schaefer400[1:10242,]
# read in color lookup table
lh_schaefer400_ct <- read.table("/cbica/projects/funcParcelSexDiff/inputData/genetics/lh.schaefer400.ct.txt")
lh_schaefer400_ct_200 <- as.data.frame(lh_schaefer400_ct[-1,])
colnames(lh_schaefer400_ct_200) <- "V1"

###
gene <- read.csv("/cbica/projects/funcParcelSexDiff/inputData/genetics/genes_20647_schaefer_V2.csv",header=FALSE)
geneMat <- readMat("/cbica/projects/funcParcelSexDiff/inputData/genetics/gene_regional_expression_zscored_Schaefer_V2.mat")
geneMat1<- as.data.frame(geneMat$gene.regional.expression)
colnames(geneMat1) <- gene$V1
geneMat2 <- cbind(lh_schaefer400_ct_200, geneMat1)
geneMat2_sort <- geneMat2[order(geneMat2$V1),]

sigX <- read.csv("/cbica/projects/funcParcelSexDiff/results/genetics/sigX.csv")
writeMat("/cbica/projects/funcParcelSexDiff/results/genetics/sigX.mat", sigX = sigX)

namelist <- sigX$gene
all_genes_vert <-data.frame(matrix(NA,nrow=length(namelist),ncol=10242))
rownames(all_genes_vert)<- namelist


for(g in 1:length(namelist)){
  gene_var <- as.character(namelist[g])
  gene_var_mat <- subset(geneMat2_sort, select= c("V1", gene_var))
  colnames(gene_var_mat) <- c("V1", "gene_var")
  gene_vert<-matrix(NA,nrow=10242,ncol=1)
  for (i in 1:200){
    for (j in 1:10242){
      if (gene_var_mat$V1[i] == lh_schaefer400[j]){
        gene_vert[j] <- gene_var_mat$gene_var[i] }}}
  matPath <- "/cbica/projects/funcParcelSexDiff/results/genetics/gene_mat/"
  writeMat(paste0(matPath, gene_var, '_vert.mat'), gene_vert = gene_vert)
  all_genes_vert[g, 1:10242] <- gene_vert
} 

all_genes_mean <- as.data.frame(colMeans(all_genes_vert[sapply(all_genes_vert, is.numeric)]))
colnames(all_genes_mean) <- "gene_avg"
matPath <- "/cbica/projects/funcParcelSexDiff/results/genetics/gene_mat/"
writeMat("/cbica/projects/funcParcelSexDiff/results/genetics/gene_mat/all_genes_mean.mat", all_genes_mean = all_genes_mean)
