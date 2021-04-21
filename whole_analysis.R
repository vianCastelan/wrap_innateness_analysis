library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("genefilter")
library("vsn")
library("pheatmap")
library("gplots")
library("genefilter")
output_dir="/mnt/BioAdHoc/Groups/vd-vijay/vcastelan/innateness_analysis/lmm_test_1000mvg_PC1/wrap"#set output dir
setwd(paste(output_dir,"/data",sep=""))

counts <- read.csv("GSE109125_Gene_count_table_filtered_justnaive_quit_ILCs.csv",
                   header=TRUE, row.names=1)
col_data <- read.csv("master_tpm_metadata_2020-12-01_justnaive_quit_ILCs.csv",
                     header=TRUE, row.names=1)

#make a DESeq Dataset as an object. This will take into account out counts and activation data, as well as experimental design (experimental variable being set as condition):
colnames(counts) <- NULL

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col_data,
                              design = ~ cell_type)
#DESeq analysis and summarize the results with res
#performs an internal normalization where geometric mean is calculated for each gene across all samples. 
#The counts for a gene in each sample is then divided by this mean.
#Differential Expression Analysis Based On The Negative Binomial (A.K.A. Gamma-Poisson) Distribution
dds <- DESeq(dds)
res <- results(dds)

#I calculate highly variable gene (HVG), since that discovery allows the detection of genes that contribute strongly to cell-to-cell variation within a homogeneous cell population

rld <- rlog(dds) # to get the log transform of gene expression
vsd <- varianceStabilizingTransformation(dds)
ntd <- normTransform(dds) #log2 transformation log2(n + 1)

pdf(file=paste(output_dir,"/whole_analysis.pdf", sep = ""))
#Heatmap of the sample-to-sample distances
#Get sample-to-sample distances.
#gives us an overview over similarities and dissimilarities between samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$all_conditions
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
ntop <- 1000
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(assay(vsd)[select, ])
pca<-prcomp(mat,scale= TRUE)
attributes(pca)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2))

#write.csv(pca$x[,], file="PCs_1000genes_quit_ILCs.csv")

data_ct<-as.data.frame(pca$x[,])
### PC1 vs PC2
#---cell type
ggplot(data_ct, aes(pca$x[,1], pca$x[,2], color=as.factor(col_data$all_conditions))) +
  geom_point(size=3) +
  ggtitle("PC1 vs PC2") +
  labs(color = "Cell type")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

### PC1 vs PC3
#---cell type
ggplot(data_ct, aes(pca$x[,1], pca$x[,3], color=as.factor(col_data$all_conditions))) +
  geom_point(size=3) +
  ggtitle("PC1 vs PC3") +
  labs(color = "Cell type")+
  xlab(paste0("PC1 ",percentVar[1],"% variance")) +
  ylab(paste0("PC3 ",percentVar[3],"% variance"))


### PC2 vs PC3
#---cell type
ggplot(data_ct, aes(pca$x[,2], pca$x[,3], color=as.factor(col_data$cell_type))) +
  geom_point(size=3) +
  ggtitle("PC2 vs PC3") +
  labs(color = "Cell type")+
  xlab(paste0("PC2 ",percentVar[2],"% variance")) +
  ylab(paste0("PC3 ",percentVar[3],"% variance"))



#########--- lmm test ---###########

library(tidyverse)
library(lme4)
library(lmerTest) #### install.packages("lmerTest") to get better pval from ANOVA-like test for Random Effects


pc_coord <- read.table(file = "PCs_1000genes_quit_ILCs.csv", sep = ",", row.names = 1, header = TRUE) ## first run with 1000 HVG and ILCs by Vianka
col_data <- read.csv(file = "master_tpm_metadata_2020-12-01_justnaive_quit_ILCs.csv") ## still has activated T cells
gene_exp <- read.csv(file = "GSE109125_Gene_count_table_filtered_justnaive_quit_ILCs.csv" , header = TRUE, check.names = FALSE) ## DESeq2 by Vianka
#col_data <- col_data
#gene_exp <- counts

col_data <- col_data[col_data$names %in% row.names(pc_coord),] ## filter col_data for samples in pc_coord
gene_exp %>% select(gene_symbol, col_data$names) -> gene_exp ## filter gen_exp for samples in col_data
####################################### single gene selection test ##########################################
gene = "Ccl5" # "Ifng" # "Tbx21" # "Sell" # GzmB #Runx3
# selecting gene expression
gene_exp[gene_exp$gene_symbol == gene,]
test_gene <- unname(t(gene_exp[gene_exp$gene_symbol == gene,][2:22])) ## remove names
attributes(test_gene) <- NULL ## remove attributes
# selecting PC coordenates
test_pc <- pc_coord$PC1
## master table
master_table <- col_data
master_table$PC1 <- test_pc
master_table$gene <- test_gene
#####################################  function for overdispertion ##########################################
## obtained from http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html ### UPDATE : this pval estimation based on overdispertion is not correct.
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
########################################## LMM with just one gene ##########################################
test_model <- lmer(gene ~ PC1 + (1|cell_type), data = master_table) ## correct formula for random effect of cell_type /// REML
summary(test_model)$coefficients[6] ## t-statistic for fixed slope (beta)
summary(test_model)$coefficients[1] ## mean intercept for random effects
summary(test_model)$coefficients[2] ## fixed slope (beta) value
#overdisp_fun(test_model)[3] ## degrees of freedom
#2*pt(summary(test_model)$coefficients[6], df = overdisp_fun(test_model)[3], lower=FALSE) ## pval for slope (beta)

### new PVAL estimation
## based on https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode
ranova(test_model)[["Pr(>Chisq)"]][[2]] ## pval from ANOVA-like table for random-effects (lmerTest package)  <---- use this PVAL
anova(test_model) ## this give you a F-statistic, but no degree of freedom to test. Can't use.
VarCorr(test_model) ## Variance Estimates
ranef(test_model) ## random effect intercepts
confint(test_model) ## confidence intervals
confint(test_model,parm="PC1",method="boot")  ## bootstrap estimation of confidence intervals ## WARNING: this also takes a long time to compute!!!!
### Add CONFINT to table
#as.data.frame(VarCorr(test_model))$vcov[1] / (as.data.frame(VarCorr(test_model))$vcov[1] + as.data.frame(VarCorr(test_model))$vcov[2]) ## Variance explained by cell_type
## plot test
print("Plot 1, gene test")
par(mar=c(5,4,4,2))
plot(gene ~PC1, col=as.factor(master_table$cell_type), pch=16, cex=1,data=master_table)
abline(a = summary(test_model)$coefficients[1], b = summary(test_model)$coefficients[2], lwd=1)
text(gene ~PC1, labels=cell_type,data=master_table,cex=0.5, font=1, pos=4)

#FROM HERE
####################################### single gene selection test ##########################################
gene2 = "Nkg7" # "Ifng" # "Tbx21" # "Sell" # GzmB #Runx3
# selecting gene expression
gene_exp[gene_exp$gene_symbol == gene2,]
test_gene2 <- unname(t(gene_exp[gene_exp$gene_symbol == gene2,][2:22])) ## remove names
attributes(test_gene2) <- NULL ## remove attributes
# selecting PC coordenates
test_pc <- pc_coord$PC1
## master table
master_table <- col_data
master_table$PC1 <- test_pc
master_table$gene2 <- test_gene2

########################################## LMM with just one gene2 ##########################################
test_model2 <- lmer(gene2 ~ PC1 + (1|cell_type), data = master_table) ## correct formula for random effect of cell_type /// REML
summary(test_model2)$coefficients[6] ## t-statistic for fixed slope (beta)
summary(test_model2)$coefficients[1] ## mean intercept for random effects
summary(test_model2)$coefficients[2] ## fixed slope (beta) value

### new PVAL estimation
## based on https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode
ranova(test_model2)[["Pr(>Chisq)"]][[2]] ## pval from ANOVA-like table for random-effects (lmerTest package)  <---- use this PVAL
anova(test_model2) ## this give you a F-statistic, but no degree of freedom to test. Can't use.
VarCorr(test_model2) ## Variance Estimates
ranef(test_model2) ## random effect intercepts
confint(test_model2) ## confidence intervals
confint(test_model2,parm="PC1",method="boot")  ## bootstrap estimation of confidence intervals ## WARNING: this also takes a long time to compute!!!!
### Add CONFINT to table
#as.data.frame(VarCorr(test_model2))$vcov[1] / (as.data.frame(VarCorr(test_model2))$vcov[1] + as.data.frame(VarCorr(test_model2))$vcov[2]) ## Variance explained by cell_type
## plot test
print("Plot 1, gene test")
par(mar=c(5,4,4,2))
plot(gene2 ~PC1, col=as.factor(master_table$cell_type), pch=16, cex=1,data=master_table)
abline(a = summary(test_model2)$coefficients[1], b = summary(test_model2)$coefficients[2], lwd=1)
text(gene2 ~PC1, labels=cell_type,data=master_table,cex=0.5, font=1, pos=4)

####################################################################################################################################################################################################################################################################################################
####################################### single gene selection test ##########################################
gene3 = "Lef1" # "Ifng" # "Tbx21" # "Sell" # GzmB #Runx3
# selecting gene expression
gene_exp[gene_exp$gene_symbol == gene3,]
test_gene3 <- unname(t(gene_exp[gene_exp$gene_symbol == gene3,][2:22])) ## remove names
attributes(test_gene3) <- NULL ## remove attributes
# selecting PC coordenates
test_pc <- pc_coord$PC1
## master table
master_table <- col_data
master_table$PC1 <- test_pc
master_table$gene3 <- test_gene3
########################################## LMM with just one gene3 ##########################################
test_model3 <- lmer(gene3 ~ PC1 + (1|cell_type), data = master_table) ## correct formula for random effect of cell_type /// REML
summary(test_model3)$coefficients[6] ## t-statistic for fixed slope (beta)
summary(test_model3)$coefficients[1] ## mean intercept for random effects
summary(test_model3)$coefficients[2] ## fixed slope (beta) value

### new PVAL estimation
## based on https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode
ranova(test_model3)[["Pr(>Chisq)"]][[2]] ## pval from ANOVA-like table for random-effects (lmerTest package)  <---- use this PVAL
anova(test_model3) ## this give you a F-statistic, but no degree of freedom to test. Can't use.
VarCorr(test_model3) ## Variance Estimates
ranef(test_model3) ## random effect intercepts
confint(test_model3) ## confidence intervals
confint(test_model3,parm="PC1",method="boot")  ## bootstrap estimation of confidence intervals ## WARNING: this also takes a long time to compute!!!!
### Add CONFINT to table
#as.data.frame(VarCorr(test_model3))$vcov[1] / (as.data.frame(VarCorr(test_model3))$vcov[1] + as.data.frame(VarCorr(test_model3))$vcov[2]) ## Variance explained by cell_type
## plot test
print("Plot 1, gene test")
par(mar=c(5,4,4,2))
plot(gene3 ~PC1, col=as.factor(master_table$cell_type), pch=16, cex=1,data=master_table)
abline(a = summary(test_model3)$coefficients[1], b = summary(test_model3)$coefficients[2], lwd=1)
text(gene3 ~PC1, labels=cell_type,data=master_table,cex=0.5, font=1, pos=4)

#UNTIL HERE
#### ---- Volcano plots ---- ####
library("plotly")
library("ggplot2")
library("gridExtra")
library(ggrepel)
library(grid)
input_file <- paste(output_dir,"/mouse_beta_innateness/beta_table_1000mvg_PC1_b2.tsv",sep ="")
diff_df <- read.delim(file = input_file,header = TRUE,sep = '\t')
colnames(diff_df)
print("beta_table read")
diff_df <- diff_df[c("gene","beta","pval")]

	#### --- volcano plot all beta, all pval --- ####
p<-ggplot(diff_df,aes(x=(beta),y=-log10(pval),text=gene)) + geom_point()+
  ggtitle("1000mvg, PC1, no ILCs, abs(beta>2), all pval") +
  xlab(paste0("Innateness score")) +
  ylab(paste0("-log10(pval)"))
p

	#### --- sub boxplot all beta table --- ####
innateness_table_score<-read.table(paste(output_dir,"/export/innateness_python.tsv",sep=""), header=TRUE, sep = "\t")
a<-ggplot(innateness_table_score, aes(x = reorder(cell_type,-innateness_score), y = innateness_score, fill=as.factor(cell_type)))+
  scale_y_reverse()+
  geom_boxplot() +
  ggtitle("Distribution abs(beta scores>2)(1000mvg, PC1)") +
  xlab(paste0("Cell type")) +
  ylab(paste0("Innateness score"))+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")
a
b<-ggplot(innateness_table_score, aes(x = reorder(cell_type,-innateness_score), y = innateness_score, fill=as.factor(cell_type)))+
  scale_y_reverse()+
  geom_violin() +
  ggtitle("Distribution abs(beta scores>2)(1000mvg, PC1)") +
  xlab(paste0("Cell type")) +
  ylab(paste0("Innateness score"))+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")
b

data<-diff_df[which(diff_df['pval'] != 0 & abs(diff_df['beta']) > 2 ),c("gene","beta","pval")]
my_tail<-data%>%arrange(beta)%>%tail(20)
my_head<-data%>%arrange(beta)%>%head(20)
grid.newpage()
grid.table(my_head)
grid.newpage()
grid.table(my_tail)
top_peaks <- data[with(data, order(beta, pval)),][1:15,]
top_peaks <- rbind(top_peaks, data[with(data, order(-beta, pval)),][1:15,])

p + geom_label_repel(data = top_peaks,aes(x=(beta),y=-log10(pval),label=gene),
                     box.padding   = 0.08,
                     point.padding = 0.3,
                     segment.color = 'grey50') + theme_classic()

dev.off()
