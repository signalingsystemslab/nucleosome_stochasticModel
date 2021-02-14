# analysis redo for nucleosome stochastic modeling paper

# 4/14/2020

library(DESeq2);library(broom);library(tidyverse);library(dplyr);library(tidyr);library(reshape2);
library(devtools);library(pheatmap);library(ksheu.library1);library(vegan);library(readxl);library(cqn)
library(CCA); library(CCP);library(ggplot2);library(GGally);library(RColorBrewer);library(RRHO);library(R.matlab)
library(clipr);library(matrixStats);library(forcats)

#########################################################################################
#---------------------------------------------atac analysis ---------
setwd("F:/atac_bmdm/4stochastic")

#these are 2 reps of IkbaM/M bmdms atac seq (peaks from reps were merged, then reads counted)
rep1 = read.delim("multicov_merged_biorep1_count.txt", header = F)
rep2 = read.delim("multicov_merged_biorep2_count.txt", header = F)
test = data.frame(peaksize = rep1$V3-rep1$V2)
ggplot(test, aes(log10(peaksize)))+geom_histogram(bins=100)+theme_bw(base_size = 18)
test = cbind(rep1[, c(1:3,4:8, 9:13)], rep2[,c(4:11)]) #all samples including 24hrs
colnames(test) = c("chr", "start", "end", "M0_1", "M2_1","M4_1", "M8_1","M24_1", "W0_1",
                   "W4_1", "W2_1","W8_1", "W24_1","M0_2", "M2_2","M4_2", "M8_2", "W0_2",
                   "W2_2", "W4_2","W8_2")  #reswap the sample swap of "W2_1", "W4_1"

#continue with the version without 24hrs----
test = cbind(rep1[, c(1:3,4:7, 9:12)], rep2[,c(4:11)])
colnames(test) = c("chr", "start", "end", "M0_1", "M2_1","M4_1", "M8_1", "W0_1",
                  "W4_1", "W2_1","W8_1", "M0_2", "M2_2","M4_2", "M8_2", "W0_2",
                  "W2_2", "W4_2","W8_2") #reswap the sample swap of "W2_1", "W4_1" 
test = test[,  c("chr", "start", "end", "M0_1", "M2_1","M4_1", "M8_1", "W0_1",
                 "W2_1","W4_1","W8_1", "M0_2", "M2_2","M4_2", "M8_2", "W0_2",
                 "W2_2", "W4_2","W8_2")]

cts = test[,-c(1:3)] #normalize only
dds <- DESeqDataSetFromMatrix(countData = cts, colData = data.frame(name=colnames(cts), type = rep("null",16)), design = ~1)
dds <- estimateSizeFactors(dds)

#########################################################
#do PCA----
vsd = vst(dds)
mat <- assay(vsd)
rownames(mat) = paste0("Peak_", seq(1:nrow(mat)))
write.table(cbind(Peak = rownames(mat), mat), "atac_allpeaks_VST_bmdmikbamutant_counts.txt", sep = "\t", row.names = F, quote = F)
PCA_from_file("atac_allpeaks_VST_bmdmikbamutant_counts.txt", center = T, scale = T)
plot_pca3("atac_allpeaks_VST_bmdmikbamutant_counts_prcomp_scores.txt", coldata$name, coldata$type, coldata$time, PCx = "PC1", PCy = "PC2")
mat <- limma::removeBatchEffect(mat, vsd$batch)
assay(vsd) <- mat
plotPCA(vsd)


# get all counts----
all <- data.frame(counts(dds, normalized=TRUE))
rownames(all) = paste0("Peak_", seq(1:nrow(cts)))

#########################################################
# get induced peaks----
cpm = all

# make FC  after averaging 2 reps
cpmFC <- cbind( ((cpm[,2:4] + 1) / (cpm.avg[,1] + 1)), (cpm[,c(6:8)]+1)/(cpm.avg[,5]+1),
                ((cpm[,10:12] + 1) / (cpm[,9] + 1)), (cpm[,c(14:16)]+1)/(cpm[,13]+1))
temp <- ifelse(cpmFC >=2, 1, 0) # threshold of 4-fold as cutoff, log2FC of 1
summary(rowSums(temp)>=3) # meet induction threshold in at least two conditions

# 4514 peaks. also meet max cpm >= 1 and basal cpm < 2
cpm.avg = data.frame(cbind(WT0 = rowMeans(cpm[,c(5,13)]), WT2 = rowMeans(cpm[,c(6,14)]), WT4 = rowMeans(cpm[,c(7,15)]),WT8 = rowMeans(cpm[,c(8,16)]),
                           MT0 =rowMeans(cpm[,c(1,9)]), MT2 = rowMeans(cpm[,c(2,10)]), MT4 = rowMeans(cpm[,c(3,11)]), MT8 = rowMeans(cpm[,c(4,12)]) ))
induced <- rowSums(temp)>=1 & cpm.avg[,1] < 2 & cpm.avg[,4] < 2 
summary(induced)
induced <- as.data.frame(cpm.avg[induced,]) # make table of induced peaks (cpm)
pheatmap(na.omit(induced), cluster_cols = F, show_rownames = F, scale = "row")

#make kmeans heatmap----
# scale the cpm dataframe
induced <- induced[order(rownames(induced)),]
scaled <- t(scale(t(induced)))

set.seed(1)
k3b <- kmeans(scaled, 2, nstart = 25)

# order data by cluster
ord_data <- scaled[order(k3b$cluster, decreasing = F), ]

# make clustering annotation
annot_r <- data.frame(row.names = rownames(scaled), cluster = factor(k3b$cluster))
table(k3b$cluster)
num <- table(k3b$cluster)

# make heatmap
pheatmap(ord_data, annotation_row = annot_r,cluster_rows = F, cluster_cols = F, show_rownames = F, 
         # gaps_row = c(num[1], num[1]+ num[2]), 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-2,seq(-1.25,1.25,length=100),2))
pheatmap(ord_data[c(1:num[1]),], annotation_row = annot_r,cluster_rows = F, cluster_cols = F, show_rownames = F, 
         # gaps_row = c(num[1], num[1]+ num[2]), 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-2,seq(-1.25,1.25,length=100),2))

##########################################################
# plot FC
inducedFC <- cpm.avg[rownames(scaled),]
# order data by cluster
ord_data <- inducedFC[order(k3b$cluster, decreasing = F), ]

# make clustering annotation
annot_r <- data.frame(row.names = rownames(scaled), cluster = factor(k3b$cluster))
table(k3b$cluster)
num <- table(k3b$cluster)

# make heatmap
pheatmap(ord_data,cluster_rows = F, cluster_cols = F, show_rownames = F, gaps_row = c(num[1]), 
         colorRampPalette(c("white", "midnightblue"))(100))
pheatmap(ord_data,cluster_rows = F, cluster_cols = F, show_rownames = F, gaps_row = c(num[1]), 
         colorRampPalette(c("white", "darkgreen"))(100), breaks=c(0,seq(1, 20,length=100),35))


##########################################################
#using DESeq2----
# get coldata
coldata = read.delim("F:/atac_bmdm/4stochastic/atac_bmdm_coldata.txt")
rownames(coldata) = coldata$name
coldata$time = as.factor(coldata$time)

# run comparison
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ (type) + (time) + (type):(time))
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds, test = "LRT", reduced = ~type)
tail(colData(dds))
resultsNames(dds) # lists the coefficients
res <- results(dds, name = c("type_WT_vs_mutant"))
res
summary(res)
# plotMA(res, ylim = c(-5, 5))
head(res[order(res$padj),], 4)

res2 <- results(dds, name = c("time_4_vs_0"))
res2
summary(res2)

#plot log2 fold changes histogram, and as heatmap
log2fc = res$log2FoldChange[res$log2FoldChange>0 & res$pvalue<0.1]
log2fc = res2$log2FoldChange[res2$log2FoldChange>0 & res2$pvalue<0.1]
hist(log2fc, breaks = 25)
summary(log2fc)
write.csv(log2fc, "F://enhancer_stochastic/ATACseq_mutant_2hr_log2fc.csv", quote = F, row.names = F)

###########################################################################
# heatmap of counts
counts = data.frame(cbind(WT0 = rowMeans(cpm[,c(5,13)]), WT2 = rowMeans(cpm[,c(6,14)]), WT4 = rowMeans(cpm[,c(7,15)]),
               MT0 =rowMeans(cpm[,c(1,9)]), MT2 = rowMeans(cpm[,c(2,10)]), MT4 = rowMeans(cpm[,c(3,11)]) ))

# mat = counts[res2$log2FoldChange>0 & res2$pvalue<0.1 & (rownames(counts) %in% rownames(induced)), ]
mat = counts[res$log2FoldChange>0 & res$pvalue<0.05 & counts[,1] < 5 & counts[,4] < 5 , ]

pheatmap(mat, 
         cluster_cols = F,show_rownames = F, show_colnames = T, clustering_method = "ward.D2",
         scale = "row", gaps_col = c(3),
         # colorRampPalette(c("blue", "white", "red"))(100),
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-4,seq(-1,1.5,length=100),4),
         border_color=NA) 
pheatmap(mat, 
         cluster_cols = F,show_rownames = F, show_colnames = T, clustering_method = "ward.D2",
         scale = "none", gaps_col = c(3),
         colorRampPalette(c("blue", "white", "red"))(100),
         # colorRampPalette((brewer.pal(n = 9, name ="Greens"))[2:9])(103),
         # breaks=c(0,seq(1,7,length=100),8),
         border_color=NA) 
do_kmeans_clustering(mat, k_clusters = 3, cluster_cols = F, show_rownames = F, show_colnames = T,
                     colseps = c(3),breaks = c(-3, seq(-1.75, 1.75, length = 100), 3) )



# -------------------------------------get proportions of max----------------------
peaks = data.frame(res[res$log2FoldChange>0&res$pvalue<0.01,])
pheatmap(cts[rownames(cts) %in% rownames(peaks),], cluster_cols = F,
         scale = "row", gaps_col = c(4,8,12),
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103), show_rownames = F,
         breaks=c(-2,seq(-1,1,length=100),2), border_color=NA) 
percentage = data.frame(cts[rownames(cts) %in% rownames(peaks), c(2,10)]) #at 2hr
# percentage = data.frame(cts[rownames(cts) %in% rownames(peaks), c(3,11)]) #at 4hr

#get top 1% of values mean
# n <- 1
# top = mean(rowMeans(percentage)[rowMeans(percentage) > quantile(rowMeans(percentage),prob=1-n/100)])

percentage = log2(percentage)/max(log2(percentage))
hist(rowMeans(percentage), breaks = 15)
write.csv(rowMeans(percentage), "F://enhancer_stochastic/ATACseq_mutant_2hr_distrib_proportionOfMax.csv", quote = F, row.names = F)


# ----------------------------------------plot a couple single cell trajectories------------------


x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories/ikbamut_10ngTNF.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories/tnf_10ng.mat")


mat.1 = x.1$trajectories
mat.1 = mat.1-mat.1[,1]
mat.1[which(mat.1 <0)] <-0
mat.1 = data.frame(mat.1); rownames(mat.1) = seq(1:nrow(mat.1))
colnames(mat.1) = seq(0, (ncol(mat.1)-1)*5, 5)
order = rownames(mat.1)[order(apply(mat.1[,c(1:20)], 1, max), decreasing = T)] #max within first 100 mins
mat.1 = mat.1[match(order, rownames(mat.1)), ]
#Reshape the data for ggplot
mat.1$cellid = seq(1:nrow(mat.1))
mat.1.m <- reshape2::melt(mat.1[,c(1:48, ncol(mat.1))],id=c("cellid"))
head(mat.1.m)
mat.1.m$type = "Mut"
#ggplot it...
ggplot(mat.1.m[grepl("^25$|^2$", mat.1.m$cellid),], aes(x = as.numeric(variable), y = value, color = as.factor(cellid)))+
  geom_line(size = 3)+theme_bw()


mat.2 = x.2$data
mat.2 = mat.2-mat.2[,1]
mat.2[which(mat.2 <0)] <-0
mat.2 = data.frame(mat.2); rownames(mat.2) = seq(1:nrow(mat.2))
colnames(mat.2) = seq(0, (ncol(mat.2)-1)*5, 5)
order = rownames(mat.2)[order(apply(mat.2[,c(1:20)], 1, max), decreasing = T)] #max within first 100 mins
mat.2 = mat.2[match(order, rownames(mat.2)), ]
#Reshape the data for ggplot
mat.2$cellid = seq(1:nrow(mat.2))
mat.2.m <- reshape2::melt(mat.2[,c(1:48, ncol(mat.2))],id=c("cellid"))
head(mat.2.m)
mat.2.m$type = "WT"
ggplot(mat.2.m[grepl("^25$|^2$", mat.2.m$cellid),], aes(x = as.numeric(variable), y = value, color = as.factor(cellid)))+
  geom_line(size = 3)+theme_bw()

tmp = rbind(mat.1.m, mat.2.m)
ggplot(tmp[grepl("^25$", tmp$cellid),], aes(x = as.numeric(variable)*5, y = value, color = (type)))+
  geom_line(size = 3)+theme_bw(base_size = 25)+xlab("Minutes")+scale_x_continuous(breaks=c(0,60,120,180,240))
write.table(tmp[grepl("^25$", tmp$cellid),], "F://enhancer_stochastic/two_sample_traces.txt", quote=F,row.names = F, sep = "\t")

pheatmap(mat.1[,c(1:48)], cluster_cols = F, cluster_rows = F,show_rownames = F,show_colnames = F,scale = "none", main = "IkBa M/M BMDM 10ngTNF",clustering_method = "ward.D2",
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
         breaks=c(-.5,seq(0,6,length=100),10))
pheatmap(mat.2[,c(1:48)], cluster_cols = F, cluster_rows = F,show_rownames = F,show_colnames = F,scale = "none", main = "IkBa M/M BMDM 10ngTNF",clustering_method = "ward.D2",
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
         breaks=c(-.5,seq(0,6,length=100),10))

