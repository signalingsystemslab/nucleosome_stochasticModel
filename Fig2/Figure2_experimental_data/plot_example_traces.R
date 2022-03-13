# ----------------------------------------plot a couple single cell trajectories------------------


x.1 = readMat("./Fig2/Figure2_experimental_data/ikbamut_10ngTNF.mat")
x.2 = readMat("./Fig2/Figure2_experimental_data/tnf_10ng.mat")


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
# write.table(tmp[grepl("^25$", tmp$cellid),], "./Fig2/Figure2_experimental_data/two_sample_traces.txt", quote=F,row.names = F, sep = "\t")

pheatmap(mat.1[,c(1:48)], cluster_cols = F, cluster_rows = F,show_rownames = F,show_colnames = F,scale = "none", main = "IkBa M/M BMDM 10ngTNF",clustering_method = "ward.D2",
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
         breaks=c(-.5,seq(0,6,length=100),10))
pheatmap(mat.2[,c(1:48)], cluster_cols = F, cluster_rows = F,show_rownames = F,show_colnames = F,scale = "none", main = "IkBa M/M BMDM 10ngTNF",clustering_method = "ward.D2",
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
         breaks=c(-.5,seq(0,6,length=100),10))

