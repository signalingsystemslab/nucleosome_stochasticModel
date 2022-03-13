
library(ggpubr);library(ksheu.library1);library(dplyr);library(MASS);library(ISLR);library(ggplot2);library(R.matlab);library(pheatmap)
library(RColorBrewer);library(reshape);library(tidyr);library(broom);library(dplyr);library(tidyverse);library(DESeq2)
library(R.matlab);library(DESeq2);library(edgeR);library(NucleoATACR)



#NUCLEOATAC peak width analysis----
# setwd("D://atac_bmdm/atac_PE/out_AH_08202019/NucleoATAC/03_nucposbed/")
# setwd("F://atac_bmdm_PE/NucleoATAC/03_nucposbed/")
setwd("./Fig4/Figure4_experimental_data/NucleoATAC_outputs")
bed1 = read.delim("NucleoATAC_4KO-0h.nucpos.bed", header = F)
bed2 = read.delim("NucleoATAC_4KO-4h.nucpos.bed", header = F)

hist(bed1$V13)
hist(bed2$V13)

hist(bed1$V5)
hist(bed2$V5)

#TSS centered analysis and combinations----
motif1 = "nfkb-all3"#"nfkb" ,"nfkb-5half", "isre", "ap1", "p53", "ascl1", "irf1", "stat1", "pu1"
motif2 = "nfkb-all3"#"nfkb" ,"nfkb-5half", "isre", "ap1", "p53", "ascl1", "irf1", "stat1", "pu1"

hist.motif1 = read.delim(paste0("NucleoATAC_4KO-0h.nucpos.bed.", motif1,".motif.size400.nohist.txt"))
colnames(hist.motif1)[ncol(hist.motif1)] = "distance2center"

hist.motif2 = read.delim(paste0("NucleoATAC_4KO-0h.nucpos.bed.", motif2,".motif.size400.nohist.txt"))
hist.all = cbind(hist.motif1, motif2 = hist.motif2[, ncol(hist.motif1)])

colnames(hist.all)[ncol(hist.all)] = "distance2center"
hist1$distance2center = gsub("\\s*\\([^\\)]+\\)","",as.character(hist1$distance2center))
hist1.numbers = unlist(str_split(hist1$distance2center, ","));hist1.numbers = as.numeric(hist1.numbers[hist1.numbers!= ""])


#NUCLEOATAC 4KO motif analysis----
# setwd("F://atac_bmdm_PE/NucleoATAC/03_nucposbed/")
setwd("./Fig4/Figure4_experimental_data/NucleoATAC_outputs")

if (1){
motif = "nfkb-all3" #"nfkb" ,"nfkb-5half", "isre", "ap1", "p53", "ascl1", "irf1", "stat1", "pu1"

hist1 = read.delim(paste0("NucleoATAC_4KO-0h.nucpos.bed.", motif,".motif.size400.nohist.txt"))
hist2 = read.delim(paste0("NucleoATAC_4KO-4h.nucpos.bed.", motif,".motif.size400.nohist.txt"))

# remove low confidence peaks
# hist1 = hist1[hist1$Peak.Score > .4,]
# hist2 = hist2[hist2$Peak.Score > .4,]
# hist3 = hist3[hist3$Peak.Score > .4,]
# hist4 = hist4[hist4$Peak.Score > .4,]


if(motif =="nfkb-all3"){
  hist1$combine <- gsub(",{2,}", ",", sub("^,+|,+$", "", do.call(paste, c(hist1[22],hist1[23],hist1[24], sep=","))))
  hist2$combine <- gsub(",{2,}", ",", sub("^,+|,+$", "", do.call(paste, c(hist2[22],hist2[23],hist2[24], sep=","))))

}



colnames(hist1)[ncol(hist1)] = "distance2center"
hist1$distance2center = gsub("\\s*\\([^\\)]+\\)","",as.character(hist1$distance2center))
hist1.numbers = unlist(str_split(hist1$distance2center, ","));hist1.numbers = as.numeric(hist1.numbers[hist1.numbers!= ""])
# hist1.numbers = as.numeric(gsub(",..*","",as.character(hist1$distance2center))); hist1.numbers = na.omit(hist1.numbers)
hist1.numbers.frame = data.frame(distace2dyad = hist1.numbers[abs(hist1.numbers)<200])

colnames(hist2)[ncol(hist2)] = "distance2center"
hist2$distance2center = gsub("\\s*\\([^\\)]+\\)","",as.character(hist2$distance2center))
hist2.numbers = unlist(str_split(hist2$distance2center, ","));hist2.numbers = as.numeric(hist2.numbers[hist2.numbers!= ""])
# hist2.numbers = as.numeric(gsub(",..*","",as.character(hist2$distance2center))); hist2.numbers = na.omit(hist2.numbers)
hist2.numbers.frame = data.frame(distace2dyad = hist2.numbers[abs(hist2.numbers)<200])


}

ggplot(data=hist1.numbers.frame, aes(hist1.numbers.frame$distace2dyad)) + geom_histogram(bins = 35)+ggtitle("4KO-0hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
ggplot(data=hist2.numbers.frame, aes(hist2.numbers.frame$distace2dyad)) + geom_histogram(bins = 35)+ggtitle("4KO-4hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)

#combine
dat <- data.frame(dens = c(hist1.numbers.frame$distace2dyad, hist2.numbers.frame$distace2dyad)
                  , lines = c(rep("0hr",length(hist1.numbers.frame$distace2dyad)), rep("4hr", length(hist2.numbers.frame$distace2dyad))))
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.3)+theme_bw(base_size = 24)+ggtitle(paste0("4KO - ", motif))
# write.table(dat, "NucleoATAC_density_v1data.txt", quote=F,sep="\t", row.names = F)


################################################################################
# plot find motif output instead----
# setwd("F://atac_bmdm_PE/NucleoATAC/03_nucposbed/")
# find1 = read.delim("NucleoATAC_4KO-0h.nucpos.bed.find.nfkb.motif.txt")
# find2 = read.delim("NucleoATAC_4KO-4h.nucpos.bed.find.nfkb.motif.txt")
# 
# hist(find1$Offset, breaks = 20)
# hist(find2$Offset, breaks = 20)
# 
# ggplot(data=find1, aes(find1$Offset)) + geom_histogram(bins = 25)+ggtitle("4KO-0hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
# ggplot(data=find2, aes(find2$Offset)) + geom_histogram(bins = 25)+ggtitle("4KO-4hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
# 
# dat <- data.frame(dens = c(find1$Offset, find2$Offset)
#                   , lines = c(rep("0hr",length(find1$Offset)), rep("4hr", length(find2$Offset))))
# ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.3)+theme_bw(base_size = 24)+ggtitle(paste0("4KO - ", motif))
# 


################################################################################
#replotting of binding motifs----
# ggplot(output_table) + geom_segment(aes(x= start_pos, y= row_pos, xend= end_pos, yend= row_pos), size=1) + 
#   # scale_y_reverse(limits = c(263, 0), expand=c(0.0001,0.0001), labels = gene_name, breaks = 1:263) +   ####expand command is done to get rid of white space
#   scale_x_continuous(limits = c(-200, 200), breaks = c(-200, -100, 0, 100, 200), 
#                      labels = c("-0.2kb", "-0.1kb", "dyad", "+0.1kb", "+0.2kb"), expand=c(0,0), position = "top") + 
#   theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
#                      axis.title.x = element_blank(), panel.grid.major = element_blank(), 
#                      panel.grid.minor = element_blank(), axis.text.x = element_text(face = "bold", size = 10), 
#                      axis.line = element_line(color = "black"),axis.text.y = element_text(face = "bold", size = 10)) + 
#   scale_color_manual(values = c("red", "black"))

##############################################################################
# using matched locations----


hist1$hist2locations = hist2[, c(25)][match(hist1$Gene.Name, hist2$Gene.Name) ]

#if not found in hist1, then ignore
hist1$hist2.ignore  = ifelse(hist1$distance2center=="", "", hist1$hist2locations)

#if in hist2location only, add 0
# hist1$distance2center = ifelse((hist1$distance2center=="") & (hist1$hist2locations!="") , runif(1, -200, 200), hist1$distance2center)

#take just the first, or splitting on comma and taking all
# hist1.numbers = as.numeric(gsub(",..*","",as.character(hist1$distance2center))); hist1.numbers = na.omit(hist1.numbers)
hist1.numbers = unlist(str_split(hist1$distance2center, ","));hist1.numbers = as.numeric(hist1.numbers[hist1.numbers!= ""])
hist1.numbers.frame = na.omit(data.frame(distace2dyad = hist1.numbers[abs(hist1.numbers)<200]))

# hist2.numbers = as.numeric(gsub(",..*","",as.character(hist1$hist2locations))); hist2.numbers = na.omit(hist2.numbers)
hist2.numbers = unlist(str_split(hist1$hist2.ignore, ","));hist2.numbers = as.numeric(hist2.numbers[hist2.numbers!= ""])
hist2.numbers.frame = na.omit(data.frame(distace2dyad = hist2.numbers[abs(hist2.numbers)<200]))


ggplot(data=hist1.numbers.frame, aes(hist1.numbers.frame$distace2dyad)) + geom_histogram(bins = 20)+ggtitle("4KO-0hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
ggplot(data=hist2.numbers.frame, aes(hist2.numbers.frame$distace2dyad)) + geom_histogram(bins = 20)+ggtitle("4KO-4hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)

#combine
dat <- data.frame(distace2dyad = c(hist1.numbers.frame$distace2dyad, hist2.numbers.frame$distace2dyad)
                  , time = c(rep("0hr",length(hist1.numbers.frame$distace2dyad)), rep("4hr", length(hist2.numbers.frame$distace2dyad))))
ggplot(dat, aes(x = distace2dyad, fill = time)) + 
  geom_density(alpha = 0.3)+
  theme_bw(base_size = 24)+ggtitle(paste0("4KO - ", motif))

ggplot(dat, aes(x = distace2dyad, fill = time)) + 
  geom_histogram(aes(fill = time),alpha=0.3, position = 'identity', bins = 15)+
  theme_bw(base_size = 24)+ggtitle(paste0("4KO - ", motif)) +scale_color_brewer(palette="Dark2")+ scale_fill_brewer(palette="Dark2")

# write.table(dat, "nucleosomes_nfkb-all3_matched.txt",quote = F, sep = "\t", row.names = F)

###########################################################################
# probability of eviction over the binned locations ----
if(1){
  bins = 25
  p2 = ggplot(dat[abs(dat$distace2dyad)<200,], aes(x = distace2dyad, fill = time)) + 
    geom_histogram(aes(fill = time),alpha=0.3, position = 'identity', bins = bins)+
    theme_bw(base_size = 24)+ggtitle(paste0("4KO - ", motif)) +scale_color_brewer(palette="Dark2")+ scale_fill_brewer(palette="Dark2")
  
  table = ggplot_build(p2)$data[[1]]
  frame = data.frame(bin1 = table$xmin[1:bins],bin2 = table$xmax[1:bins], 
                     div = (table$y[1:bins] / table$y[(bins+1):(bins*2)]),
                     diff = (table$y[1:bins] - table$y[(bins+1):(bins*2)]),
                     prob = (table$y[1:bins]-table$y[(bins+1):(bins*2)])/table$y[1:bins])
  frame$bin = rowMeans(frame[,c(1:2)])
  # write.table(frame, "nucleosomes_nfkb-all3_matched_binned_probabilities.txt",quote = F, sep = "\t", row.names = F)
  
  
  frame$logratio = log10(table$y[1:bins])/log10(table$y[(bins+1):(bins*2)])
  frame$logprob = log10(frame$prob)
  
  ggplot(frame, aes(bin, prob))+geom_line(size = 2)+theme_bw(base_size = 20)#+ylim(c(0, max(frame$div)))
  
  ggplot(frame, aes(bin, prob))+geom_point(size = 2)+theme_bw(base_size = 20)+
    stat_smooth(aes(bin, prob), method = loess, span = 0.1, se = T, level = .8)+xlab("distance2dyad")
  
  # ggplot(frame, aes(bin, logprob))+geom_line(size = 2)+theme_bw(base_size = 20)#+ylim(c(0, max(frame$div)))
  # ggplot(frame, aes(bin, div))+geom_density(aes(y =div),size = 2)+theme_bw(base_size = 20)
  
   
  # ggplot(frame, aes(bin, log10(div)))+geom_line(size = 2)+theme_bw(base_size = 20)
  # ggplot(frame, aes(bin, logratio))+geom_line(size = 2)+theme_bw(base_size = 20)
  
  
}

