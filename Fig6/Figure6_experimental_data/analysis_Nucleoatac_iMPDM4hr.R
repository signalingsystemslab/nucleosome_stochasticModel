
library(ggpubr);library(ksheu.library1);library(dplyr);library(MASS);library(ISLR);library(ggplot2);library(R.matlab);library(pheatmap)
library(RColorBrewer);library(reshape);library(tidyr);library(broom);library(dplyr);library(tidyverse);library(DESeq2)
library(R.matlab);library(DESeq2);library(edgeR);library(NucleoATACR)



#NUCLEOATAC peak width analysis----
# setwd("D://atac_bmdm/atac_PE/out_AH_08202019/NucleoATAC/03_nucposbed/")
# setwd("F://atac_bmdm_PE/NucleoATAC/03_nucposbed/")
# bed1 = read.delim("NucleoATAC_4KO-0h.nucpos.bed", header = F)
# bed2 = read.delim("NucleoATAC_4KO-4h.nucpos.bed", header = F)

# setwd("F://enhancer_stochastic/ATAC_iMPDMs_4hrs/03_nucposbed/")
setwd("./Fig6/Figure6_experimental_data/NucleoATAC_outputs")
bed1 = read.delim("NucleoATAC_Unstim1.nucpos.bed", header = F)
bed2 = read.delim("NucleoATAC_LPS.nucpos.bed", header = F)


hist(bed1$V13)
hist(bed2$V13)

hist(bed1$V5)
hist(bed2$V5)

###################################################
#TSS centered analysis and combinations----
motif1 = "irf3"#"nfkb" ,"nfkb-5half", "isre", "ap1", "p53", "ascl1", "irf1", "stat1", "pu1"
motif2 = "irf3"#"nfkb" ,"nfkb-5half", "isre", "ap1", "p53", "ascl1", "irf1", "stat1", "pu1"
stim = "Unstim1"
hist.motif1 = read.delim(paste0("NucleoATAC_",stim,".nucpos.bed.", motif1,".motif.size400.nohist.txt"))
colnames(hist.motif1)[ncol(hist.motif1)] = "distance2center"

hist.motif2 = read.delim(paste0("NucleoATAC_",stim,".nucpos.bed.", motif2,".motif.size400.nohist.txt"))
colnames(hist.motif2)[ncol(hist.motif2)] = "distance2center"
hist.all = cbind(hist.motif1, motif2 = hist.motif2[, ncol(hist.motif1)])
colnames(hist.all)[ncol(hist.all)] = "distance2center"

hist.motif1$distance2center = gsub("\\s*\\([^\\)]+\\)","",as.character(hist.motif1$distance2center))
hist.motif2$distance2center = gsub("\\s*\\([^\\)]+\\)","",as.character(hist.motif2$distance2center))

hist1.numbers = unlist(str_split(hist.motif1$distance2center, ","));hist1.numbers = as.numeric(hist1.numbers[hist1.numbers!= ""])
hist2.numbers = unlist(str_split(hist.motif2$distance2center, ","));hist2.numbers = as.numeric(hist2.numbers[hist2.numbers!= ""])
hist(hist1.numbers)
hist(hist2.numbers)

frame = data.frame(position = c(hist1.numbers, hist2.numbers), 
                   motif = c(rep("nfkb",length(hist1.numbers)),rep("oct4",length(hist2.numbers)) ))
ggplot(frame, aes(position, color = motif, fill = motif))+geom_density(size = 2, alpha = 0.3)+theme_bw(base_size = 14)


#NUCLEOATAC LPS motif analysis----
# setwd("F://atac_bmdm_PE/NucleoATAC/03_nucposbed/")
# setwd("F://enhancer_stochastic/ATAC_iMPDMs_4hrs/03_nucposbed/")
setwd("./Fig6/Figure6_experimental_data/NucleoATAC_outputs")

if (1){
motif = "irf3"#"nfkb-all3" #"nfkb" ,"nfkb-5half", "isre", "ap1", "p53", "ascl1", "irf1", "stat1", "pu1"
stim = "LPS"

hist1 = read.delim(paste0("NucleoATAC_Unstim1.nucpos.bed.", motif,".motif.size400.nohist.txt"))
hist2 = read.delim(paste0("NucleoATAC_",stim,".nucpos.bed.", motif,".motif.size400.nohist.txt"))

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

ggplot(data=hist1.numbers.frame, aes(hist1.numbers.frame$distace2dyad)) + geom_histogram(bins = 35)+ggtitle("Unstim1")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
ggplot(data=hist2.numbers.frame, aes(hist2.numbers.frame$distace2dyad)) + geom_histogram(bins = 35)+ggtitle(paste0(stim,"-4hr"))+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)

#combine
dat <- data.frame(dens = c(hist1.numbers.frame$distace2dyad, hist2.numbers.frame$distace2dyad)
                  , lines = c(rep("0hr",length(hist1.numbers.frame$distace2dyad)), rep("4hr", length(hist2.numbers.frame$distace2dyad))))
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.3)+theme_bw(base_size = 24)+ggtitle(paste0(stim, " - ", motif))
# write.table(dat, "NucleoATAC_density_v1data.txt", quote=F,sep="\t", row.names = F)


##############################################################################
# using matched locations----


hist1$hist2locations = hist2[, ncol(hist2)][match(hist1$Gene.Name, hist2$Gene.Name) ]

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


ggplot(data=hist1.numbers.frame, aes(hist1.numbers.frame$distace2dyad)) + geom_histogram(bins = 20)+ggtitle("LPS-0hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
ggplot(data=hist2.numbers.frame, aes(hist2.numbers.frame$distace2dyad)) + geom_histogram(bins = 20)+ggtitle("LPS-4hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)

#combine
dat <- data.frame(distace2dyad = c(hist1.numbers.frame$distace2dyad, hist2.numbers.frame$distace2dyad)
                  , time = c(rep("0hr",length(hist1.numbers.frame$distace2dyad)), rep("4hr", length(hist2.numbers.frame$distace2dyad))))
ggplot(dat, aes(x = distace2dyad, fill = time)) + 
  geom_density(alpha = 0.3)+
  theme_bw(base_size = 24)+ggtitle(paste0(stim," - ", motif))

ggplot(dat, aes(x = distace2dyad, fill = time)) + 
  geom_histogram(aes(fill = time),alpha=0.3, position = 'identity', bins = 15)+
  theme_bw(base_size = 24)+ggtitle(paste0(stim," - ", motif)) +scale_color_brewer(palette="Dark2")+ scale_fill_brewer(palette="Dark2")

# write.table(dat, "nucleosomes_nfkb-all3_matched_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
# write.table(dat, "nucleosomes_isre_matched_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
# write.table(dat, "nucleosomes_irf3_matched_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)

###########################################################################
# probability of eviction over the binned locations ----
if(1){
  bins = 25
  p2 = ggplot(dat[abs(dat$distace2dyad)<200,], aes(x = distace2dyad, fill = time)) + 
    geom_histogram(aes(fill = time),alpha=0.3, position = 'identity', bins = bins)+
    theme_bw(base_size = 24)+ggtitle(paste0(stim, " - ", motif)) +scale_color_brewer(palette="Dark2")+ scale_fill_brewer(palette="Dark2")
  
  table = ggplot_build(p2)$data[[1]]
  frame = data.frame(bin1 = table$xmin[1:bins],bin2 = table$xmax[1:bins], 
                     div = (table$y[1:bins] / table$y[(bins+1):(bins*2)]),
                     diff = (table$y[1:bins] - table$y[(bins+1):(bins*2)]),
                     prob = (table$y[1:bins]-table$y[(bins+1):(bins*2)])/table$y[1:bins])
  frame$bin = rowMeans(frame[,c(1:2)])
  # write.table(frame, "nucleosomes_nfkb-all3_matched_binned_probabilities_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
  # write.table(frame, "nucleosomes_isre_matched_binned_probabilities_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
  # write.table(frame, "nucleosomes_irf3_matched_binned_probabilities_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
  
  
  frame$logratio = log10(table$y[1:bins])/log10(table$y[(bins+1):(bins*2)])
  frame$logprob = log10(frame$prob)
  
  ggplot(frame, aes(bin, prob))+geom_line(size = 2)+theme_bw(base_size = 20)#+ylim(c(0, max(frame$div)))
  
  ggplot(frame, aes(bin, prob))+geom_point(size = 2)+theme_bw(base_size = 20)+
    # stat_smooth(aes(bin, prob), method = loess, span = 0.1, se = T, level = .8)+
    xlab("distance2dyad")+ggtitle(paste0(stim,"-4hrs-", motif))+ylim(0.6,1)
  
  # ggplot(frame, aes(bin, logprob))+geom_line(size = 2)+theme_bw(base_size = 20)#+ylim(c(0, max(frame$div)))
  # ggplot(frame, aes(bin, div))+geom_density(aes(y =div),size = 2)+theme_bw(base_size = 20)
  
   
  # ggplot(frame, aes(bin, log10(div)))+geom_line(size = 2)+theme_bw(base_size = 20)
  # ggplot(frame, aes(bin, logratio))+geom_line(size = 2)+theme_bw(base_size = 20)
  
  
}













##############################################################################
# loop through all combos----


for (s in c("LPS", "TNF", "IFNb","CpG", "P3C", "PIC")){
  for (m in c("nfkb-all3", "isre", "ap1", "irf3", "irf1", "stat1", "pu1", "oct4")){
    print(s)
    print(m)
    
    if (1){
      motif = m #"nfkb-all3" #"nfkb" ,"nfkb-5half", "isre", "ap1", "p53", "ascl1", "irf1", "stat1", "pu1"
      stim = s
      
      hist1 = read.delim(paste0("NucleoATAC_Unstim1.nucpos.bed.", motif,".motif.size400.nohist.txt"))
      hist2 = read.delim(paste0("NucleoATAC_",stim,".nucpos.bed.", motif,".motif.size400.nohist.txt"))
      
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
    
    ggplot(data=hist1.numbers.frame, aes(hist1.numbers.frame$distace2dyad)) + geom_histogram(bins = 35)+ggtitle("Unstim1")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
    ggplot(data=hist2.numbers.frame, aes(hist2.numbers.frame$distace2dyad)) + geom_histogram(bins = 35)+ggtitle(paste0(stim,"-4hr"))+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
    
    #combine
    dat <- data.frame(dens = c(hist1.numbers.frame$distace2dyad, hist2.numbers.frame$distace2dyad)
                      , lines = c(rep("0hr",length(hist1.numbers.frame$distace2dyad)), rep("4hr", length(hist2.numbers.frame$distace2dyad))))
    ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.3)+theme_bw(base_size = 24)+ggtitle(paste0(stim, " - ", motif))
    # write.table(dat, "NucleoATAC_density_v1data.txt", quote=F,sep="\t", row.names = F)
    
    # using matched locations----
    
    
    hist1$hist2locations = hist2[, ncol(hist2)][match(hist1$Gene.Name, hist2$Gene.Name) ]
    
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
    
    
    ggplot(data=hist1.numbers.frame, aes(hist1.numbers.frame$distace2dyad)) + geom_histogram(bins = 20)+ggtitle("LPS-0hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
    ggplot(data=hist2.numbers.frame, aes(hist2.numbers.frame$distace2dyad)) + geom_histogram(bins = 20)+ggtitle("LPS-4hr")+xlab("NFkB binding - Relative to dyad")+theme_bw(base_size = 24)
    
    #combine
    dat <- data.frame(distace2dyad = c(hist1.numbers.frame$distace2dyad, hist2.numbers.frame$distace2dyad)
                      , time = c(rep("0hr",length(hist1.numbers.frame$distace2dyad)), rep("4hr", length(hist2.numbers.frame$distace2dyad))))
    ggplot(dat, aes(x = distace2dyad, fill = time)) + 
      geom_density(alpha = 0.3)+
      theme_bw(base_size = 24)+ggtitle(paste0(stim," - ", motif))
    
    ggplot(dat, aes(x = distace2dyad, fill = time)) + 
      geom_histogram(aes(fill = time),alpha=0.3, position = 'identity', bins = 15)+
      theme_bw(base_size = 24)+ggtitle(paste0(stim," - ", motif)) +scale_color_brewer(palette="Dark2")+ scale_fill_brewer(palette="Dark2")
    
    # write.table(dat, "nucleosomes_nfkb-all3_matched_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
    # write.table(dat, "nucleosomes_isre_matched_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
    
    
    # probability of eviction over the binned locations ----
    if(1){
      bins = 25
      p2 = ggplot(dat[abs(dat$distace2dyad)<200,], aes(x = distace2dyad, fill = time)) + 
        geom_histogram(aes(fill = time),alpha=0.3, position = 'identity', bins = bins)+
        theme_bw(base_size = 24)+ggtitle(paste0(stim, " - ", motif)) +scale_color_brewer(palette="Dark2")+ scale_fill_brewer(palette="Dark2")
      
      table = ggplot_build(p2)$data[[1]]
      frame = data.frame(bin1 = table$xmin[1:bins],bin2 = table$xmax[1:bins], 
                         div = (table$y[1:bins] / table$y[(bins+1):(bins*2)]),
                         diff = (table$y[1:bins] - table$y[(bins+1):(bins*2)]),
                         prob = (table$y[1:bins]-table$y[(bins+1):(bins*2)])/table$y[1:bins])
      frame$bin = rowMeans(frame[,c(1:2)])
      # write.table(frame, "nucleosomes_nfkb-all3_matched_binned_probabilities_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
      # write.table(frame, "nucleosomes_isre_matched_binned_probabilities_LPS4hrs.txt",quote = F, sep = "\t", row.names = F)
      
      
      frame$logratio = log10(table$y[1:bins])/log10(table$y[(bins+1):(bins*2)])
      frame$logprob = log10(frame$prob)
      
      ggplot(frame, aes(bin, prob))+geom_line(size = 2)+theme_bw(base_size = 20)#+ylim(c(0, max(frame$div)))
      
      ggplot(frame, aes(bin, prob))+geom_point(size = 2)+theme_bw(base_size = 20)+
        # stat_smooth(aes(bin, prob), method = loess, span = 0.1, se = T, level = .8)+
        xlab("distance2dyad")+ggtitle(paste0(stim,"-4hrs-", motif))+ylim(0.6,1)
      
      ggsave(paste0(s,"_4hrs_",m, ".png"))
      
      # ggplot(frame, aes(bin, logprob))+geom_line(size = 2)+theme_bw(base_size = 20)#+ylim(c(0, max(frame$div)))
      # ggplot(frame, aes(bin, div))+geom_density(aes(y =div),size = 2)+theme_bw(base_size = 20)
      
      
      # ggplot(frame, aes(bin, log10(div)))+geom_line(size = 2)+theme_bw(base_size = 20)
      # ggplot(frame, aes(bin, logratio))+geom_line(size = 2)+theme_bw(base_size = 20)
      
      
    }
  
  }
}


