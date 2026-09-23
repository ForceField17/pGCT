# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(DESeq2)
library(stringr)
library("pheatmap")
library("RColorBrewer")
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(preprocessCore)
library(clusterProfiler)
library(enrichplot)
library(ggbeeswarm)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#raw data preprocessing
Sample_features <- read.table("../../TGCT_subtyping/Info_temp.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- gsub("-",".",samples$AID)
samples$AID <- rownames(samples) 
#samples$Histology <- samples$Dominant.His
################
counts1 <- read.table("matrix_TCGA-TGCT.txt", head=T, row.names = "A_CpG_ID")
nrow(counts1)

samples <- samples[which(samples$AID %in% colnames(counts1)),]
counts1 <- counts1[, samples$AID]
samples <- samples[which(samples$AID %in% colnames(counts1)),]




keep <- rowSums(!is.na(counts1)) >= ncol(counts1) 

nrow(counts1)
expMatrix <- counts1[keep,]
nrow(expMatrix)


#PCA
library(genefilter)
library("sva")
library(ggrepel)

rv <- rowVars(expMatrix)
ntop <- 1000000
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,  length(rv)))]
pca <- prcomp(t(expMatrix[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup <- c("AID","Histology","Dominant.His","Age")
if (!all(intgroup %in% colnames(samples))) {
  stop("the argument 'intgroup' should specify columns of colData(dds)")
}
intgroup.df <- as.data.frame(samples[, intgroup, drop = FALSE])
group <- if (length(intgroup) > 1){ factor(apply(intgroup.df, 1, paste, collapse = " : ")) }else{ samples[[intgroup]]}

d <- data.frame(pca$x,intgroup.df)
xx <- d[which(d$Dominant.His=="YST"),]
myPCA <- ggplot()+theme_classic()
myPCA <- myPCA + geom_point(data = d, aes(x =PC1, y = PC2,color=Histology),  shape = 21,alpha=0.7,size = 1.5,stroke=0.7) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) #+ coord_fixed() 
myPCA <- myPCA + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=12,face='plain',color='black'),legend.key.width=unit(1,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='right',
                       legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                       axis.text.x=element_text(size=12,face='plain',color='black'),axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
myPCA <- myPCA + geom_text_repel(data = xx, aes(x =PC1, y = PC2,  label=AID),size=1) +
  scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))
figure_1<-rbind(ggplotGrob(myPCA ),size="first")
ggsave(file="./figs/PCA_Beta_new.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)
########



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#######################33

L1 <- read.table("../450k_L1pro_list.bed",header = F)
L1_list <- unique(L1$V5[which(L1$V1 %in% rownames(expMatrix))])

P_YST_vs_G <- rep(NA,length(L1_list))
Diff_YST_vs_G <- rep(NA,length(L1_list))
P_YST_vs_O <- rep(NA,length(L1_list))
Diff_YST_vs_O <- rep(NA,length(L1_list))
P_YST_vs_A <- rep(NA,length(L1_list))
Diff_YST_vs_A <- rep(NA,length(L1_list))
for(i in c(1:length(L1_list))){
  L1HS <- L1[which(L1$V5==L1_list[i]),]
  myTable <- counts1[which(rownames(counts1) %in% L1HS$V1),]
  L1HSmethy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
  colnames(L1HSmethy) <- c("AID","methy")
  myTable1 <- merge(L1HSmethy,samples)
  myTable1$His <- myTable1$Histology
  myTable1$His[which(myTable1$Histology!="YST") ] <- "aOther"
  myTable1$His[which(myTable1$Histology=="Germinoma") ] <- "Germinoma"
  G1 <- myTable1$methy[which(myTable1$His=="aOther")]
  G2 <- myTable1$methy[which(myTable1$His=="Germinoma")]
  G3 <- myTable1$methy[which(myTable1$His=="YST")]
  G0 <- myTable1$methy[which(myTable1$His!="YST")]
  P_YST_vs_G[i] <- wilcox.test(G3,G2)$p.value
  P_YST_vs_O[i] <- wilcox.test(G3,G1)$p.value
  Diff_YST_vs_G[i] <- mean(G3) - mean(G2)
  Diff_YST_vs_O[i] <- mean(G3) - mean(G1)
  P_YST_vs_A[i] <- wilcox.test(G3,G0)$p.value
  Diff_YST_vs_A[i] <- mean(G3) - mean(G0)
}

final_data  <- data.frame(L1_list,Diff_YST_vs_G,Diff_YST_vs_O,P_YST_vs_G,P_YST_vs_O,P_YST_vs_A,Diff_YST_vs_A)
############################################
final_data$FDR_YST_vs_G <- p.adjust(final_data$P_YST_vs_G,method = "fdr")
final_data$FDR_YST_vs_O <- p.adjust(final_data$P_YST_vs_O,method = "fdr")
final_data$FDR_YST_vs_A <- p.adjust(final_data$P_YST_vs_A,method = "fdr",n = nrow(final_data))

resultsCNV <- final_data 

library(ggrepel)
myGene <- c("L1HS")
name <- resultsCNV[which(resultsCNV$L1_list %in% myGene),]

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_point(data = resultsCNV , aes(x = Diff_YST_vs_G, y = -log10(FDR_YST_vs_G)),fill="grey90",color="grey10",alpha=0.7,size=2, shape = 21,stroke=0.5) 
NAN_plot <- NAN_plot + geom_hline(aes(yintercept= -log10(0.05) ),color='black',size=.5,linetype='dashed') 
NAN_plot <- NAN_plot + geom_vline(aes(xintercept=0),color='black',size=.5,linetype=2)
NAN_plot <- NAN_plot + geom_point(data = name, aes(x = Diff_YST_vs_G, y = -log10(FDR_YST_vs_G)),fill="#c51b8a",size=2,shape=21,alpha=0.7,stroke=0.5)#,fill="#009BFF") 

NAN_plot <- NAN_plot + geom_text_repel(data = resultsCNV  , aes(x = Diff_YST_vs_G, y = -log10(FDR_YST_vs_G),label=L1_list),
                                       segment.color = "black",size =2,color ='black',min.segment.length = 0.4 , box.padding = 0.2)

NAN_plot <- NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,2,1,2),'lines'),
                             plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                             legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=12,face='plain'),legend.title=element_text(size=12,face='plain'),
                             axis.text.x=element_text(size=14,face='plain',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),
                             axis.title.x=element_text(size=16,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-0.28,0.38),breaks=seq(-80,80,0.1))+ scale_y_continuous(expand=c(0,0),limits=c(-0.01,2.8),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('Difference of average methylation beta value\n[YST - Germinoma]')+ylab('-log10(FDR)')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_YST_vs_G_fdr.pdf", plot=AllTumorure_2,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)


NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_point(data = resultsCNV , aes(x = Diff_YST_vs_G, y = -log10(P_YST_vs_G)),fill="grey90",color="grey10",alpha=0.7,size=2, shape = 21,stroke=0.5) 
NAN_plot <- NAN_plot + geom_hline(aes(yintercept= -log10(0.05/(nrow(final_data)/(sum(final_data$FDR_YST_vs_G < 0.05)+1)) ) ),color='black',size=.5,linetype='dashed') 
NAN_plot <- NAN_plot + geom_vline(aes(xintercept=0),color='black',size=.5,linetype=2)
NAN_plot <- NAN_plot + geom_point(data = name, aes(x = Diff_YST_vs_G, y = -log10(P_YST_vs_G)),fill="#c51b8a",size=2,shape=21,alpha=0.7,stroke=0.5)#,fill="#009BFF") 

NAN_plot <- NAN_plot + geom_text_repel(data = resultsCNV  , aes(x = Diff_YST_vs_G, y = -log10(P_YST_vs_G),label=L1_list),
                                       segment.color = "black",size =2,color ='black',min.segment.length = 0.4 , box.padding = 0.2)

NAN_plot <- NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,2,1,2),'lines'),
                             plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                             legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=12,face='plain'),legend.title=element_text(size=12,face='plain'),
                             axis.text.x=element_text(size=14,face='plain',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),
                             axis.title.x=element_text(size=16,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-0.28,0.38),breaks=seq(-80,80,0.1))+ scale_y_continuous(expand=c(0,0),limits=c(-0.01,3.8),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('Difference of average methylation beta value\n[YST - Germinoma]')+ylab('-log10(P)')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_YST_vs_G.pdf", plot=AllTumorure_2,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)



##########
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_point(data = resultsCNV , aes(x = Diff_YST_vs_O, y = -log10(P_YST_vs_O)),fill="grey90",color="grey10",alpha=0.7,size=2, shape = 21,stroke=0.5) 
NAN_plot <- NAN_plot + geom_hline(aes(yintercept= -log10(0.05/(nrow(final_data)/(sum(final_data$FDR_YST_vs_O < 0.05)+1)) ) ),color='black',size=.5,linetype='dashed') 
NAN_plot <- NAN_plot + geom_vline(aes(xintercept=0),color='black',size=.5,linetype=2)
NAN_plot <- NAN_plot + geom_point(data = name, aes(x = Diff_YST_vs_O, y = -log10(P_YST_vs_O)),fill="#c51b8a",size=2,shape=21,alpha=0.7,stroke=0.5)#,fill="#009BFF") 

NAN_plot <- NAN_plot + geom_text_repel(data = resultsCNV  , aes(x = Diff_YST_vs_O, y = -log10(P_YST_vs_O),label=L1_list),
                                       segment.color = "black",size =2,color ='black',min.segment.length = 0.4 , box.padding = 0.2)

NAN_plot <- NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,2,1,2),'lines'),
                             plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                             legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=12,face='plain'),legend.title=element_text(size=12,face='plain'),
                             axis.text.x=element_text(size=14,face='plain',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),
                             axis.title.x=element_text(size=16,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-0.38,0.28),breaks=seq(-80,80,0.1))+ scale_y_continuous(expand=c(0,0),limits=c(-0.01,3.8),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('Difference of average methylation beta value\n[YST - otherGCT]')+ylab('P-value')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_YST_vs_O.pdf", plot=AllTumorure_2,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)




##########
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_point(data = resultsCNV , aes(x = Diff_YST_vs_A, y = -log10(P_YST_vs_A)),fill="grey90",color="grey10",alpha=0.7,size=2, shape = 21,stroke=0.5) 
NAN_plot <- NAN_plot + geom_hline(aes(yintercept= -log10(0.05/(nrow(final_data)/(sum(final_data$FDR_YST_vs_A < 0.05)+1)) ) ),color='black',size=.5,linetype='dashed') 
NAN_plot <- NAN_plot + geom_vline(aes(xintercept=0),color='black',size=.5,linetype=2)
NAN_plot <- NAN_plot + geom_point(data = name, aes(x = Diff_YST_vs_A, y = -log10(P_YST_vs_A)),fill="#c51b8a",size=2,shape=21,alpha=0.7,stroke=0.5)#,fill="#009BFF") 

NAN_plot <- NAN_plot + geom_text_repel(data = resultsCNV  , aes(x = Diff_YST_vs_A, y = -log10(P_YST_vs_A),label=L1_list),
                                       segment.color = "black",size =2,color ='black',min.segment.length = 0.4 , box.padding = 0.2)

NAN_plot <- NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,2,1,2),'lines'),
                             plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                             legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=12,face='plain'),legend.title=element_text(size=12,face='plain'),
                             axis.text.x=element_text(size=14,face='plain',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),
                             axis.title.x=element_text(size=16,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-0.38,0.28),breaks=seq(-80,80,0.1))+ scale_y_continuous(expand=c(0,0),limits=c(-0.01,3.8),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('Difference of average methylation beta value\n[YST - nonYST]')+ylab('P-value')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_YST_vs_A.pdf", plot=AllTumorure_2,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)

