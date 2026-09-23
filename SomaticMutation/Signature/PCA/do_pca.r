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

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

#raw data preprocessing
aaa <- read.table("../AFdistribution/results/Summary_all_SomaticAlter.WGS.txt",header = F)
Sample_features <- read.table("../../../lib/Clinical.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$AID
samples <- samples[which(samples$Seq=="WGS" & (samples$AID %in%  aaa$V1[which(aaa$V3>=200)])),] 
samples$Site[which(samples$Site!="Intracranial")] <- "Extracranial"

SNV <- read.table("../SBS96_matrix_allWGS.txt",header = F)
colnames(SNV) <- c("AID","SBS","Freq")
SBS96 <- unique(SNV$SBS)


MatrixSBS <- matrix(rep(0,length(SBS96) * nrow(samples)), nrow = length(SBS96), ncol = nrow(samples))
for(i in 1:nrow(samples)){
  for(j in 1:length(SBS96)){
      tmp <- SNV$Freq[which(SNV$AID==samples$AID[i] & SNV$SBS==SBS96[j])] 
      if(length(tmp) > 0){
        MatrixSBS[j,i] <- tmp
      }
  }
}
rownames(MatrixSBS) <- SBS96
colnames(MatrixSBS) <- samples$AID

counts1 <- MatrixSBS
keep <- rowSums(counts1 >= 8) >= 5
counts <- counts1#[keep,]
nrow(counts1)
nrow(counts)
counts <- log2(counts+1)

library(preprocessCore)



library(genefilter)
library(ggrepel)



scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

xxx <- t(counts)
xxx <- scale_rows(xxx)
pca <- prcomp(xxx)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup <- c("AID", "Histology","Site")
if (!all(intgroup %in% colnames(samples))) {
  stop("the argument 'intgroup' should specify columns of colData(dds)")
}
intgroup.df <- as.data.frame(samples[, intgroup, drop = FALSE])
group <- if (length(intgroup) > 1){ factor(apply(intgroup.df, 1, paste, collapse = " : ")) }else{ samples[[intgroup]]}

d <- data.frame(pca$x,intgroup.df)

myPCA <- ggplot()+theme_classic()
myPCA <- myPCA + geom_point(data = d, aes(x =PC1, y = PC2,  fill = Histology,shape=Site),alpha=0.75,size = 2.6,stroke=0.3) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
  scale_shape_manual(values = c(22,21))+ 
  scale_fill_manual(values = c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))+
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) #
myPCA <- myPCA + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=12,face='plain',color='black'),legend.key.width=unit(1,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='right',
                       legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                       axis.text.x=element_text(size=12,face='plain',color='black'),axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
myPCA <- myPCA + coord_fixed() + 
  scale_x_continuous(expand=c(0,0),limits = c(-3.8,6.4),breaks = seq(-8,8,2)) + scale_y_continuous(expand=c(0,0),limits = c(-4.8,5.3),breaks = seq(-8,8,2)) 
#myPCA <- myPCA + geom_text_repel(data = d, aes(x =PC1, y = PC2,  label=AID),size=2)
figure_1<-rbind(ggplotGrob(myPCA ),size="first")
ggsave(file="./PCA.pdf", plot=figure_1,bg = 'white', width = 11, height = 8, units = 'cm', dpi = 600)
########
