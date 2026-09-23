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
#samples$Histology <- samples$RNA_cluster
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
intgroup <- c("AID","Histology","Dominant.His","Age","RNA_cluster","Methy_cluster")
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
ggsave(file="./figs/PCA_Beta.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)
########



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#######################33

L1 <- read.table("../450k_L1pro_list.bed",header = F)


################# All L1
L1HS <- L1#[which(L1$V5=="L1HS"),]
table(L1HS$V3)

myTable <- counts1[which(rownames(counts1) %in% L1HS$V1),]

L1HSmethy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
colnames(L1HSmethy) <- c("AID","methy")


##################
myTable1 <- merge(L1HSmethy,samples)
#myTable1 <- myTable1[which(myTable1$Histology %in% c("YST","IMT","Teratoma","Germinoma")),]
myTable1$Histology <- myTable1$Histology
myTable1$His <- myTable1$Histology
myTable1$His[which(myTable1$Histology!="YST") ] <- "aOther"
myTable1$His[which(myTable1$Histology=="Germinoma") ] <- "Germinoma"

NAN_plot <- ggplot(data=myTable1,aes(x=His,y=methy)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable1,aes(fill=Histology),shape=21,width = 0.3,size=2,alpha=0.75,stroke=0.1, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=median, aes(ymax = ..y.., ymin = ..y..),color="black",geom = "errorbar" ,width = 0.19,size=0.6)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0.2,1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("Avargea beta-value of\nCpG sites in L1pro") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))
NAN_plot <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))

the <- compare_means(methy ~ His,  data = myTable1,paired = F ,method = "wilcox.test")
the 
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 2,
  y.position = c( 0.88,0.92,0.96)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="figs/beta_L1.pdf", plot=figure_2,bg = 'white', width =7, height = 8, units = 'cm', dpi = 600)

tmp <- myTable1[,c("case","methy")]
colnames(tmp) <- c("case","methy_L1pro")
write.table(tmp,file = "methy_L1pro.txt",row.names = F,quote = F,sep = "\t")

###############3
################# All L1
L1HS <- L1[which(L1$V5=="L1PA2"),]
table(L1HS$V3)

myTable <- counts1[which(rownames(counts1) %in% L1HS$V1),]

L1HSmethy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
colnames(L1HSmethy) <- c("AID","methy")


myTable1 <- merge(L1HSmethy,samples)
#myTable1 <- myTable1[which(myTable1$Histology %in% c("YST","IMT","Teratoma","Germinoma")),]
myTable1$Histology <- myTable1$Histology
potential_ysts <- c("TCGA_YU_A94I","TCGA_2G_AAGT","TCGA_2G_AAFV","TCGA_YU_AA4L","TCGA_YU_AA61","TCGA_W4_A7U3","TCGA_2G_AAH4")
#myTable1$Histology[which(myTable1$case %in% potential_ysts)] <- "YST"
myTable1$His <- myTable1$Histology
myTable1$His[which(myTable1$Histology!="YST") ] <- "aOther"
myTable1$His[which(myTable1$Histology=="Germinoma") ] <- "Germinoma"


NAN_plot <- ggplot(data=myTable1,aes(x=His,y=methy)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable1,aes(fill=Histology),shape=21,width = 0.3,size=2,alpha=0.75,stroke=0.1, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=median, aes(ymax = ..y.., ymin = ..y..),color="black",geom = "errorbar" ,width = 0.19,size=0.6)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0.2,1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("Avargea beta-value of\nCpG sites in L1PA2pro") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))
NAN_plot <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))

the <- compare_means(methy ~ His,  data = myTable1,paired = F ,method = "wilcox.test")
the 
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 2,
  y.position = c( 0.88,0.92,0.96)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="figs/beta_L1PA2.pdf", plot=figure_2,bg = 'white', width =7, height = 8, units = 'cm', dpi = 600)

tmp <- myTable1[,c("case","methy")]
colnames(tmp) <- c("case","methy_L1PA2pro")
write.table(tmp,file = "methy_L1PA2pro.txt",row.names = F,quote = F,sep = "\t")

###############3
################# All L1
L1HS <- L1[which(L1$V5=="HAL1"),]
table(L1HS$V3)

myTable <- counts1[which(rownames(counts1) %in% L1HS$V1),]

L1HSmethy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
colnames(L1HSmethy) <- c("AID","methy")


myTable1 <- merge(L1HSmethy,samples)
myTable1$His <- myTable1$Histology
myTable1$His[which(myTable1$Histology!="YST") ] <- "aOther"
myTable1$His[which(myTable1$Histology=="Germinoma") ] <- "Germinoma"

NAN_plot <- ggplot(data=myTable1,aes(x=His,y=methy)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable1,aes(fill=Histology),shape=21,width = 0.3,size=2,alpha=0.75,stroke=0.1, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=median, aes(ymax = ..y.., ymin = ..y..),color="black",geom = "errorbar" ,width = 0.19,size=0.6)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0.1,1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("Avargea beta-value of\nCpG sites in HAL1pro") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))
NAN_plot <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))

the <- compare_means(methy ~ His,  data = myTable1,paired = F ,method = "wilcox.test")
the 
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 2,
  y.position = c( 0.88,0.92,0.96)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="figs/beta_HAL1.pdf", plot=figure_2,bg = 'white', width =7, height = 8, units = 'cm', dpi = 600)

tmp <- myTable1[,c("case","methy")]
colnames(tmp) <- c("case","methy_HAL1pro")
write.table(tmp,file = "methy_HAL1pro.txt",row.names = F,quote = F,sep = "\t")

###############3
################# All L1
L1HS <- L1[which(L1$V5=="L1HS"),]
table(L1HS$V3)

myTable <- counts1[which(rownames(counts1) %in% L1HS$V1),]

L1HSmethy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
colnames(L1HSmethy) <- c("AID","methy")


myTable1 <- merge(L1HSmethy,samples)
#myTable1 <- myTable1[which(myTable1$Histology %in% c("YST","IMT","Teratoma","Germinoma")),]
myTable1$His <- myTable1$Histology
myTable1$His[which(myTable1$Histology!="YST") ] <- "aOther"
myTable1$His[which(myTable1$Histology=="Germinoma") ] <- "Germinoma"


NAN_plot <- ggplot(data=myTable1,aes(x=His,y=methy)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable1,aes(fill=Histology),shape=21,width = 0.3,size=2,alpha=0.75,stroke=0.1, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=median, aes(ymax = ..y.., ymin = ..y..),color="black",geom = "errorbar" ,width = 0.19,size=0.6)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0.2,1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("Avargea beta-value of\nCpG sites in L1HSpro") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))
NAN_plot <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))

the <- compare_means(methy ~ His,  data = myTable1,paired = F ,method = "wilcox.test")
the 
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 2,
  y.position = c( 0.88,0.92,0.96)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="figs/beta_L1HS.pdf", plot=figure_2,bg = 'white', width =7, height = 8, units = 'cm', dpi = 600)


tmp <- myTable1[,c("case","methy")]
colnames(tmp) <- c("case","methy_L1HSpro")
write.table(tmp,file = "methy_L1HSpro.txt",row.names = F,quote = F,sep = "\t")











##################33
TE <- read.table("../../PRDM2/GSEA/TEexpression/data/All_TE_list.txt",header = F)
colnames(TE) <- c("name","family","type")
TE <- TE[which(TE$family=="L1"),]

TEcounts <- read.table("All_TE_expression_TGCT.txt",header = T,row.names = "Geneid")
Info <- samples
rownames(Info) <- gsub("-","_",Info$case)

keep <- rowSums(TEcounts >= 20) >= 5
TEcounts1 <- TEcounts[keep,]
nrow(TEcounts1 )
nrow(TEcounts)

#TEcounts <- TEcounts1[which(rownames(TEcounts1) %in% TE$name),rownames(Info)] 
TEcounts <- TEcounts1[,rownames(Info)] 
dds=DESeqDataSetFromMatrix(countData = TEcounts,colData = Info,design = ~ Histology)

vsd <- vst(dds ,blind = TRUE,nsub=100)
head(assay(vsd), 3)
expMatrix_TGCT <- data.frame(assay(vsd))

rv <- rowVars(expMatrix_TGCT)
ntop <- 1000000
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,  length(rv)))]
pca <- prcomp(t(expMatrix_TGCT[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup <- c("AID","Histology","Dominant.His","Age")
if (!all(intgroup %in% colnames(Info))) {
  stop("the argument 'intgroup' should specify columns of colData(dds)")
}
intgroup.df <- as.data.frame(Info[, intgroup, drop = FALSE])
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
ggsave(file="./figs/PCA_Exp.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)
########




tmp <- data.frame(colnames(expMatrix_TGCT),t(expMatrix_TGCT["L1HS",]))
colnames(tmp) <- c("case","exp_L1HS")
write.table(tmp,file = "exp_L1HS.txt",row.names = F,quote = F,sep = "\t")

tmp <- data.frame(colnames(expMatrix_TGCT),t(expMatrix_TGCT["L1PA2",]))
colnames(tmp) <- c("case","exp_L1PA2")
write.table(tmp,file = "exp_L1PA2.txt",row.names = F,quote = F,sep = "\t")

tmp <- data.frame(colnames(expMatrix_TGCT),t(expMatrix_TGCT["HAL1",]))
colnames(tmp) <- c("case","exp_HAL1")
write.table(tmp,file = "exp_HAL1.txt",row.names = F,quote = F,sep = "\t")



L1exp <- data.frame(colnames(expMatrix_TGCT),t(expMatrix_TGCT["L1HS",]))
colnames(L1exp) <- c("case","exp")

L1HS <- L1[which(L1$V5=="L1HS"),]
myTable <- counts1[which(rownames(counts1) %in% L1HS$V1),]
L1HSmethy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
colnames(L1HSmethy) <- c("AID","methy")
myTable1 <- merge(L1HSmethy,samples)
myTable1$case <- gsub("-","_",myTable1$case)
qq <- merge(myTable1,L1exp)


cor.test(formula = ~ exp + methy, data = qq,method = "pearson")
cor.test(formula = ~ exp + methy, data = qq,method = "spearman",exact= F,continuity=T,conf.level = 0.95)

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = methy, y = exp),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = qq,  aes(x = methy, y = exp,fill = Histology),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1HS expression")+xlab("Average beta-value of L1HSpro")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(0.48,0.92),breaks = seq(-1,1,0.1))+scale_y_continuous(expand=c(0,0),limits = c(9.8,15.2),breaks=seq(0,50,1))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="figs/linear_L1HS_methy_vs_exp.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ exp + methy, data = qq,method = "spearman")

############3
xx <- qq[which(qq$Histology=="YST"),]
cor.test(formula = ~ exp + methy, data = xx,method = "spearman")


###############################
qq <- qq[which(qq$RT!="NA"),]

cor.test(formula = ~ exp + RT, data = qq,method = "spearman",exact= F,continuity=T,conf.level = 0.95)
cor.test(formula = ~ methy + RT, data = qq,method = "spearman",exact= F,continuity=T,conf.level = 0.95)

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = methy, y = RT),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = qq,  aes(x = methy, y = RT,fill = Histology),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1HS RT")+xlab("Average beta-value of L1HSpro")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(0.48,0.92),breaks = seq(-1,1,0.1))+scale_y_continuous(expand=c(0,0),limits = c(0,24),breaks=seq(0,50,5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="figs/linear_L1HS_methy_vs_RT.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ RT + methy, data = qq,method = "spearman")


##################################33
#qq$Histology <- qq$RNA_cluster
new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = exp, y = RT),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = qq,  aes(x = exp, y = RT,fill = Histology,shape=chr1p36),alpha=0.75,size=2.2)  
new.plot <- new.plot + geom_label_repel(data = qq,  aes(x = exp, y = RT,label = AID),alpha=0.9,size=1)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1HS RT")+xlab("RNA expresssion of L1HS")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(10,15),breaks=seq(0,50,1))+scale_y_continuous(expand=c(0,0),limits = c(-1,18),breaks=seq(0,50,5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))+ scale_shape_manual(values = c(21,22))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="figs/linear_L1HS_exp_vs_RT_all.pdf", plot=plot7,bg = 'white', width =8, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ RT + exp, data = qq,method = "spearman")

############3
xx <- qq[which(qq$RNA_cluster=="YST"),]
cor.test(formula = ~ RT + exp, data = xx,method = "spearman")
###################333

xx <- qq[which(qq$chr1p36=="Loss"),]
cor.test(formula = ~ RT + exp, data = xx,method = "spearman")
cor.test(formula = ~ RT + methy, data = xx,method = "spearman")


new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = xx, aes(x = exp, y = RT),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = xx,  aes(x = exp, y = RT,fill = Histology),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1HS RT")+xlab("RNA expresssion of L1HS")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(10.8,14.8),breaks=seq(0,50,1))+scale_y_continuous(expand=c(0,0),limits = c(-1,18),breaks=seq(0,50,5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="figs/linear_L1HS_exp_vs_RT_1p36loss.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ RT + exp, data = xx,method = "spearman")


xx <- qq[which(qq$chr1p36!="Loss"),]
cor.test(formula = ~ RT + exp, data = xx,method = "spearman")
cor.test(formula = ~ RT + methy, data = xx,method = "spearman")


new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = xx, aes(x = exp, y = RT),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = xx,  aes(x = exp, y = RT,fill = Histology),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1HS RT")+xlab("RNA expresssion of L1HS")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(10.8,14.8),breaks=seq(0,50,1))+scale_y_continuous(expand=c(0,0),limits = c(-1,18),breaks=seq(0,50,5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="figs/linear_L1HS_exp_vs_RT_1p36noloss.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ RT + exp, data = xx,method = "spearman")








new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = xx, aes(x = methy, y = RT),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = xx,  aes(x = methy, y = RT,fill = Histology),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1HS RT")+xlab("methy of L1HS")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(0.58,0.92),breaks = seq(-1,1,0.1))+scale_y_continuous(expand=c(0,0),limits = c(-1,18),breaks=seq(0,50,5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="figs/linear_L1HS_methy_vs_RT_1p36loss.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ RT + methy, data = xx,method = "spearman")


###################333
xx <- qq[which(qq$chr1p36!="Loss"),]
cor.test(formula = ~ RT + exp, data = xx,method = "spearman")
cor.test(formula = ~ RT + methy, data = xx,method = "spearman")
###################333




qq <- merge(myTable1,L1exp)
table1 <- qq
table1$His <- "YST"
table1$His[which(table1$Histology != "YST")] <- "aOther"
table1$His[which(table1$Histology == "Germinoma")] <- "Germinoma"


NAN_plot <- ggplot(data=table1,aes(x=His,y=exp)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=Histology),shape=21,width = 0.3,size=2,alpha=0.75,stroke=0.1, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=median, aes(ymax = ..y.., ymin = ..y..),color="black",geom = "errorbar" ,width = 0.19,size=0.6)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(10,16),breaks = seq(0,100,1)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("mRNA expression of L1HS") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))
NAN_plot <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))

the <- compare_means(exp ~ His,  data = table1,paired = F ,method = "wilcox.test")
the 
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 2,
  y.position = c( 14.8,15.2,15.6)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="figs/exp_L1HS.pdf", plot=figure_2,bg = 'white', width =7, height = 8, units = 'cm', dpi = 600)



L1exp <- data.frame(colnames(expMatrix_TGCT),t(expMatrix_TGCT["L1HS",]))
colnames(L1exp) <- c("case","exp")
myTable <- merge(myTable1,L1exp)
myTable$His <- "YST"
myTable$His[which(table1$Histology != "YST")] <- "aOther"
myTable$His[which(table1$Histology == "Germinoma")] <- "Germinoma"
myTable$rank <- 1
myTable$rank[which(myTable$His=="YST")] <- 3
myTable$rank[which(myTable$His=="Germinoma")] <- 2
library(ggpubr)
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(data=myTable ,aes(y=exp ,x=reorder(His,rank)),width=0.8,size=0.5,alpha=0,outlier.shape = NA)
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable ,aes(y=exp ,x=reorder(His,rank),fill=Histology),shape=21,width = 0.3,size=1.6,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0)
#geom_beeswarm(dodge.width=0.8,cex=2,data=table5,aes(y=Freq ,x=Type,fill=chr1p36),size=2.6,shape=21,alpha=0.7,stroke=0.5, stat = "identity")
NAN_plot <- NAN_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                            axis.text.x=element_text(size=11,angle=45,vjust = 1,hjust = 1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot+ggtitle(NULL)+ylab("L1HS expression")+xlab(NULL)
NAN_plot <- NAN_plot +  scale_y_continuous(expand=c(0,0),limits=c(9.6,15.4),breaks = seq(0,100,1))
NAN_plot <- NAN_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA097",YST="#FF99FF"))
plotxxx<-cbind(ggplotGrob(NAN_plot ),size="first")
ggsave(file=paste0("figs/exp_L1HS_final.pdf"), plot=plotxxx,bg = 'white', width = 6, height = 9, units = 'cm', dpi = 600)
the <- compare_means(exp ~ His,  data = myTable,method =  "wilcox")
the 





############

L1exp <- data.frame(colnames(expMatrix_TGCT),t(expMatrix_TGCT["L1PA2",]))
colnames(L1exp) <- c("case","exp")

L1PA2 <- L1[which(L1$V5=="L1PA2"),]
myTable <- counts1[which(rownames(counts1) %in% L1PA2$V1),]
L1PA2methy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
colnames(L1PA2methy) <- c("AID","methy")
myTable1 <- merge(L1PA2methy,samples)
myTable1$case <- gsub("-","_",myTable1$case)
qq <- merge(myTable1,L1exp)


cor.test(formula = ~ exp + methy, data = qq,method = "pearson")
cor.test(formula = ~ exp + methy, data = qq,method = "spearman",exact= F,continuity=T,conf.level = 0.95)

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = methy, y = exp),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = qq,  aes(x = methy, y = exp,fill = Histology),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1PA2,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1PA2 expression")+xlab("Average beta-value of L1PA2pro")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(0.48,0.92),breaks = seq(-1,1,0.1))+scale_y_continuous(expand=c(0,0),limits = c(11.8,17.2),breaks=seq(0,50,1))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="figs/linear_L1PA2_methy_vs_exp.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ exp + methy, data = qq,method = "spearman")


############3
xx <- qq[which(qq$Histology=="YST"),]
cor.test(formula = ~ exp + methy, data = xx,method = "spearman")


#table1 <- myTable[which(myTable$Histology == "YST" & myTable$chr1p36!="Not.Available"),]
table1 <- qq
table1$His <- "YST"
table1$His[which(table1$Histology != "YST")] <- "aOther"
table1$His[which(table1$Histology == "Germinoma")] <- "Germinoma"


#potential_ysts <- c("TCGA_YU_A94I","TCGA_2G_AAGT","TCGA_2G_AAFV","TCGA_YU_AA4L","TCGA_YU_AA61","TCGA_W4_A7U3","TCGA_2G_AAH4")
#table1$His[which(table1$case %in% potential_ysts)] <- "YST"

NAN_plot <- ggplot(data=table1,aes(x=His,y=exp)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=Histology),shape=21,width = 0.3,size=2,alpha=0.75,stroke=0.1, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=median, aes(ymax = ..y.., ymin = ..y..),color="black",geom = "errorbar" ,width = 0.19,size=0.6)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(12,18),breaks = seq(0,100,1)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("mRNA expression of L1PA2") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))
NAN_plot <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))

the <- compare_means(exp ~ His,  data = table1,paired = F ,method = "wilcox.test")
the 
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 2,
  y.position = c( 16.8,17.2,17.6)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="figs/exp_L1PA2.pdf", plot=figure_2,bg = 'white', width =7, height = 8, units = 'cm', dpi = 600)




##############################
L1exp <- data.frame(colnames(expMatrix_TGCT),t(expMatrix_TGCT["HAL1",]))
colnames(L1exp) <- c("case","exp")

HAL1 <- L1[which(L1$V5=="HAL1"),]
myTable <- counts1[which(rownames(counts1) %in% HAL1$V1),]
HAL1methy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
colnames(HAL1methy) <- c("AID","methy")
myTable1 <- merge(HAL1methy,samples)
myTable1$case <- gsub("-","_",myTable1$case)
qq <- merge(myTable1,L1exp)


cor.test(formula = ~ exp + methy, data = qq,method = "pearson")
cor.test(formula = ~ exp + methy, data = qq,method = "spearman",exact= F,continuity=T,conf.level = 0.95)

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = methy, y = exp),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = qq,  aes(x = methy, y = exp,fill = Histology),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = HAL1,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("HAL1 expression")+xlab("Average beta-value of HAL1pro")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(0.12,0.92),breaks = seq(-1,1,0.1))+scale_y_continuous(expand=c(0,0),limits = c(12.1,14.9),breaks=seq(0,50,0.5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="figs/linear_HAL1_methy_vs_exp.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ exp + methy, data = qq,method = "spearman")


############3
xx <- qq[which(qq$Final_His=="YST"),]
cor.test(formula = ~ exp + methy, data = xx,method = "spearman")


#table1 <- myTable[which(myTable$Histology == "YST" & myTable$chr1p36!="Not.Available"),]
table1 <- qq
table1$His <- "YST"
table1$His[which(table1$Histology != "YST")] <- "aOther"
table1$His[which(table1$Histology == "Germinoma")] <- "Germinoma"


#potential_ysts <- c("TCGA_YU_A94I","TCGA_2G_AAGT","TCGA_2G_AAFV","TCGA_YU_AA4L","TCGA_YU_AA61","TCGA_W4_A7U3","TCGA_2G_AAH4")
#table1$His[which(table1$case %in% potential_ysts)] <- "YST"

NAN_plot <- ggplot(data=table1,aes(x=His,y=exp)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=Histology),shape=21,width = 0.3,size=2,alpha=0.75,stroke=0.1, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=median, aes(ymax = ..y.., ymin = ..y..),color="black",geom = "errorbar" ,width = 0.19,size=0.6)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(12,16),breaks = seq(0,100,1)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("mRNA expression of HAL1") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))
NAN_plot <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST=gg_color_hue(5)[5],EC="grey40"))

the <- compare_means(exp ~ His,  data = table1,paired = F ,method = "wilcox.test")
the 
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 2,
  y.position = c( 14.8,15.2,15.6)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="figs/exp_HAL1.pdf", plot=figure_2,bg = 'white', width =7, height = 8, units = 'cm', dpi = 600)







myTable <- myTable1

L1_list <- unique(TE$name[which(TE$name %in%  rownames(expMatrix_TGCT))])

P_YST_vs_G <- rep(NA,length(L1_list))
Diff_YST_vs_G <- rep(NA,length(L1_list))
P_YST_vs_O <- rep(NA,length(L1_list))
Diff_YST_vs_O <- rep(NA,length(L1_list))
P_YST_vs_A <- rep(NA,length(L1_list))
Diff_YST_vs_A <- rep(NA,length(L1_list))
for(i in c(1:length(L1_list))){
  L1exp <- data.frame(colnames(expMatrix_TGCT),t(expMatrix_TGCT[L1_list[i],]))
  colnames(L1exp) <- c("case","exp")
  myTable1 <- merge(myTable,L1exp)
  myTable1$His <- myTable1$Histology
  myTable1$His[which(myTable1$Histology!="YST") ] <- "aOther"
  myTable1$His[which(myTable1$Histology=="Germinoma") ] <- "Germinoma"
  G1 <- myTable1$exp[which(myTable1$His=="aOther")]
  G2 <- myTable1$exp[which(myTable1$His=="Germinoma")]
  G3 <- myTable1$exp[which(myTable1$His=="YST")]
  G0 <- myTable1$exp[which(myTable1$His!="YST")]
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

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-2,2),breaks=seq(-80,80,0.5))+ scale_y_continuous(expand=c(0,0),limits=c(-0.01,2.8),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('Difference of average methylation beta value\n[YST - Germinoma]')+ylab('-log10(FDR)')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_Exp_YST_vs_G_fdr.pdf", plot=AllTumorure_2,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)

##############3
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_point(data = resultsCNV , aes(x = Diff_YST_vs_A, y = -log10(FDR_YST_vs_A)),fill="grey90",color="grey10",alpha=0.7,size=2, shape = 21,stroke=0.5) 
NAN_plot <- NAN_plot + geom_hline(aes(yintercept= -log10(0.05) ),color='black',size=.5,linetype='dashed') 
NAN_plot <- NAN_plot + geom_vline(aes(xintercept=0),color='black',size=.5,linetype=2)
NAN_plot <- NAN_plot + geom_point(data = name, aes(x = Diff_YST_vs_G, y = -log10(FDR_YST_vs_A)),fill="#c51b8a",size=2,shape=21,alpha=0.7,stroke=0.5)#,fill="#009BFF") 

NAN_plot <- NAN_plot + geom_text_repel(data = resultsCNV  , aes(x = Diff_YST_vs_G, y = -log10(FDR_YST_vs_A),label=L1_list),
                                       segment.color = "black",size =2,color ='black',min.segment.length = 0.4 , box.padding = 0.2)

NAN_plot <- NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,2,1,2),'lines'),
                             plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                             legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=12,face='plain'),legend.title=element_text(size=12,face='plain'),
                             axis.text.x=element_text(size=14,face='plain',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),
                             axis.title.x=element_text(size=16,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-2,2),breaks=seq(-80,80,0.5))+ scale_y_continuous(expand=c(0,0),limits=c(-0.01,2.8),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('Difference of average methylation beta value\n[YST - All]')+ylab('-log10(FDR)')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_Exp_YST_vs_A_fdr.pdf", plot=AllTumorure_2,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)



