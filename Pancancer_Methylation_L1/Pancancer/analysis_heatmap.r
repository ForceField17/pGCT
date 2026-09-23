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
  hcl(h = hues, l = 65, c = 100)[1:n]
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

TCGA <- c("TCGA-ACC" ,"TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-DLBC",
          "TCGA-ESCA","TCGA-GBM" ,"TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LAML",
          "TCGA-LGG" ,"TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-MESO","TCGA-OV","TCGA-PAAD",
          "TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-THCA",
          "TCGA-THYM","TCGA-UCEC","TCGA-UCS" ,"TCGA-UVM" ,"TCGA-TGCT")

L1 <- read.table("../450k_L1pro_list.bed",header = F)
L1 <- L1[,c("V1","V5")]
colnames(L1) <- c("AID","type")


### All L1pro
All_L1methy <- c()
for(i in c(1:length(TCGA))){
  myTable <- read.table(paste0("./data/matrix_",TCGA[i],".txt"),header = T, row.names = "A_CpG_ID")
  L1methy <- as.data.frame(rowMeans((myTable),na.rm = T));L1methy$AID <- rownames(L1methy)
  colnames(L1methy) <- c("methy","AID")
  L1methy <- L1methy[which(L1methy$methy!="NaN"),]
  L1methy <- merge(L1methy,L1)
  
  tmpTable <- data.frame(aggregate(L1methy$methy, list(L1methy$type), FUN=mean))
  colnames(tmpTable) <- c("LINE1","mean_methy")
  tmpTable$Cancer <- TCGA[i]
  All_L1methy <- rbind(All_L1methy,tmpTable)
}


###################
samples <- read.table("../../Transcriptome_subtyping/RNA_TCGA/Final_Genetic_profiles_2025.txt",sep = "\t",row.names = "SampleID",header = T)
samples$AID <- gsub("-",".",rownames(samples))
Germinoma <- samples$AID[which(samples$Histology=="Germinoma")]
YST <- samples$AID[which(samples$Histology=="YST")]
OtherGCT <- samples$AID[which(samples$Histology!="YST" & samples$Histology!="Germinoma" )]
#################3
myTable <- read.table("../TGCT/matrix_TCGA-TGCT.txt",header = T, row.names = "A_CpG_ID")

myTable1 <- myTable[,YST]
L1methy <- as.data.frame(rowMeans((myTable1),na.rm = T));L1methy$AID <- rownames(L1methy)
colnames(L1methy) <- c("methy","AID")
L1methy <- L1methy[which(L1methy$methy!="NaN"),]
L1methy <- merge(L1methy,L1)
tmpTable <- data.frame(aggregate(L1methy$methy, list(L1methy$type), FUN=mean))
colnames(tmpTable) <- c("LINE1","mean_methy")
tmpTable$Cancer <- "YST"
All_L1methy <- rbind(All_L1methy,tmpTable)

myTable1 <- myTable[,Germinoma]
L1methy <- as.data.frame(rowMeans((myTable1),na.rm = T));L1methy$AID <- rownames(L1methy)
colnames(L1methy) <- c("methy","AID")
L1methy <- L1methy[which(L1methy$methy!="NaN"),]
L1methy <- merge(L1methy,L1)
tmpTable <- data.frame(aggregate(L1methy$methy, list(L1methy$type), FUN=mean))
colnames(tmpTable) <- c("LINE1","mean_methy")
tmpTable$Cancer <- "Germinoma"
All_L1methy <- rbind(All_L1methy,tmpTable)

myTable1 <- myTable[,OtherGCT]
L1methy <- as.data.frame(rowMeans((myTable1),na.rm = T));L1methy$AID <- rownames(L1methy)
colnames(L1methy) <- c("methy","AID")
L1methy <- L1methy[which(L1methy$methy!="NaN"),]
L1methy <- merge(L1methy,L1)
tmpTable <- data.frame(aggregate(L1methy$methy, list(L1methy$type), FUN=mean))
colnames(tmpTable) <- c("LINE1","mean_methy")
tmpTable$Cancer <- "OtherGCT"
All_L1methy <- rbind(All_L1methy,tmpTable)


myTable1 <- All_L1methy
myRank <- data.frame(aggregate(myTable1$mean_methy, list(myTable1$Cancer), FUN=mean))
colnames(myRank) <- c("Cancer","rank")
new_table1 <- merge(myTable1,myRank)
Age <- read.table("../L1_age.txt",header = T)
new_table1 <- merge(new_table1,Age)

F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table1,aes(x=reorder(Cancer,rank),y=reorder(LINE1,L1_age) ,fill=mean_methy),color='black',width=1,height=1,size=0.25,stat='identity')
F1A.pLot<-F1A.pLot+ scale_fill_gradientn(name=NULL,colours=c('#29419E','#698DC9',"grey90",'#F67172','#DC2B18'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=9,face='plain'),legend.key.width=unit(0.2,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='right',axis.ticks = element_blank(),
                         legend.direction='vertical',legend.text=element_text(size=9,face='plain'),axis.line = element_blank(),
                         axis.text.x =element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=9,face='italic',color='black'),
                         axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1A.pLot

figure<-cbind(ggplotGrob(F1A.pLot),size="last")
ggsave(file="Heatmap_test.pdf", plot=figure,bg = 'white', width =18, height = 24, units = 'cm', dpi = 600)

#########33
library(reshape2)
library(tidyr)

# 转换为矩阵
matrix_result <- acast(All_L1methy, LINE1 ~ Cancer, value.var = "mean_methy")
cohort <- myRank$Cancer[order(myRank$rank,decreasing = T)]
Age <- Age[which(Age$LINE1 %in% rownames(matrix_result)),]
matrix_result <- matrix_result[Age$LINE1[order(Age$L1_age,decreasing = T)],c(cohort[grep("TCGA",cohort)],"OtherGCT","YST","Germinoma")]


pdf(file = paste0("./heatmap_all.pdf"),width=12,height=6)
pheatmap(t(matrix_result),cutree_cols= 4,cutree_rows= 2,cluster_cols = F,cluster_rows = F,treeheight_col = 30,
         color = colorRampPalette(c("#000099","#0000CC","#0000FF","#3333FF","#6666FF","#9999FF","#CCCCFF", "grey95","#FFCCCC","#FF9999","#FF6666","#FF3333", "#FF0000","#CC0000","#990000"))(50),
         annotation_names_row=F, gaps_col = c(3,14,48,52,57), 
         annotation_names_col=F, gaps_row = c(33),
         show_rownames = T,annotation_legend = T,cellheight=9.5,
         scale = "column")
dev.off()



##################test

myTable1 <- All_L1methy
myRank <- data.frame(aggregate(myTable1$mean_methy, list(myTable1$Cancer), FUN=mean))
colnames(myRank) <- c("Cancer","rank")
new_table1 <- merge(myTable1,myRank)
Age <- read.table("../L1_age.txt",header = T)
new_table1 <- merge(new_table1,Age)

L1HS <- myTable1[which(myTable1$LINE1=="L1HS"),]
xxx <- merge(L1HS,myRank)
xxx$L1HS_frac <- xxx$mean_methy/xxx$rank


# 转换为矩阵
matrix_result <- acast(All_L1methy, LINE1 ~ Cancer, value.var = "mean_methy")
cohort <- myRank$Cancer[order(myRank$rank,decreasing = T)]
Age <- Age[which(Age$LINE1 %in% rownames(matrix_result)),]
matrix_result <- matrix_result[Age$LINE1[order(Age$L1_age,decreasing = T)],c(cohort[grep("TCGA",cohort)],"OtherGCT","YST","Germinoma")]


pdf(file = paste0("./heatmap_all.pdf"),width=12,height=6)
pheatmap(t(matrix_result),cutree_cols= 4,cutree_rows= 2,cluster_cols = F,cluster_rows = F,treeheight_col = 30,
         color = colorRampPalette(c("#000099","#0000CC","#0000FF","#3333FF","#6666FF","#9999FF","#CCCCFF", "grey95","#FFCCCC","#FF9999","#FF6666","#FF3333", "#FF0000","#CC0000","#990000"))(50),
         annotation_names_row=F, gaps_col = c(3,14,48,52,57), 
         annotation_names_col=F, gaps_row = c(33),
         show_rownames = T,annotation_legend = T,cellheight=9.5,
         scale = "column")
dev.off()



