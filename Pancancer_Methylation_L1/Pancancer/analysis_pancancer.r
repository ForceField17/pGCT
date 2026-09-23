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

TCGA <- c("TCGA-ACC" ,"TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-DLBC",
          "TCGA-ESCA","TCGA-GBM" ,"TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LAML",
          "TCGA-LGG" ,"TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-MESO","TCGA-OV","TCGA-PAAD",
          "TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-THCA",
          "TCGA-THYM","TCGA-UCEC","TCGA-UCS" ,"TCGA-UVM" ,"TCGA-TGCT")

L1 <- read.table("../450k_L1pro_list.bed",header = F)
L1HS <- L1[which(L1$V5=="L1HS"),]

#myTable <- counts1[which(rownames(counts1) %in% L1HS$V1),]
#L1HSmethy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
#colnames(L1HSmethy) <- c("AID","methy")


### All L1pro
All_L1methy <- c()
for(i in c(1:length(TCGA))){
  myTable <- read.table(paste0("./data/matrix_",TCGA[i],".txt"),header = T, row.names = "A_CpG_ID")
  L1methy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
  colnames(L1methy) <- c("AID","methy")
  L1methy$Cancer <- TCGA[i]
  All_L1methy <- rbind(All_L1methy,L1methy)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

##################
Sample_features <- read.table("../../TGCT_subtyping/Info_temp.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- gsub("-",".",samples$AID)
samples$AID <- rownames(samples) 
Germinoma <- samples$AID[which(samples$Histology=="Germinoma")]
YST <- samples$AID[which(samples$Histology=="YST")]

myTable1 <- All_L1methy
myTable1$His <- "others"
myTable1$His[which(myTable1$Cancer=="TCGA-TGCT")] <- "OtherGCT"
myTable1$His[which(myTable1$AID %in% myTable1[YST,]$AID)] <- "YST"
myTable1$His[which(myTable1$AID %in% myTable1[Germinoma,]$AID)] <- "Germinoma"
myTable1[YST,]

myRank <- data.frame(aggregate(myTable1$methy, list(myTable1$Cancer), FUN=mean))
colnames(myRank) <- c("Cancer","rank")
myTable1 <- merge(myTable1,myRank)

NAN_plot <- ggplot(data=myTable1,aes(x=reorder(Cancer,-rank) ,y=methy)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable1,aes(fill=His),shape=21,width = 0.3,size=1.2,alpha=0.75,stroke=0.01, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=mean, aes(ymax = ..y.., ymin = ..y..),color="blue",geom = "errorbar" ,width = 0.3,size=0.6)
#NAN_plot <- NAN_plot + geom_point(aes(x=reorder(Cancer,rank) ,y=rank),color="blue",shape=9)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0.2,1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("Avargea beta-value of\nCpG sites in L1pro") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF13FF",EC="grey40",OtherGCT="#a6dba0",others="#e7d4e8"))
NAN_plot1 <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF13FF",EC="grey40",OtherGCT="#a6dba0",others="#e7d4e8"))

the <- compare_means(methy ~ Cancer,  data = myTable1,paired = F ,method = "wilcox.test")
the 

#################
qq <- myTable1[which(myTable1$Cancer == "TCGA-TGCT"),]
qq$rank <- 1; qq$rank[which(qq$His=="Germinoma")] <- 3;qq$rank[which(qq$His=="YST")] <- 2
NAN_plot <- ggplot(data=qq,aes(x=reorder(His,rank) ,y=methy)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=qq,aes(fill=His),shape=21,width = 0.3,size=1.2,alpha=0.75,stroke=0.01, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=mean, aes(ymax = ..y.., ymin = ..y..),color="blue",geom = "errorbar" ,width = 0.3,size=0.6)
#NAN_plot <- NAN_plot + geom_point(aes(x=reorder(Cancer,rank) ,y=rank),color="blue",shape=9)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,0.2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_blank())
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0.2,1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("Avargea beta-value of\nCpG sites in L1pro") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF13FF",EC="grey40",OtherGCT="#a6dba0",others="#e7d4e8"))
NAN_plot2 <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF13FF",EC="grey40",OtherGCT="#a6dba0",others="#e7d4e8"))

the <- compare_means(methy ~ His,  data = qq,paired = F ,method = "wilcox.test")
the 

figure<-cbind(ggplotGrob(NAN_plot1),ggplotGrob(NAN_plot2),size="last")
panels <- figure$layout$t[grep("panel", figure$layout$name)]

figure$widths[7]  <- unit(33,'null')
figure$widths[20] <- unit(3,'null')

ggsave(file="beta_L1pro.pdf", plot=figure,bg = 'white', width =28, height = 8, units = 'cm', dpi = 600)






##############################################33L1HSpro
All_L1HSmethy <- c()
for(i in c(1:length(TCGA))){
  counts1 <- read.table(paste0("./data/matrix_",TCGA[i],".txt"),header = T, row.names = "A_CpG_ID")
  myTable <- counts1[which(rownames(counts1) %in% L1HS$V1),]
  L1methy <- data.frame(colnames(myTable),rowMeans(t(myTable),na.rm = T))
  colnames(L1methy) <- c("AID","methy")
  L1methy$Cancer <- TCGA[i]
  All_L1HSmethy <- rbind(All_L1HSmethy,L1methy)
}

##################

myTable1 <- All_L1HSmethy
myTable1$His <- "others"
myTable1$His[which(myTable1$Cancer=="TCGA-TGCT")] <- "OtherGCT"
myTable1$His[which(myTable1$AID %in% myTable1[YST,]$AID)] <- "YST"
myTable1$His[which(myTable1$AID %in% myTable1[Germinoma,]$AID)] <- "Germinoma"
myTable1[YST,]

myRank <- data.frame(aggregate(myTable1$methy, list(myTable1$Cancer), FUN=mean))
colnames(myRank) <- c("Cancer","rank")
myTable1 <- merge(myTable1,myRank)

NAN_plot <- ggplot(data=myTable1,aes(x=reorder(Cancer,-rank) ,y=methy)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable1,aes(fill=His),shape=21,width = 0.3,size=1.2,alpha=0.75,stroke=0.01, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=mean, aes(ymax = ..y.., ymin = ..y..),color="blue",geom = "errorbar" ,width = 0.3,size=0.6)
#NAN_plot <- NAN_plot + geom_point(aes(x=reorder(Cancer,rank) ,y=rank),color="blue",shape=9)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0.2,1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("Avargea beta-value of\nCpG sites in L1HSpro") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF13FF",EC="grey40",OtherGCT="#a6dba0",others="#e7d4e8"))
NAN_plot1 <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF13FF",EC="grey40",OtherGCT="#a6dba0",others="#e7d4e8"))

the <- compare_means(methy ~ Cancer,  data = myTable1,paired = F ,method = "wilcox.test")
the 

#################
qq <- myTable1[which(myTable1$Cancer == "TCGA-TGCT"),]
qq$rank <- 1; qq$rank[which(qq$His=="Germinoma")] <- 2;qq$rank[which(qq$His=="YST")] <- 3
NAN_plot <- ggplot(data=qq,aes(x=reorder(His,rank) ,y=methy)) + theme_classic() 
NAN_plot <- NAN_plot + geom_quasirandom(data=qq,aes(fill=His),shape=21,width = 0.3,size=1.2,alpha=0.75,stroke=0.01, varwidth = T,dodge.width = 0)
NAN_plot <- NAN_plot + stat_summary(fun=mean, aes(ymax = ..y.., ymin = ..y..),color="blue",geom = "errorbar" ,width = 0.3,size=0.6)
#NAN_plot <- NAN_plot + geom_point(aes(x=reorder(Cancer,rank) ,y=rank),color="blue",shape=9)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,0.2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),
                            axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_blank())
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0.2,1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot 
NAN_plot<- NAN_plot +ylab("Avargea beta-value of\nCpG sites in L1HSpro") +xlab(NULL)
NAN_plot <- NAN_plot +   scale_fill_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF13FF",EC="grey40",OtherGCT="#a6dba0",others="#e7d4e8"))
NAN_plot2 <- NAN_plot +   scale_color_manual(name=NULL,values=c(Germinoma=gg_color_hue(5)[1],IMT=gg_color_hue(5)[2],Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF13FF",EC="grey40",OtherGCT="#a6dba0",others="#e7d4e8"))

the <- compare_means(methy ~ His,  data = qq,paired = F ,method = "wilcox.test")
the 

figure<-cbind(ggplotGrob(NAN_plot1),ggplotGrob(NAN_plot2),size="last")
panels <- figure$layout$t[grep("panel", figure$layout$name)]

figure$widths[7]  <- unit(33,'null')
figure$widths[20] <- unit(3,'null')

ggsave(file="beta_L1HSpro.pdf", plot=figure,bg = 'white', width =28, height = 8, units = 'cm', dpi = 600)



