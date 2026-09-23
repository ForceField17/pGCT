# pGCT
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
library(oncoprint)
library(reshape2)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



#########################
table0 <- read.table("./results/Final_PointMut.txt",sep="\t",  head=T,row.names = "AID")
table0 <- table0[,which(colnames(table0) != "X7p" & colnames(table0) != "X7q"  & colnames(table0) != "X20p" & colnames(table0) != "LRP2" & colnames(table0)!="BRCA1")]
temp <- melt(as.matrix(table0 ))
#temp$value[which(temp$Var1 %in% c("C002","F026","F028","F049","F052","SJ02","SJ06","SJ08","SJ09") & temp$Var2=="KIT")] <- "amplification"
temp <- temp[which(temp$value!="WT"),]
temp$His <- "SomaticMut"
temp$His[which(temp$Var1=="SJ18" & temp$Var2=="PTEN")] <- "GermlineMut_and_SomaticLOH"
temp$His[which(temp$Var1=="C002" & temp$Var2=="ATM")] <- "GermlineMut_and_SomaticLOH"
temp$His[which(temp$Var1=="SJ17" & temp$Var2=="CBL")] <- "GermlineMut_and_SomaticLOH"
temp$His[which(temp$Var1=="SJ20" & temp$Var2=="CHEK2")] <- "GermlineMut_and_SomaticLOH"
temp$His[which(temp$Var1=="F049" & temp$Var2=="BLM")] <- "Unknown"

table1 <- temp
temp <- table1[which(table1$His=="SomaticMut"),]
myTable <- as.data.frame(table(temp$Var2))
myTable$His <- "SomaticMut"

temp <- table1[which(table1$His=="GermlineMut_and_SomaticLOH"),]
table2 <- as.data.frame(table(temp$Var2))
table2$His <- "GermlineMut_and_SomaticLOH"
myTable <- rbind(myTable,table2)

temp <- table1[which(table1$His=="Unknown"),]
table2 <- as.data.frame(table(temp$Var2))
table2$His <- "Unknown"
myTable <- rbind(myTable,table2)

temp <- melt(as.matrix(table0 ))
#temp$value[which(temp$Var1 %in% c("C002","F026","F028","F049","F052","SJ02","SJ06","SJ08","SJ09") & temp$Var2=="KIT")] <- "amplification"
temp <- temp[which(temp$value!="WT"),]
temp <- as.data.frame(table(temp$Var2));colnames(temp) <- c("Var1","Count"); temp$Rank <- c(1:nrow(temp))
myTable <- merge(myTable,temp,1,1)

myTable$PcentCount <- myTable$Count/nrow(table0)*100
myTable$PcentFreq <- myTable$Freq/nrow(table0)*100
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot +
  geom_bar(data=myTable,aes(x=reorder(Var1,Rank),y=PcentFreq,fill=His),alpha=1,color="grey20",size=0.4,width=0.6,stat='identity',position=position_stack()) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='top',legend.text=element_text(size=10,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,angle=90,face='italic',vjust=0.5,hjust=1,color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='italic',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_discrete() + scale_y_continuous(expand=c(0,0),limits=c(0,16),breaks = seq(0,100,4) )
NAN_plot<- NAN_plot +ylab("Percentage") +xlab(NULL)#+coord_flip()
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values=c("#fdb863","#998ec3","white"))
#NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot3 <- NAN_plot
length1 <- length(unique(myTable$Var1))

#####################################33
####
Sample_features <- read.table("../../lib/Clinical.txt",sep = "\t",header = T)
WGD <- read.table('./results/Final_WGD.txt',header = T)

tmp <- merge(Sample_features,WGD,1,1);rownames(tmp) <- tmp$AID; tmp <- tmp[,-1]
tmp <- tmp[,c("Seq", "Histology", "Domin.His", "Subtype", "Age", "Site", "Gender", "Karyotype", "WGD")]
Clinical <- melt(as.matrix(tmp))
colnames(Clinical) <- c('AID','gene','Info')
#####
Arm <- read.table('./results/Final_Arm.txt',header = T,row.names = "AID")
Focal <- read.table('./results/Final_Focal.txt',header = T,row.names = "AID")
CNA <- cbind(Arm,Focal); CNA <- CNA[,c("chrXq","chr22q11","chr21q","chr20q","chr18q21","chr16p","chr13q","chr12p","chr11q24","chr8q13","chr8p","chr6q25","chr4q","chr4q12","chr3p23","chr2p24","chr1q21","chr1p36")]

CNV <- melt(as.matrix(CNA))
colnames(CNV) <- c('AID','gene','CNV')
CNV <- CNV[which(CNV$CNV=="Gain"),]
CNV$His <- "SomaticMut"
CNV$His[which(CNV$AID %in% c("NG5","NG7","SJ14") & CNV$gene=="chr21q")] <- "GermlineMut"
CNV$His[which(CNV$AID %in% c("NG8","NG12") & CNV$gene=="chrXq")] <- "GermlineMut"
CNV$His[which(CNV$AID %in% c("C010","C010r") & CNV$gene=="chrXq")] <- "GermlineMut"

table1 <- CNV
temp <- table1[which(table1$His=="SomaticMut"),]
myTable <- as.data.frame(table(temp$gene))
myTable$His <- "SomaticMut"

temp <- table1[which(table1$His!="SomaticMut"),]
table2 <- as.data.frame(table(temp$gene))
table2$His <- "GermlineMut"
myTable <- rbind(myTable,table2)

temp <- melt(as.matrix(CNA ))
temp <- temp[which(temp$value=="Gain"),]
temp <- as.data.frame(table(temp$Var2));colnames(temp) <- c("Var1","Count"); temp$Rank <- c(1:nrow(temp))
myTable <- merge(myTable,temp,1,1)

myTable$PcentCount <- myTable$Count/nrow(CNA)*100
myTable$PcentFreq <- myTable$Freq/nrow(CNA)*100

#######
CNV <- melt(as.matrix(CNA))
colnames(CNV) <- c('AID','gene','CNV')
CNV <- CNV[which(CNV$CNV=="Loss"),]
CNV$His <- "SomaticMut"
CNV$His[which(CNV$AID %in% c("C010","C010r") & CNV$gene=="chrXq")] <- "GermlineMut"

table1 <- CNV
temp <- table1[which(table1$His=="SomaticMut"),]
myTable2 <- as.data.frame(table(temp$gene))
myTable2$His <- "SomaticMut"

temp <- table1[which(table1$His!="SomaticMut"),]
table2 <- as.data.frame(table(temp$gene))
table2$His <- "GermlineMut"
myTable2 <- rbind(myTable2,table2)

temp <- melt(as.matrix(CNA ))
temp <- temp[which(temp$value=="Loss"),]
temp <- as.data.frame(table(temp$Var2));colnames(temp) <- c("Var1","Count"); temp$Rank <- c(1:nrow(temp))
myTable2 <- merge(myTable2,temp,1,1)

myTable2$PcentCount <- -myTable2$Count/nrow(CNA)*100
myTable2$PcentFreq <- -myTable2$Freq/nrow(CNA)*100
#######
myTable <- rbind(myTable,myTable2)
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_bar(data=myTable,aes(x=reorder(Var1,-Rank),y=PcentFreq,fill=His),alpha=1,color="grey20",size=0.4,width=0.6,stat='identity',position=position_stack()) +
  geom_hline(yintercept = 0,size=0.5,color="black")
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                            plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='top',legend.text=element_text(size=10,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,angle=90,face='plain',vjust=0.5,hjust=1,color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_discrete() + scale_y_continuous(expand=c(0,0),limits=c(-35,35),breaks = seq(-60,60,15) )
NAN_plot<- NAN_plot +ylab("Percentage") +xlab(NULL)#+coord_flip()
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values=c("#e08214","#998ec3"))
#NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot4 <- NAN_plot

length2 <- length(unique(myTable$Var1))


figure<-cbind(ggplotGrob(NAN_plot3),ggplotGrob(NAN_plot4),size="last")
panels <- figure$layout$t[grep("panel", figure$layout$name)]

figure$widths[7]  <- unit(length1/32,'null')
figure$widths[20] <- unit(length2/32,'null')


ggsave(file="Summary_percentage_SG.pdf", plot=figure,bg = 'white', width =28, height = 8, units = 'cm', dpi = 600)

