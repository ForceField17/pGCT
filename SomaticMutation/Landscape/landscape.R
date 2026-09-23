# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

####
Sample_features <- read.table("../../lib/Clinical.txt",sep = "\t",header = T)
WGD <- read.table('./results/Final_WGD.txt',header = T)

tmp <- merge(Sample_features,WGD,1,1);rownames(tmp) <- tmp$AID; tmp <- tmp[,-1]
#tmp <- tmp[,c("Seq", "Histology", "Domin.His", "Subtype", "Age", "Site", "Gender", "GInstability", "Karyotype", "WGD")]
tmp <- tmp[,c("Seq", "Histology", "Domin.His", "Subtype", "Age", "Site", "Gender", "Karyotype", "WGD")]
Clinical <- melt(as.matrix(tmp))
colnames(Clinical) <- c('AID','gene','Info')
#####

#####
PointMut<- read.table('./results/Final_PointMut.txt',header = T,row.names = "AID"); #PointMut<- PointMut[,rev(colnames(PointMut))]
PointMut<- PointMut[,c("CHEK2","ATM","BLM","ARID4B","ARID1A","ARID1B","BCORL1","BCOR","PIK3CD","PIK3CA","PIK3CB","AKT3","MTOR","PTEN","SOS1","NF1","USP28","FGFR2","CBL" ,"NRAS","RRAS2","KRAS","KIT")] 
Arm <- read.table('./results/Final_Arm.txt',header = T,row.names = "AID")
Focal <- read.table('./results/Final_Focal.txt',header = T,row.names = "AID")
CNA <- cbind(Arm,Focal); CNA <- CNA[,c("chrXq","chr22q11","chr21q","chr20q","chr18q21","chr16p","chr13q","chr12p","chr11q24","chr8q13","chr8p","chr6q25","chr4q12","chr4q","chr3p23","chr2p24","chr1q21","chr1p36")]
LOH <- read.table('./results/Final_LOH.txt',header = T,row.names = "AID")

Mut <- melt(as.matrix(PointMut))
colnames(Mut) <- c('AID','gene','Point')
LOH <- melt(as.matrix(LOH))
colnames(LOH) <- c('AID','gene','LOH')
CNV <- melt(as.matrix(CNA))
colnames(CNV) <- c('AID','gene','CNV')
KITamp <- CNV[which(CNV$gene=="chr4q12" & CNV$CNV=="Gain"),]
KITamp$gene <- "KIT"
#####

###Germinoma
#####################################################################################################################################################################
sample <- as.data.frame(Sample_features[which(Sample_features$Histology=="Germinoma"),])
sortList <- read.table("../../lib/sort_Germinoma.txt",header = T)
sample <- merge(sample,sortList,1,1); scaleW1 <- nrow(sample)

mut1 <- Mut[which(Mut$AID %in% sample$AID),]
mut2 <- LOH[which(LOH$AID %in% sample$AID),]
mut3 <- CNV[which(CNV$AID %in% sample$AID),]; mut3 <- mut3[which(mut3$gene %in% c("chrXq","chrXp","chr22q11","chr21q","chr20q","chr18q21","chr16p","chr13q","chr12p","chr11q24","chr9q","chr9p","chr8q","chr8q13","chr8p","chr4q","chr3p23","chr2p24","chr1q","chr1q21","chr6q25","chr1p36")),]
geneCNV <- KITamp[which(KITamp$AID %in% sample$AID),]

new_table1 <- merge(mut1,sample,1,1)
new_table2 <- merge(mut2,sample,1,1);new_table2 <- new_table2[which(new_table2$LOH=="LOH"),]
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table1,aes(x=reorder(AID,order),y=gene,fill=Point),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+geom_point(data = geneCNV,aes(x=AID,y=gene),shape=17,color='#CD0000',size=2.6)
F1A.pLot<-F1A.pLot+geom_text(data = new_table2,aes(x=AID,y=gene),label="L",color='black',size=3.5)
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',missense='#66c2a5',inframeIndel='#a6d854',frameshift='#fc8d62',splice="#619CFF",new_splice="#619CFF",nonsense='#df65b0',f_NA='grey'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part1 <- F1A.pLot
Part1
gene.scale1 = 1 + (ncol(PointMut) - 14)/14

new_table3 <- merge(mut3,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table3,aes(x=reorder(AID,order),y=gene,fill=CNV),color='white',width=1,height=1,size=0.4,stat='identity')

F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Gain='#CD0000',Loss='#1400CE',f_NA='grey'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part2 <- F1A.pLot
Part2
gene.scale2 = 1 + (ncol(CNA) -1 - 14)/14

##############
Info <- Clinical[which(Clinical$AID %in% sample$AID),]
data1 <- Info[which(Info$gene %in% c("Domin.His","Site","Gender","Race","Seq")),]
data2 <- Info[which(Info$gene %in% c("WGD","Karyotype","GInstability")),]
data3 <- Info[which(Info$gene %in% c("Age")),]

new_table4 <- merge(data1,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table4,aes(x=reorder(AID,order),y=gene,fill=Info),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Germinoma="#FFA49C",IMT="#CBCE00",Teratoma="#1aa3d9",YST="#FF9FFF",EmbryonalC="#BDBAD7",f_NA='grey',Mix='grey',NGGCT='grey',No='grey',NANA='grey95',
                                    Asian="#ffd92f",China="#ffd92f" ,Japan="#ffd92f", Black="#386cb0", White= "#80cdc1", Indian="#b3de69", White.Asian="#b3de69" , White.Black="#b3de69", White.Indian="#b3de69",WES="#93C6E7",WGS="#54278f",
                                    Intracranial="#9970ab",Extracranial="#a6dba0",Ovary="#a6dba0",Pelvis="#a6dba0",Testis="#a6dba0",Male='#b3cde3',Female='#fccde5',Normal='grey95',Trisomy21="#F6BA6F",XO="#6DA9E4",XXY='#ADE4DB'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part4 <- F1A.pLot
Part4
gene.scale4 = 1 + (length(unique(new_table4$gene)) - 14)/14

new_table5 <- merge(data2,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table5,aes(x=reorder(AID,order),y=gene,fill=Info),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Mut="#dd3497",Haploid="#3A8891",f_NA='grey',Normal='grey95',Trisomy21="#9970ab",XO="#00B4FF",XXY='#5aae61',NANA='grey95',XY="#b15928"))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part3 <- F1A.pLot
Part3
gene.scale3 = 1 + (length(unique(new_table5$gene)) - 14)/14*0.99

new_table6 <- merge(data3,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_hline(yintercept = 3,color="#f46d43",linetype=2)
F1A.pLot<-F1A.pLot+geom_point(data = new_table6,aes(x=reorder(AID,order),y=as.numeric(as.character(Info))),shape=16,color='black',size=2)
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='grey98',color='transparent',size=1),plot.margin=unit(c(0.45,0.1,1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),axis.line.x = element_blank(),
                         axis.text.x=element_text(size=10,angle = 90,vjust = 0.5,hjust = 1,face="plain",color="black"),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete()+xlab(NULL)+ylab(NULL)+scale_y_continuous(expand=c(0,0),limits=c(-1,21),breaks = seq(0,20,10) )
Part5 <- F1A.pLot
Part5
gene.scale5 = 1 + (3 - 14)/14

figure_1 <- rbind(ggplotGrob(Part1),ggplotGrob(Part2),ggplotGrob(Part3),ggplotGrob(Part4),ggplotGrob(Part5),size="first")
panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
figure_1$heights[panels][1] <- unit(gene.scale1,'null')
figure_1$heights[panels][2] <- unit(gene.scale2,'null')
figure_1$heights[panels][3] <- unit(gene.scale3,'null')
figure_1$heights[panels][4] <- unit(gene.scale4,'null')
figure_1$heights[panels][5] <- unit(gene.scale5,'null')
ggsave(file="test_Germinoma.pdf", plot=figure_1,bg = 'white', width = 20, height = 30, units = 'cm', dpi = 600)
#####################################################################################################################################################################


###YST
#####################################################################################################################################################################
sample <- as.data.frame(Sample_features[which(Sample_features$Histology=="YST"),])
sortList <- read.table("../../lib/sort_YST.txt",header = T)
sample <- merge(sample,sortList,1,1); scaleW2 <- nrow(sample)

mut1 <- Mut[which(Mut$AID %in% sample$AID),]
mut2 <- LOH[which(LOH$AID %in% sample$AID),]
mut3 <- CNV[which(CNV$AID %in% sample$AID),]; mut3 <- mut3[which(mut3$gene %in% c("chrXq","chrXp","chr22q11","chr21q","chr20q","chr18q21","chr16p","chr13q","chr12p","chr11q24","chr9q","chr9p","chr8q","chr8q13","chr8p","chr4q","chr3p23","chr2p24","chr1q","chr1q21","chr6q25","chr1p36")),]
geneCNV <- KITamp[which(KITamp$AID %in% sample$AID),]

new_table1 <- merge(mut1,sample,1,1)
new_table2 <- merge(mut2,sample,1,1);new_table2 <- new_table2[which(new_table2$LOH=="LOH"),]
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table1,aes(x=reorder(AID,order),y=gene,fill=Point),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+geom_point(data = geneCNV,aes(x=AID,y=gene),shape=17,color='#CD0000',size=2.6)
F1A.pLot<-F1A.pLot+geom_text(data = new_table2,aes(x=AID,y=gene),label="L",color='black',size=3.5)
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',missense='#66c2a5',inframeIndel='#a6d854',frameshift='#fc8d62',splice="#619CFF",new_splice="#619CFF",nonsense='#df65b0',f_NA='grey'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part1 <- F1A.pLot
Part1
gene.scale1 = 1 + (ncol(PointMut) - 14)/14

new_table3 <- merge(mut3,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table3,aes(x=reorder(AID,order),y=gene,fill=CNV),color='white',width=1,height=1,size=0.4,stat='identity')

F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Gain='#CD0000',Loss='#1400CE',f_NA='grey'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part2 <- F1A.pLot
Part2
gene.scale2 = 1 + (ncol(CNA) -1 - 14)/14

##############
Info <- Clinical[which(Clinical$AID %in% sample$AID),]
data1 <- Info[which(Info$gene %in% c("Domin.His","Site","Gender","Race","Seq")),]
data2 <- Info[which(Info$gene %in% c("WGD","Karyotype","GInstability")),]
data3 <- Info[which(Info$gene %in% c("Age")),]

new_table4 <- merge(data1,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table4,aes(x=reorder(AID,order),y=gene,fill=Info),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Germinoma="#FFA49C",IMT="#CBCE00",Teratoma="#1aa3d9",YST="#FF9FFF",EmbryonalC="#BDBAD7",f_NA='grey',Mix='grey',NGGCT='grey',No='grey',NANA='grey95',
                                                        Asian="#ffd92f",China="#ffd92f" ,Japan="#ffd92f", Black="#386cb0", White= "#80cdc1", Indian="#b3de69", White.Asian="#b3de69" , White.Black="#b3de69", White.Indian="#b3de69",WES="#93C6E7",WGS="#54278f",
                                                        Intracranial="#9970ab",Extracranial="#a6dba0",Ovary="#a6dba0",Pelvis="#a6dba0",Testis="#a6dba0",Male='#b3cde3',Female='#fccde5',Normal='grey95',Trisomy21="#f46d43",XO="#f46d43",XXY='#f46d43'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part4 <- F1A.pLot
Part4
gene.scale4 = 1 + (length(unique(new_table4$gene)) - 14)/14

new_table5 <- merge(data2,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table5,aes(x=reorder(AID,order),y=gene,fill=Info),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Mut="#dd3497",Haploid="#3A8891",f_NA='grey',Normal='grey95',Trisomy21="#9970ab",XO="#00B4FF",XXY='#5aae61',NANA='grey95',XY="#b15928"))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part3 <- F1A.pLot
Part3
gene.scale3 = 1 + (length(unique(new_table5$gene)) - 14)/14*0.99

new_table6 <- merge(data3,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_hline(yintercept = 3,color="#f46d43",linetype=2)
F1A.pLot<-F1A.pLot+geom_point(data = new_table6,aes(x=reorder(AID,order),y=as.numeric(as.character(Info))),shape=16,color='black',size=2)
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='grey98',color='transparent',size=1),plot.margin=unit(c(0.45,0.1,1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_text(size=10,angle = 90,vjust = 0.5,hjust = 1,face="plain",color="black"),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete()+xlab(NULL)+ylab(NULL)+scale_y_continuous(expand=c(0,0),limits=c(-1,21),breaks = seq(0,20,10) )
Part5 <- F1A.pLot
Part5
gene.scale5 = 1 + (3 - 14)/14

figure_2 <- rbind(ggplotGrob(Part1),ggplotGrob(Part2),ggplotGrob(Part3),ggplotGrob(Part4),ggplotGrob(Part5),size="first")
panels <- figure_2$layout$t[grep("panel", figure_2$layout$name)]
figure_2$heights[panels][1] <- unit(gene.scale1,'null')
figure_2$heights[panels][2] <- unit(gene.scale2,'null')
figure_2$heights[panels][3] <- unit(gene.scale3,'null')
figure_2$heights[panels][4] <- unit(gene.scale4,'null')
figure_2$heights[panels][5] <- unit(gene.scale5,'null')
ggsave(file="test_YST.pdf", plot=figure_2,bg = 'white', width = 12, height = 30, units = 'cm', dpi = 600)
#####################################################################################################################################################################


###MT
#####################################################################################################################################################################
sample <- as.data.frame(Sample_features[which(Sample_features$Histology=="Teratoma" | Sample_features$Histology=="IMT"),])
sortList <- read.table("../../lib/sort_Teratomas2.txt",header = T)
sample <- merge(sample,sortList,1,1); scaleW3 <- nrow(sample)

mut1 <- Mut[which(Mut$AID %in% sample$AID),]
mut2 <- LOH[which(LOH$AID %in% sample$AID),]
mut3 <- CNV[which(CNV$AID %in% sample$AID),]; mut3 <- mut3[which(mut3$gene %in% c("chrXq","chrXp","chr22q11","chr21q","chr20q","chr18q21","chr16p","chr13q","chr12p","chr11q24","chr9q","chr9p","chr8q","chr8q13","chr8p","chr4q","chr3p23","chr2p24","chr1q","chr1q21","chr6q25","chr1p36")),]
geneCNV <- KITamp[which(KITamp$AID %in% sample$AID),]

new_table1 <- merge(mut1,sample,1,1)
new_table2 <- merge(mut2,sample,1,1);new_table2 <- new_table2[which(new_table2$LOH=="LOH"),]
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table1,aes(x=reorder(AID,order),y=gene,fill=Point),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+geom_point(data = geneCNV,aes(x=AID,y=gene),shape=17,color='#CD0000',size=2.6)
F1A.pLot<-F1A.pLot+geom_text(data = new_table2,aes(x=AID,y=gene),label="L",color='black',size=3.5)
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',missense='#66c2a5',inframeIndel='#a6d854',frameshift='#fc8d62',splice="#619CFF",new_splice="#619CFF",nonsense='#df65b0',f_NA='grey'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part1 <- F1A.pLot
Part1
gene.scale1 = 1 + (ncol(PointMut) - 14)/14

new_table3 <- merge(mut3,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table3,aes(x=reorder(AID,order),y=gene,fill=CNV),color='white',width=1,height=1,size=0.4,stat='identity')

F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Gain='#CD0000',Loss='#1400CE',f_NA='grey'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part2 <- F1A.pLot
Part2
gene.scale2 = 1 + (ncol(CNA) -1 - 14)/14

##############
Info <- Clinical[which(Clinical$AID %in% sample$AID),]
data1 <- Info[which(Info$gene %in% c("Domin.His","Site","Gender","Race","Seq")),]
data2 <- Info[which(Info$gene %in% c("WGD","Karyotype","GInstability")),]
data3 <- Info[which(Info$gene %in% c("Age")),]

new_table4 <- merge(data1,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table4,aes(x=reorder(AID,order),y=gene,fill=Info),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Germinoma="#FFA49C",IMT="#CBCE00",Teratoma="#1aa3d9",YST="#FF9FFF",EmbryonalC="#BDBAD7",f_NA='grey',Mix='grey',NGGCT='grey',No='grey',NANA='grey95',
                                                        Asian="#ffd92f",China="#ffd92f" ,Japan="#ffd92f", Black="#386cb0", White= "#80cdc1", Indian="#b3de69", White.Asian="#b3de69" , White.Black="#b3de69", White.Indian="#b3de69",WES="#93C6E7",WGS="#54278f",
                                                        Intracranial="#9970ab",Extracranial="#a6dba0",Ovary="#a6dba0",Pelvis="#a6dba0",Testis="#a6dba0",Male='#b3cde3',Female='#fccde5',Normal='grey95',Trisomy21="#f46d43",XO="#f46d43",XXY='#f46d43'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part4 <- F1A.pLot
Part4
gene.scale4 = 1 + (length(unique(new_table4$gene)) - 14)/14

new_table5 <- merge(data2,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table5,aes(x=reorder(AID,order),y=gene,fill=Info),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Mut="#dd3497",Haploid="#3A8891",f_NA='grey',Normal='grey95',Trisomy21="#9970ab",XO="#00B4FF",XXY='#5aae61',NANA='grey95',XY="#b15928"))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part3 <- F1A.pLot
Part3
gene.scale3 = 1 + (length(unique(new_table5$gene)) - 14)/14*0.99

new_table6 <- merge(data3,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_hline(yintercept = 3,color="#f46d43",linetype=2)
F1A.pLot<-F1A.pLot+geom_point(data = new_table6,aes(x=reorder(AID,order),y=as.numeric(as.character(Info))),shape=16,color='black',size=2)
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='grey98',color='transparent',size=1),plot.margin=unit(c(0.45,0.1,1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_text(size=10,angle = 90,vjust = 0.5,hjust = 1,face="plain",color="black"),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete()+xlab(NULL)+ylab(NULL)+scale_y_continuous(expand=c(0,0),limits=c(-1,21),breaks = seq(0,20,10) )
Part5 <- F1A.pLot
Part5
gene.scale5 = 1 + (3 - 14)/14

figure_3 <- rbind(ggplotGrob(Part1),ggplotGrob(Part2),ggplotGrob(Part3),ggplotGrob(Part4),ggplotGrob(Part5),size="first")
panels <- figure_3$layout$t[grep("panel", figure_3$layout$name)]
figure_3$heights[panels][1] <- unit(gene.scale1,'null')
figure_3$heights[panels][2] <- unit(gene.scale2,'null')
figure_3$heights[panels][3] <- unit(gene.scale3,'null')
figure_3$heights[panels][4] <- unit(gene.scale4,'null')
figure_3$heights[panels][5] <- unit(gene.scale5,'null')
ggsave(file="test_MT.pdf", plot=figure_3,bg = 'white', width = 12, height = 30, units = 'cm', dpi = 600)
#####################################################################################################################################################################


###Mix
#####################################################################################################################################################################
sample <- as.data.frame(Sample_features[which(Sample_features$Histology=="Mix" | Sample_features$Histology=="EmbryonalC"),])
sortList <- read.table("../../lib/sort_Mix.txt",header = T)
sample <- merge(sample,sortList,1,1); scaleW5 <- nrow(sample)

mut1 <- Mut[which(Mut$AID %in% sample$AID),]
mut2 <- LOH[which(LOH$AID %in% sample$AID),]
mut3 <- CNV[which(CNV$AID %in% sample$AID),]; mut3 <- mut3[which(mut3$gene %in% c("chrXq","chrXp","chr22q11","chr21q","chr20q","chr18q21","chr16p","chr13q","chr12p","chr11q24","chr9q","chr9p","chr8q","chr8q13","chr8p","chr4q","chr3p23","chr2p24","chr1q","chr1q21","chr6q25","chr1p36")),]
geneCNV <- KITamp[which(KITamp$AID %in% sample$AID),]

new_table1 <- merge(mut1,sample,1,1)
new_table2 <- merge(mut2,sample,1,1);new_table2 <- new_table2[which(new_table2$LOH=="LOH"),]
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table1,aes(x=reorder(AID,order),y=gene,fill=Point),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+geom_point(data = geneCNV,aes(x=AID,y=gene),shape=17,color='#CD0000',size=2.6)
F1A.pLot<-F1A.pLot+geom_text(data = new_table2,aes(x=AID,y=gene),label="L",color='black',size=3.5)
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',missense='#66c2a5',inframeIndel='#a6d854',frameshift='#fc8d62',splice="#619CFF",new_splice="#619CFF",nonsense='#df65b0',f_NA='grey'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part1 <- F1A.pLot
Part1
gene.scale1 = 1 + (ncol(PointMut) - 14)/14

new_table3 <- merge(mut3,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table3,aes(x=reorder(AID,order),y=gene,fill=CNV),color='white',width=1,height=1,size=0.4,stat='identity')

F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Gain='#CD0000',Loss='#1400CE',f_NA='grey'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part2 <- F1A.pLot
Part2
gene.scale2 = 1 + (ncol(CNA) -1 - 14)/14

##############
Info <- Clinical[which(Clinical$AID %in% sample$AID),]
data1 <- Info[which(Info$gene %in% c("Domin.His","Site","Gender","Race","Seq")),]
data2 <- Info[which(Info$gene %in% c("WGD","Karyotype","GInstability")),]
data3 <- Info[which(Info$gene %in% c("Age")),]

new_table4 <- merge(data1,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table4,aes(x=reorder(AID,order),y=gene,fill=Info),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Germinoma="#FFA49C",IMT="#CBCE00",Teratoma="#1aa3d9",YST="#FF9FFF",EmbryonalC="#BDBAD7",f_NA='grey',Mix='grey',NGGCT='grey',No='grey',NANA='grey95',
                                                        Asian="#ffd92f",China="#ffd92f" ,Japan="#ffd92f", Black="#386cb0", White= "#80cdc1", Indian="#b3de69", White.Asian="#b3de69" , White.Black="#b3de69", White.Indian="#b3de69",WES="#93C6E7",WGS="#54278f",
                                                        Intracranial="#9970ab",Extracranial="#a6dba0",Ovary="#a6dba0",Pelvis="#a6dba0",Testis="#a6dba0",Male='#b3cde3',Female='#fccde5',Normal='grey95',Trisomy21="#f46d43",XO="#f46d43",XXY='#f46d43'))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part4 <- F1A.pLot
Part4
gene.scale4 = 1 + (length(unique(new_table4$gene)) - 14)/14

new_table5 <- merge(data2,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table5,aes(x=reorder(AID,order),y=gene,fill=Info),color='white',width=1,height=1,size=0.4,stat='identity')
F1A.pLot<-F1A.pLot+scale_fill_manual(name=NULL,values=c(WT='grey95',Mut="#dd3497",Haploid="#3A8891",f_NA='grey',Normal='grey95',Trisomy21="#9970ab",XO="#00B4FF",XXY='#5aae61',NANA='grey95',XY="#b15928"))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)
Part3 <- F1A.pLot
Part3
gene.scale3 = 1 + (length(unique(new_table5$gene)) - 14)/14*0.99

new_table6 <- merge(data3,sample,1,1)
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_hline(yintercept = 3,color="#f46d43",linetype=2)
F1A.pLot<-F1A.pLot+geom_point(data = new_table6,aes(x=reorder(AID,order),y=as.numeric(as.character(Info))),shape=16,color='black',size=2)
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='grey98',color='transparent',size=1),plot.margin=unit(c(0.45,0.1,1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_text(size=10,angle = 90,vjust = 0.5,hjust = 1,face="plain",color="black"),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete()+xlab(NULL)+ylab(NULL)+scale_y_continuous(expand=c(0,0),limits=c(-1,21),breaks = seq(0,20,10) )
Part5 <- F1A.pLot
Part5
gene.scale5 = 1 + (3 - 14)/14

figure_5 <- rbind(ggplotGrob(Part1),ggplotGrob(Part2),ggplotGrob(Part3),ggplotGrob(Part4),ggplotGrob(Part5),size="first")
panels <- figure_5$layout$t[grep("panel", figure_5$layout$name)]
figure_5$heights[panels][1] <- unit(gene.scale1,'null')
figure_5$heights[panels][2] <- unit(gene.scale2,'null')
figure_5$heights[panels][3] <- unit(gene.scale3,'null')
figure_5$heights[panels][4] <- unit(gene.scale4,'null')
figure_5$heights[panels][5] <- unit(gene.scale5,'null')
ggsave(file="test_Mix.pdf", plot=figure_5,bg = 'white', width = 12, height = 30, units = 'cm', dpi = 600)
#####################################################################################################################################################################



figure<-cbind(figure_1,figure_2,figure_3,figure_5)

panels <- figure$layout$t[grep("panel", figure$layout$name)]

figure$widths[7]  <- unit(scaleW1/34,'null')
figure$widths[20] <- unit(scaleW2/34,'null')
figure$widths[33] <- unit(scaleW3/34,'null')
figure$widths[46] <- unit(scaleW5/34,'null')
ggsave(file="somatic9.pdf", plot=figure,bg = 'white', width =50, height = 32, units = 'cm', dpi = 600)





