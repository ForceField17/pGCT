# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(stringr)
library(fgsea)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
library(karyoploteR)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}


YST <- c("C005","F007","F017","F026","F029","F036","F037","F038","F045","F055","F056","SJ01","SJ10","SJ25","SJ32")
myTable <- read.table("./SBSyst_distribution.txt",header = F)
myTable$rank <- as.numeric(sub("X","23",sub("chr","",myTable$V3)))
myTable$type <- "Others"
myTable$type[which(myTable$V2 %in% YST)] <- "YST"

####plot
F1a.plot<-ggplot(data=myTable,aes(x=reorder(V3,rank),y=as.integer(V4)/1000000,color=type))
F1a.plot<-F1a.plot+ theme_classic()#geom_boxplot(data=myTable,aes(x=Dominant.His,y=Exp,fill=Dominant.His),width=0.3,alpha=0.2,outlier.shape = NA)+
F1a.plot<-F1a.plot+ geom_quasirandom(data=myTable,aes(x=reorder(V3,rank),y=-as.integer(V4)/1000000,color=type),width = 0.4,size=0.7,alpha=0.9,varwidth = T,bandwidth = 0.001,position = "identity")
#F1a.plot<-F1a.plot+ geom_quasirandom(data=table1,aes(x=reorder(V1,rank),y=as.integer(V2),color=type),width = 0.5,size=0.5,alpha=0.9,varwidth = T,bandwidth = 0.01,position = "identity")
#F1a.plot<-F1a.plot+ geom_quasirandom(data=table2,aes(x=reorder(V1,rank),y=as.integer(V2),color=type),width = 0.2,size=0.5,alpha=0.7,varwidth = T,bandwidth = 0.001,position = "identity")
F1a.plot<-F1a.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),axis.line.x = element_blank(),
                         plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="top",legend.text=element_text(size=14,hjust=0,face='bold'),
                         axis.text.x=element_text(size=12,angle=35,vjust=1,hjust=1,color="black",face = "plain"),axis.text.y=element_text(size=12,face='plain',color='black'),
                         axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab(NULL)+ylab("Chromosome coordinate (Mb)")+ scale_y_continuous(expand=c(0,0),limits=c(-250,1),breaks = seq(-250,0,50))#+coord_flip()
#the <- compare_means(Exp ~ DNMs,  data = myTable)
#my_comparisons <- list( c("chr1", "chr4"),c("chr22","YST"),c("Teratoma","YST"),c("IMT","Teratoma"),c("Germinoma","IMT"),c("Germinoma","Teratoma"))
#F1a.plot<-F1a.plot+stat_compare_means(comparisons = my_comparisons)

F1a.plot<-F1a.plot + scale_fill_manual(name=NULL,values=c('#fc8d62','#66c2a5'))
                                       
F1a.plot<-F1a.plot + scale_color_manual(name=NULL,values=c('#fc8d62','#66c2a5'))
                                       

figure_2<-rbind(ggplotGrob(F1a.plot),size="last")

ggsave(file=paste0("chromosome_SBSyst.pdf"), plot=figure_2,bg = 'white', width = 40, height = 16, units = 'cm', dpi = 600)




pdf(file = paste0("chromosome_SBSyst_chromosome.pdf"),width=7,height=20, onefile = T)

table1 <- myTable[which(myTable$V9=="A[C-G]G"),]
df1 <- data.frame(table1$V3,table1$V4,table1$V4+1,rep("+",nrow(table1)),c(1:nrow(table1)))
colnames(df1) <- c("chr","start","end","strand","gene_id")
gr1 <- GenomicRanges::makeGRangesFromDataFrame(df1)
gr1

table2 <- myTable[which(myTable$V9=="A[C-T]G"),]
df2 <- data.frame(table2$V3,table2$V4,table2$V4+1,rep("+",nrow(table2)),c(1:nrow(table2)))
colnames(df2) <- c("chr","start","end","strand","gene_id")
gr2 <- GenomicRanges::makeGRangesFromDataFrame(df2)
gr2

kp <- plotKaryotype(genome = "hg38",chromosomes="canonical",cytobands=NULL,plot.type=2)
kp <- kpPlotDensity(kp, gr1,ymax=3,col="#fdb863",lwd=0.5,data.panel=1,window.size=1000000)
kp <- kpPlotDensity(kp, gr2,ymax=3,col="#80cdc1",lwd=0.5,data.panel=2,window.size=1000000)
dev.off() 




pdf(file = paste0("chromosome_SBSyst_chromosome_histology.pdf"),width=12,height=16, onefile = T)

table1 <- myTable[which(myTable$type=="YST"),]
df1 <- data.frame(table1$V3,table1$V4,table1$V4+1,rep("+",nrow(table1)),c(1:nrow(table1)))
colnames(df1) <- c("chr","start","end","strand","gene_id")
gr1 <- GenomicRanges::makeGRangesFromDataFrame(df1)
gr1

table2 <- myTable[which(myTable$type!="YST"),]
df2 <- data.frame(table2$V3,table2$V4,table2$V4+1,rep("+",nrow(table2)),c(1:nrow(table2)))
colnames(df2) <- c("chr","start","end","strand","gene_id")
gr2 <- GenomicRanges::makeGRangesFromDataFrame(df2)
gr2

kp <- plotKaryotype(genome = "hg38",chromosomes="canonical",cytobands=NULL,plot.type=2)
kp <- kpPlotDensity(kp, gr1,ymax=3,col="#FF9FFF",lwd=0.5,data.panel=1,window.size=1000000)
kp <- kpPlotDensity(kp, gr2,ymax=3,col="#80cdc1",lwd=0.5,data.panel=2,window.size=1000000)
dev.off() 




YST <- c("C005","F007","F017","F026","F029","F036","F037","F038","F045","F055","F056","SJ01","SJ10","SJ25","SJ32")
myTable <- read.table("./SBSall_distribution.txt",header = F)
myTable$rank <- as.numeric(sub("X","23",sub("chr","",myTable$V3)))
myTable$type <- "NonYST"
myTable$type[which(myTable$V2 %in% YST)] <- "YST"
myTable$SBS <- myTable$V9
myTable$SBS[which(!(myTable$V9 %in% c("A[C-G]G","A[C-T]G")))] <- "Others"

xxx <- myTable[,c("SBS","type")]
a<-table(xxx)


table2 <- read.table("box_SBSyst.txt",header = T)

table2 <- table2 %>%
  arrange(-as.integer(His),-as.integer(SBSgct)) %>%
  mutate(Pos = cumsum(Freq) - 0.5*Freq) %>%
  group_by(His) %>%
  mutate(frac=Freq/sum(Freq))
table2

table2$rank <- c(4:1)

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=reorder(His,rank),y=frac,fill=SBSgct),alpha=0.95,width=0.65, color = "black",stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.005),breaks = seq(0,1,0.2))
#revp_plot<-revp_plot+geom_text(data = table2,aes(x=reorder(Site,-Freq),y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,2,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="top",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=25,vjust=1,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Propotion of tumors")+scale_x_discrete(position = "bottom")#+coord_flip()
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(YES='#ff7f00',NO='grey90'))  #c(Male='#41b6c4',Female='#e78ac3'))
revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="Fisher_SBSgct.pdf", plot=figure_2,bg = 'white', width = 5, height = 10, units = 'cm', dpi = 600)
Mu.FEtest <- cbind(c(0,0),c(0,0))
Mu.FEtest[1,1] <- table2$Freq[1]
Mu.FEtest[1,2] <- table2$Freq[2]
Mu.FEtest[2,1] <- table2$Freq[3]
Mu.FEtest[2,2] <- table2$Freq[4]
fisher.test(Mu.FEtest,alternative ="two.sided")



table2 <- read.table("box_SBS18.txt",header = T)

table2 <- table2 %>%
  arrange(-as.integer(His),-as.integer(SBS18)) %>%
  mutate(Pos = cumsum(Freq) - 0.5*Freq) %>%
  group_by(His) %>%
  mutate(frac=Freq/sum(Freq))
table2

table2$rank <- c(4:1)

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=reorder(His,rank),y=frac,fill=SBS18),alpha=0.95,width=0.65, color = "black",stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.005),breaks = seq(0,1,0.2))
#revp_plot<-revp_plot+geom_text(data = table2,aes(x=reorder(Site,-Freq),y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,2,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="top",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=25,vjust=1,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Propotion of tumors")+scale_x_discrete(position = "bottom")#+coord_flip()
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(YES='#de77ae',NO='grey90'))  #c(Male='#41b6c4',Female='#e78ac3'))
revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="Fisher_SBS18.pdf", plot=figure_2,bg = 'white', width = 5, height = 10, units = 'cm', dpi = 600)
Mu.FEtest <- cbind(c(0,0),c(0,0))
Mu.FEtest[1,1] <- table2$Freq[1]
Mu.FEtest[1,2] <- table2$Freq[2]
Mu.FEtest[2,1] <- table2$Freq[3]
Mu.FEtest[2,2] <- table2$Freq[4]
fisher.test(Mu.FEtest,alternative ="two.sided")