 # WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


Clinical <- read.table("../../../lib/Clinical.txt",sep = "\t",header = T)
PointMut <- read.table("../../Landscape/results/Final_PointMut.txt",header=T)
Arm <- read.table("../../Landscape/results/Final_Arm.txt",header=T)
Focal <- read.table("../../Landscape/results/Final_Focal.txt",header=T)
WGD <- read.table("../../Landscape/results/Final_WGD.txt",header=T)
Number <- read.table("../../Landscape/results/Final_Number.txt",header=T)

PointMut$MAPK <- "WT"; 
PointMut$MAPK[which(PointMut$KIT!="WT" | PointMut$KRAS!="WT" | PointMut$RRAS2!="WT" | PointMut$NRAS!="WT" |
                      PointMut$CBL!="WT" | PointMut$FGFR2!="WT" | PointMut$NF1!="WT" | PointMut$SOS1!="WT")] <- "Mut"
PointMut$PI3K <- "WT"; 
PointMut$PI3K[which(PointMut$PTEN!="WT" | PointMut$MTOR!="WT" | PointMut$AKT3!="WT" | PointMut$PIK3CA!="WT" |
                      PointMut$PIK3CB!="WT" | PointMut$PIK3CD!="WT")] <- "Mut"
PointMut$DNArepair <- "WT"; 
PointMut$DNArepair[which(PointMut$BRCA1!="WT" | PointMut$BLM!="WT" | PointMut$ATM!="WT" | PointMut$CHEK2!="WT")] <- "Mut"

Interaction <- merge(Clinical,PointMut)
Interaction <- merge(Interaction,Arm)
Interaction <- merge(Interaction,Focal)
Interaction <- merge(Interaction,WGD)
Interaction <- merge(Interaction,Number)


Interaction <- Interaction[which(Interaction$Histology=="YST"),]
AID <- Interaction$AID



##################################
x <- c()
for(i in c(11:38)){
  if(length(which(Interaction[,i] != "WT")) >= 2){
    x <- c(x,i)
  }
}
Mut1 <- Interaction[,x];Mut<- Mut1; Mut[,] <- 0
Mut[Mut1!="WT"] <- 1 

##################################
x <- c()
for(i in c(39:62)){
  if(length(which(Interaction[,i] == "Gain")) >= 2){
    x <- c(x,i)
  }
}
Gain1 <- Interaction[,x];Gain<- Gain1; Gain[,] <- 0
Gain[Gain1=="Gain"] <- 1 ; colnames(Gain) <- paste0(colnames(Gain),".Gain")

##################################
x <- c()
for(i in c(39:62)){
  if(length(which(Interaction[,i] == "Loss")) >= 2){
    x <- c(x,i)
  }
}
Loss1 <- Interaction[,x];Loss<- Loss1; Loss[,] <- 0
Loss[Loss1=="Loss"] <- 1 ; colnames(Loss) <- paste0(colnames(Loss),".Loss")
##################

##################
WGD <- rep(0,nrow(Interaction))
WGD[which(Interaction$WGD == "Mut")] = 1


High_CNA <- rep(NA,nrow(Interaction))
all<-Interaction$TotalCNA[which(Interaction$TotalCNA>=0)]
cutoff <- max(summary(all)[[4]],summary(all)[[3]])
High_CNA[which(Interaction$TotalCNA > cutoff)] = 1
High_CNA[which(Interaction$TotalCNA <= cutoff & Interaction$TotalCNA >= 0)] = 0

High_TE <- rep(NA,nrow(Interaction))
all<-Interaction$TotalTE[which(Interaction$TotalTE>=0)]
cutoff <- max(summary(all)[[4]],summary(all)[[3]])
High_TE[which(Interaction$TotalTE > cutoff)] = 1
High_TE[which(Interaction$TotalTE <= cutoff & Interaction$TotalTE >= 0)] = 0

High_SV <- rep(NA,nrow(Interaction))
all<-Interaction$TotalSV[which(Interaction$TotalSV>=0)]
cutoff <- max(summary(all)[[4]],summary(all)[[3]])
High_SV[which(Interaction$TotalSV > cutoff)] = 1
High_SV[which(Interaction$TotalSV <= cutoff & Interaction$TotalSV >= 0)] = 0

High_SNV <- rep(NA,nrow(Interaction))
all<-Interaction$TotalSNV[which(Interaction$TotalSNV>=0)]
cutoff <- max(summary(all)[[4]],summary(all)[[3]])
High_SNV[which(Interaction$TotalSNV > cutoff)] = 1
High_SNV[which(Interaction$TotalSNV <= cutoff &  Interaction$TotalSNV >= 0)] = 0

SBS18_or_GCT <- rep(NA,nrow(Interaction))
SBS18_or_GCT[which(Interaction$SBS18_or_SBSgct == "Mut")] = 1
SBS18_or_GCT[which(Interaction$SBS18_or_SBSgct == "WT")] = 0


Genome <- data.frame(AID,WGD)

GenomicAlter <- cbind(Genome,Mut,Gain,Loss)
GenomicAlter <- GenomicAlter[,which(!(colnames(GenomicAlter) %in% c("chrXp.Gain","chr1q.Gain","MAPK")))]
library(reshape2)
myTable <- melt(GenomicAlter); myTable <- myTable[which(myTable$value==1),]
part1 <- as.data.frame(table(myTable$variable)); colnames(part1) <- c("Mut","Freq")

myTable <- melt(GenomicAlter); myTable <- myTable[which(myTable$value==1 | myTable$value==0),]
part2 <- as.data.frame(table(myTable$variable)); colnames(part2) <- c("Mut","Case")
myTable <- merge(part1,part2,1,1)
myTable$Frac <- myTable$Freq/myTable$Case

#step3: draw the figure
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}
plot_1.data<- myTable[which(myTable$Freq>=3),]
plot_2.data <- plot_1.data[which(plot_1.data$Frac < 0.5),]
plot_3.data <- plot_1.data[which(plot_1.data$Frac >= 0.5),]
coMu.plot<-ggplot()+theme_classic()
coMu.plot<-coMu.plot+geom_bar(data = plot_1.data,aes(x=reorder(Mut,-Frac),y=Frac),alpha=0.7,fill="transparent",size=0.6,width=0.5,stat='identity',position=position_dodge2()) +
  ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot<-coMu.plot+geom_bar(data = plot_2.data,aes(x=reorder(Mut,-Frac),y=Frac),alpha=0.9,fill="grey40",size=0.6,width=0.5,stat='identity',position=position_dodge2()) 
coMu.plot<-coMu.plot+geom_bar(data = plot_3.data,aes(x=reorder(Mut,-Frac),y=Frac),alpha=0.9,fill="#2166ac",size=0.6,width=0.5,stat='identity',position=position_dodge2()) 

#coMu.plot<-coMu.plot+geom_point(data = plot_2.data,aes(x=reorder(VarB,qValue),y=Dotsize,size=OddsR,color=Dotcolor),shape=16)+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot<-coMu.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,2,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='top',
                           legend.text=element_text(size=12,face='plain'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_text(size=14,vjust=0.5,hjust=0,face='plain'),
                           axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='plain',color='black'),
                           axis.title.y=element_text(size=20,hjust=0.5,vjust=2,face='plain',color='black'),#axis.line = element_blank(),
                           strip.text = element_text(size=18,face='bold',vjust=0.5,hjust=0.5),strip.background = element_rect(colour="black", fill=gg_color_hue(3)))
coMu.plot<-coMu.plot+scale_fill_manual(name=NULL,values =c('#c51b7d','#4d9221'))
coMu.plot<-coMu.plot+scale_x_discrete(position = "bottom")+scale_y_continuous(expand = c(0,0),limits = c(0,0.68))

plot1 <- coMu.plot

figure_1<-rbind(ggplotGrob(plot1),size="last")
ggsave(file="./ranking_YST_PPT.pdf", plot=figure_1,bg = 'white', width = 12.5, height = 8, units = 'cm', dpi = 600)










