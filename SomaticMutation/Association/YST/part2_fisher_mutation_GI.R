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
Clinical <- Clinical[which(Clinical$Seq=="WGS" & !(Clinical$AID %in% c("F028","F030","F032","F033","F035","F049"))),]
#Clinical <- Clinical[which(Clinical$Seq=="WGS" ),]

PointMut <- read.table("../../Landscape/results/Final_PointMut.txt",header=T)
Arm <- read.table("../../Landscape/results/Final_Arm.txt",header=T)
Focal <- read.table("../../Landscape/results/Final_Focal.txt",header=T)
WGD <- read.table("../../Landscape/results/Final_WGD.txt",header=T)
SV <- read.table("../../Landscape/results/Final_chromoplexy.txt",header=T)

Signature <- read.table("../../Landscape/results/Final_contribution.txt",header=T)

PointMut$MAPK <- "WT"; 
PointMut$MAPK[which(PointMut$KIT!="WT" | PointMut$KRAS!="WT" | PointMut$RRAS2!="WT" | PointMut$NRAS!="WT" |
                      PointMut$CBL!="WT" | PointMut$FGFR2!="WT" | PointMut$NF1!="WT" | PointMut$SOS1!="WT" | PointMut$USP28!="WT")] <- "Mut"
PointMut$PI3K <- "WT"; 
PointMut$PI3K[which(PointMut$PTEN!="WT" | PointMut$MTOR!="WT" | PointMut$AKT3!="WT" | PointMut$PIK3CA!="WT" |
                      PointMut$PIK3CB!="WT" | PointMut$PIK3CD!="WT")] <- "Mut"
PointMut$Chromatin <- "WT"; 
PointMut$Chromatin[which(PointMut$BCOR!="WT" | PointMut$BCORL1!="WT" | PointMut$ARID1A!="WT" | PointMut$ARID1B!="WT" | PointMut$ARID4B!="WT")] <- "Mut"

PointMut$DNArepair <- "WT"; 
PointMut$DNArepair[which(PointMut$BRCA1!="WT" | PointMut$BLM!="WT" | PointMut$ATM!="WT" | PointMut$CHEK2!="WT")] <- "Mut"

Interaction <- merge(Clinical,PointMut)
Interaction <- merge(Interaction,Arm)
Interaction <- merge(Interaction,Focal)
Interaction <- merge(Interaction,WGD)
Signature <- merge(Signature,SV)



Interaction <- Interaction[which(Interaction$Histology=="YST"),]
AID <- Interaction$AID

#######
#############

##################################
x <- c()
for(i in c(11:39)){
  if(length(which(Interaction[,i] != "WT")) >= 2){
    x <- c(x,i)
  }
}
Mut1 <- Interaction[,x];Mut<- Mut1; Mut[,] <- 0
Mut[Mut1!="WT"] <- 1 

##################################
x <- c()
for(i in c(40:63)){
  if(length(which(Interaction[,i] == "Gain")) >= 2){
    x <- c(x,i)
  }
}
Gain1 <- Interaction[,x];Gain<- Gain1; Gain[,] <- 0
Gain[Gain1=="Gain"] <- 1 ; colnames(Gain) <- paste0(colnames(Gain),".Gain")

##################################
x <- c()
for(i in c(40:63)){
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

Genome <- data.frame(AID,WGD)

GenomicAlter <- cbind(Genome,Mut,Gain,Loss)

#step2: conduct chisq.test for every-pair features and summarize the pvalue results into a table

GenomicAlter <- GenomicAlter[,c( "AID", "WGD","chrXq.Gain", "chr22q11.Gain", "chr21q.Gain", "chr20q.Gain","chr16p.Loss",  "chr12p.Gain", 
                                 "chr8p.Loss","chr6q25.Loss","chr4q.Loss","chr3p23.Gain","chr2p24.Gain","chr1q21.Gain", "chr1p36.Loss","RRAS2")]

scale_vect = function(x){
  m = mean(x, na.rm = T)
  s = sd(x, na.rm = T)
  return((x - m) / s)
}

Signature <- Signature[,c("AID","SBSgct","SBS18","ComplexSV")]
temp <- Signature[which(Signature$AID %in% Clinical$AID),]
Signature <- data.frame(temp$AID,rep("No",nrow(temp)),rep("No",nrow(temp)),rep("No",nrow(temp)) )
colnames(Signature) <- colnames(temp)
Signature$SBSgct[which(temp$SBSgct>0)] <- "Yes"
Signature$SBS18[which(temp$SBS18>0)] <- "Yes"
Signature$ComplexSV[which(temp$ComplexSV!="No")] <- "Yes"
#Signature$Chromoplexy[which(temp$Chromoplexy=="Chromoplexy")] <- "Yes"


nfeaA <- ncol(Signature)
nfeaB <- ncol(GenomicAlter)
myTable <- data.frame(rep(0,(nfeaA-1)*(nfeaB-1)),rep(0,(nfeaA-1)*(nfeaB-1)),rep(0,(nfeaA-1)*(nfeaB-1)),rep(0,(nfeaA-1)*(nfeaB-1)),rep(0,(nfeaA-1)*(nfeaB-1)),rep(0,(nfeaA-1)*(nfeaB-1)) )
colnames(myTable) <- c("VarA","VarB","OR","pValue","Loss","OddsR")
n=1
for(i in c(2:nfeaA)){
  for(j in c(2:nfeaB)){
    tmpA <- Signature[,c(1,i)]
    tmpB <- GenomicAlter[,c(1,j)]
    tmpC <- merge(tmpA,tmpB); colnames(tmpC) <- c("ID","Var1","Var2")
    tmpC <- tmpC[which(tmpC$Var1>=0 & (!is.na(tmpC$Var2)) ),]
    tmpD <- data.frame(c(length(which(tmpC$Var1=="Yes" & tmpC$Var2==1)),length(which(tmpC$Var1=="Yes" & tmpC$Var2==0)) ),c(length(which(tmpC$Var1=="No" & tmpC$Var2==1)),length(which(tmpC$Var1=="No" & tmpC$Var2==0))))
    colnames(tmpD) <- c("Yes","No")
    testing <- fisher.test(tmpD,alternative="two.sided")
    if(tmpD[1,1] == 0 | tmpD[1,2] == 0 | tmpD[2,1] == 0 | tmpD[2,2] == 0){
      OddsR <- ((tmpD[1,1]+0.5)*(tmpD[2,2]+0.5)) / ((tmpD[1,2]+0.5)*(tmpD[2,1]+0.5))
    }
    else{
      OddsR <- (tmpD[1,1]*tmpD[2,2]) / (tmpD[1,2]*tmpD[2,1])
    }
    myTable$OddsR[n] <- testing$estimate
    myTable$OR[n] <- log2(OddsR)
    myTable$pValue[n] <- testing$p.value
    myTable$Loss[n] <- length(which(tmpC$Var2==1))
    myTable$Neutral[n] <- length(which(tmpC$Var2==0))
    myTable$VarA[n] <- colnames(Signature)[i]
    myTable$VarB[n] <- colnames(GenomicAlter)[j]

    n<-n+1
  }
}

#myTable$qValue <- p.adjust(myTable$pValue,method = "bonferroni")
myTable$fdr <- p.adjust(myTable$pValue,method = "fdr")



myTable$Dotsize <- -log10(myTable$pValue)

myTable$Dotcolor <- "NS"
myTable$Dotcolor[which(myTable$OR > 0 & myTable$pValue <= 0.05)] <- "Cooccur" 
myTable$Dotcolor[which(myTable$OR < 0 & myTable$pValue <= 0.05)] <- "Exclusive"  


#step3: draw the figure
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}
plot_1.data<-myTable
plot_1.data$orderID<-c(1:nrow(plot_1.data))

#plot_1.data$Dotsize[which(plot_1.data$pValue >= 0.1)] <- 0.25
plot_2.data <- plot_1.data[which(plot_1.data$pValue >= 0.1),]
plot_3.data <- plot_1.data[which(plot_1.data$pValue < 0.1),]
plot_4.data <- plot_1.data[which(plot_1.data$pValue < 0.05),]


coMu.plot<-ggplot()+theme_classic()
coMu.plot<-coMu.plot+geom_vline(xintercept = c(1:42),linetype=2,color="grey",size=0.3)
coMu.plot<-coMu.plot+geom_hline(yintercept = plot_1.data$VarB,linetype=2,color="grey",size=0.3)
coMu.plot<-coMu.plot+geom_point(data = plot_1.data,aes(x=reorder(VarA,orderID),y=reorder(VarB,orderID),fill=OR),size=0,shape=NA)+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot<-coMu.plot+geom_point(data = plot_2.data,aes(x=reorder(VarA,orderID),y=reorder(VarB,orderID),color=OR),size=2,shape=15,stroke=0)+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot<-coMu.plot+geom_point(data = plot_3.data,aes(x=reorder(VarA,orderID),y=reorder(VarB,orderID),size=Dotsize,color=OR),shape=16,stroke=0)+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot<-coMu.plot+geom_point(data = plot_4.data,aes(x=reorder(VarA,orderID),y=reorder(VarB,orderID),size=Dotsize),color="black",shape=1,stroke=1)+ylab(NULL)+xlab(NULL)+ggtitle(NULL)

coMu.plot<-coMu.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(1,6,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.2,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='top',
                           legend.text=element_text(size=12,face='plain'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_text(size=14,vjust=0.5,hjust=0,face='plain'),
                           axis.text.x=element_text(size=14,angle=45,vjust=0,hjust=0,face='plain',color='black'),axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='plain',color='black'),
                           axis.title.y=element_text(size=20,hjust=0.5,vjust=2,face='plain',color='black'),axis.line = element_blank(),
                           strip.text = element_text(size=18,face='bold',vjust=0.5,hjust=0.5),strip.background = element_rect(colour="black", fill=gg_color_hue(3)))

coMu.plot<-coMu.plot+scale_size(name   = "Fisher's test\n   p-value", range = c(2.9,4),trans = "exp",  breaks = c( 0.7,1,1.3,1.6,2,3,4),labels = expression(0.2,0.1,0.05,0.025,10^-2,10^-3, 10^-4),guide = guide_legend(nrow=2))
#coMu.plot<-coMu.plot+scale_size(name   = "Fisher's test\n   q-value", range = c(0.5,8),trans = "hms",guide = guide_legend(nrow=3))

#coMu.plot<-coMu.plot+scale_fill_gradientn(name=NULL,colours=c("#018571","#80cdc1","#f7f7f7","#dfc27d","#a6611a"),limits=c(-1.75,1.75),breaks=c(-1.6,-0.8,0,0.8,1.6))
#coMu.plot<-coMu.plot+scale_color_gradientn(name=NULL,colours=c("#018571","#80cdc1","#f7f7f7","#dfc27d","#a6611a"),limits=c(-6.8,6.8),breaks=seq(-6,6,3))
coMu.plot<-coMu.plot+scale_color_gradientn(name="Odds ratio\n(log2)",colours=c('#4d9221','#b8e186','#f7f7f7','#f1b6da','#c51b7d'),limits=c(-6.8,6.8),breaks=seq(-6,6,3))
coMu.plot<-coMu.plot+scale_fill_gradientn(name="Odds ratio\n(log2)",colours=c('#4d9221','#b8e186','#f7f7f7','#f1b6da','#c51b7d'),limits=c(-6.8,6.8),breaks=seq(-6,6,3))

coMu.plot<-coMu.plot+scale_x_discrete(position = "top")



#'#5e3c99','#b2abd2',"grey95",'#fdb863','#e66101'





plot1 <- coMu.plot

figure_1<-rbind(ggplotGrob(plot1),size="last")
ggsave(file="./Fisher_Mutation_Mutation_YST_part2.pdf", plot=figure_1,bg = 'white', width = 9.3, height = 16, units = 'cm', dpi = 600)





