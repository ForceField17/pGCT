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


Clinical <- read.table("../../lib/Clinical.txt",sep = "\t",header = T)


PointMut <- read.table("../Landscape/results/Final_PointMut.txt",header=T)
Arm <- read.table("../Landscape/results/Final_Arm.txt",header=T)
Focal <- read.table("../Landscape/results/Final_Focal.txt",header=T)
WGD <- read.table("../Landscape/results/Final_WGD.txt",header=T)
Number <- read.table("../Landscape/results/Final_Number.txt",header=T)

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
Interaction <- merge(Interaction,Number)

AID <- Interaction$AID

#######
Germinoma <- rep(NA,nrow(Interaction))
Germinoma[which(Interaction$Histology == "Germinoma")] = 1
Germinoma[which(Interaction$Histology != "Germinoma" & Interaction$Domin.His != "Germinoma" )] = 0

YST <- rep(NA,nrow(Interaction))
YST[which(Interaction$Histology == "YST")] = 1
YST[which(Interaction$Histology != "YST" & Interaction$Domin.His != "YST" )] = 0

Mix <- rep(0,nrow(Interaction))
Mix[which(Interaction$Histology == "Mix")] = 1


NGGCT <- rep(NA,nrow(Interaction))
NGGCT[which(Interaction$Subtype == "Germinoma" | Interaction$Subtype == "MixGE")] = 0
NGGCT[which(Interaction$Subtype == "NGGCT" & Interaction$Domin.His != "Mix" & Interaction$Domin.His != "NGGCT")] = 1

Teratomas <- rep(NA,nrow(Interaction))
Teratomas[which(Interaction$Histology == "Teratoma" | Interaction$Histology == "IMT")] = 1
Teratomas[which(Interaction$Histology != "Teratoma" & Interaction$Histology != "IMT" & Interaction$Domin.His != "Teratoma" & Interaction$Domin.His != "IMT" )] = 0



Histology <- data.frame(AID,Mix,Germinoma,Teratomas,YST)
#############

#############
Infancy <- rep(0,nrow(Interaction))
Infancy[which(Interaction$Age < 3)] = 1

Preschooler <- rep(0,nrow(Interaction))
Preschooler[which(Interaction$Age >= 3 & Interaction$Age < 8)] = 1

Preteen <- rep(0,nrow(Interaction))
Preteen[which(Interaction$Age >= 8 & Interaction$Age < 13)] = 1

Adolescence <- rep(0,nrow(Interaction))
Adolescence[which(Interaction$Age >= 13 & Interaction$Age <= 20)] = 1

Child <- rep(0,nrow(Interaction))
Child[which(Interaction$Age >= 3 & Interaction$Age <= 20)] = 1

Infant <- rep(0,nrow(Interaction))
Infant[which( Interaction$Age < 3)] = 1


Age <- data.frame(AID,Infancy,Preschooler,Preteen,Adolescence)
#Age <- data.frame(AID,Infant,Child)
#############

#############
Male <- rep(0,nrow(Interaction))
Male[which(Interaction$Gender == "Male")] = 1

Female <- rep(0,nrow(Interaction))
Female[which(Interaction$Gender == "Female")] = 1

Intracranial <- rep(NA,nrow(Interaction))
Intracranial[which(Interaction$Site == "Intracranial")] = 1
Intracranial[which(Interaction$Site != "Intracranial" & Interaction$Site != "No")] = 0

Extracranial <- rep(NA,nrow(Interaction))
Extracranial[which(Interaction$Site != "Intracranial" & Interaction$Site != "No")] = 1
Extracranial[which(Interaction$Site == "Intracranial")] = 0

Info <- data.frame(AID,Male,Female,Intracranial,Extracranial)
#############

##################################
x <- c()
for(i in c(11:39)){
  if(length(which(Interaction[,i] != "WT")) >= 4){
     x <- c(x,i)
  }
}
Mut1 <- Interaction[,x];Mut<- Mut1; Mut[,] <- 0
Mut[Mut1!="WT"] <- 1 

#Mut$KIT[which(Interaction$chr4q12=="Gain")] <- 1
##################################
x <- c()
for(i in c(40:63)){
  if(length(which(Interaction[,i] == "Gain")) >= 5){
    x <- c(x,i)
  }
}
Gain1 <- Interaction[,x];Gain<- Gain1; Gain[,] <- 0
Gain[Gain1=="Gain"] <- 1 ; colnames(Gain) <- paste0(colnames(Gain),".Gain")

##################################
x <- c()
for(i in c(40:63)){
  if(length(which(Interaction[,i] == "Loss")) >= 5){
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
cutoff <- max(summary(all)[[4]],summary(all)[[5]])
High_CNA[which(Interaction$TotalCNA > cutoff)] = 1
High_CNA[which(Interaction$TotalCNA <= cutoff & Interaction$TotalCNA >= 0)] = 0

High_TE <- rep(NA,nrow(Interaction))
all<-Interaction$TotalTE[which(Interaction$TotalTE>=0)]
cutoff <- max(summary(all)[[4]],summary(all)[[5]])
High_TE[which(Interaction$TotalTE > cutoff)] = 1
High_TE[which(Interaction$TotalTE <= cutoff & Interaction$TotalTE >= 0)] = 0

High_SV <- rep(NA,nrow(Interaction))
all<-Interaction$TotalSV[which(Interaction$TotalSV>=0)]
cutoff <- max(summary(all)[[4]],summary(all)[[5]])
High_SV[which(Interaction$TotalSV > cutoff)] = 1
High_SV[which(Interaction$TotalSV <= cutoff & Interaction$TotalSV >= 0)] = 0

High_SNV <- rep(NA,nrow(Interaction))
all<-Interaction$TotalSNV[which(Interaction$TotalSNV>=0)]
cutoff <- max(summary(all)[[4]],summary(all)[[5]])
High_SNV[which(Interaction$TotalSNV > cutoff)] = 1
High_SNV[which(Interaction$TotalSNV <= cutoff &  Interaction$TotalSNV >= 0)] = 0

SBS18_or_GCT <- rep(NA,nrow(Interaction))
SBS18_or_GCT[which(Interaction$SBS18_or_SBSgct == "Mut")] = 1
SBS18_or_GCT[which(Interaction$SBS18_or_SBSgct == "WT")] = 0

SBS39 <- rep(NA,nrow(Interaction))
SBS39[which(Interaction$SBS39 == "Mut")] = 1
SBS39[which(Interaction$SBS39 == "WT")] = 0

SBS18 <- rep(NA,nrow(Interaction))
SBS18[which(Interaction$SBS18 == "Mut")] = 1
SBS18[which(Interaction$SBS18 == "WT")] = 0

SBSgct <- rep(NA,nrow(Interaction))
SBSgct[which(Interaction$SBSgct == "Mut")] = 1
SBSgct[which(Interaction$SBSgct == "WT")] = 0

Genome <- data.frame(AID,WGD,High_SNV,High_CNA,High_TE,High_SV,SBS18,SBSgct,SBS39)


ClinicalInfo <- merge(Histology,Age); ClinicalInfo<- merge(ClinicalInfo ,Info)
GenomicAlter <- cbind(Genome,Mut,Gain,Loss)
GenomicAlter$chrX.alter <- 0
GenomicAlter$chrX.alter[which(GenomicAlter$chrXq.Gain == 1 | GenomicAlter$chrXq.Loss == 1)] <- 1
#step2: conduct chisq.test for every-pair features and summarize the pvalue results into a table

ClinicalInfo <- ClinicalInfo[,c( "AID","YST" ,"Teratomas","Germinoma","Mix","Infancy","Preschooler","Preteen","Adolescence","Extracranial","Intracranial","Male","Female")]

GenomicAlter <- GenomicAlter[,c("AID","KIT","KRAS","RRAS2","MAPK","PI3K","chr1p36.Loss","chr1q21.Gain", 
                                "chr2p24.Gain","chr3p23.Gain","chr4q12.Gain","chr4q.Loss","chr6q25.Loss","chr8p.Loss","chr8p.Gain","chr8q13.Gain","chr11q24.Loss","chr12p.Gain",
                                "chr13q.Loss","chr16p.Loss","chr18q21.Loss","chr20q.Gain", "chr21q.Gain","chr22q11.Gain","chrXq.Gain","chrXq.Loss", "WGD")]


nfeaA <- ncol(ClinicalInfo)
nfeaB <- ncol(GenomicAlter)
myTable <- data.frame(rep(0,(nfeaA-1)*(nfeaB-1)),rep(0,(nfeaA-1)*(nfeaB-1)),rep(0,(nfeaA-1)*(nfeaB-1)),rep(0,(nfeaA-1)*(nfeaB-1)) )
colnames(myTable) <- c("VarA","VarB","OddsR","pValue")
n=1
for(i in c(2:nfeaA)){
  for(j in c(2:nfeaB)){
    tmpA <- ClinicalInfo[,c(1,i)]
    tmpB <- GenomicAlter[,c(1,j)]
    tmpC <- merge(tmpA,tmpB); colnames(tmpC) <- c("ID","Var1","Var2")
    
    Mu.FEtest <- cbind(c(0,0),c(0,0))
    Mu.FEtest[1,1] <- length(which(tmpC$Var1 == 1 & tmpC$Var2 == 1))
    Mu.FEtest[1,2] <- length(which(tmpC$Var1 == 1 & tmpC$Var2 == 0))
    Mu.FEtest[2,1] <- length(which(tmpC$Var1 == 0 & tmpC$Var2 == 1))
    Mu.FEtest[2,2] <- length(which(tmpC$Var1 == 0 & tmpC$Var2 == 0))
    if(Mu.FEtest[1,1] == 0 | Mu.FEtest[1,2] == 0 | Mu.FEtest[2,1] == 0 | Mu.FEtest[2,2] == 0){
      OddsR <- ((Mu.FEtest[1,1]+0.5)*(Mu.FEtest[2,2]+0.5)) / ((Mu.FEtest[1,2]+0.5)*(Mu.FEtest[2,1]+0.5))
    }
    else{
      OddsR <- (Mu.FEtest[1,1]*Mu.FEtest[2,2]) / (Mu.FEtest[1,2]*Mu.FEtest[2,1])
    }
    #fisher exact test
    pValue <- fisher.test(Mu.FEtest,alternative ="two.sided")$p.value
    #OddsR  <- fisher.test(Mu.FEtest,alternative ="two.sided")$estimate[[1]]
    
    myTable$VarA[n] <- colnames(ClinicalInfo)[i]
    myTable$VarB[n] <- colnames(GenomicAlter)[j]
    myTable$OddsR[n] <- log2(OddsR)
    myTable$pValue[n] <- pValue
    n<-n+1
  }
}

#myTable$qValue <- p.adjust(myTable$pValue,method = "bonferroni")
myTable$qValue <- p.adjust(myTable$pValue,method = "fdr")


myTable$Dotsize <- -log10(myTable$qValue)
myTable$Dotsize[which(myTable$qValue > 0.05)] <- 0

myTable$Dotcolor <- "NS"
myTable$Dotcolor[which(myTable$OddsR > 0 & myTable$qValue <= 0.1)] <- "Cooccur" 
myTable$Dotcolor[which(myTable$OddsR < 0 & myTable$qValue <= 0.1)] <- "Exclusive" 


#step3: draw the figure
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}
plot_1.data<-myTable
plot_1.data$orderID<-c(1:nrow(plot_1.data))
plot_2.data <- plot_1.data[which(plot_1.data$Dotsize > 0),]
plot_3.data <- plot_1.data[which(plot_1.data$Dotsize == 0),]
coMu.plot<-ggplot()+theme_classic()
coMu.plot<-coMu.plot+geom_point(data = plot_1.data,aes(x=reorder(VarB,orderID),y=reorder(VarA,orderID),fill=OddsR),size=0,shape=NA)+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot<-coMu.plot+geom_vline(xintercept = c(1:26),linetype=2,color="grey",size=0.3)
coMu.plot<-coMu.plot+geom_hline(yintercept = plot_1.data$VarA,linetype=2,color="grey",size=0.3)
coMu.plot<-coMu.plot+geom_point(data = plot_2.data,aes(x=reorder(VarB,orderID),y=reorder(VarA,orderID),size=Dotsize,fill=OddsR),shape=21)+ylab(NULL)+xlab(NULL)+ggtitle(NULL)

coMu.plot<-coMu.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(1,6,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.8,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='top',
                           legend.text=element_text(size=12,face='plain'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_text(size=14,vjust=0.5,hjust=0,face='plain'),
                           axis.text.x=element_text(size=14,angle=45,vjust=0,hjust=0,face='plain',color='black'),axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='plain',color='black'),
                           axis.title.y=element_text(size=20,hjust=0.5,vjust=2,face='plain',color='black'),axis.line = element_blank(),
                           strip.text = element_text(size=18,face='bold',vjust=0.5,hjust=0.5),strip.background = element_rect(colour="black", fill=gg_color_hue(3)))
#,legend.direction="horizontal"
coMu.plot<-coMu.plot+scale_size(name   = "Fisher's test\n   q-value", range = c(2,8),transform = "exp",  breaks = c( 1,1.3,2,3,4,5,6,7),labels = expression(0.1,0.05,10^-2,10^-3, 10^-4,10^-5, 10^-6,10^-7),guide = guide_legend(nrow=2))
#coMu.plot<-coMu.plot+scale_size(name   = "Fisher's test\n   q-value", range = c(0.5,8),trans = "hms",guide = guide_legend(nrow=3))

coMu.plot<-coMu.plot+scale_fill_gradientn(name="Odds ratio\n(log2)",colours=c('#4d9221','#b8e186','white','#f1b6da','#c51b7d'),limits=c(-7,7),breaks=seq(-6,6,3))
coMu.plot<-coMu.plot+scale_x_discrete(position = "top")

plot1 <- coMu.plot

figure_1<-rbind(ggplotGrob(plot1),size="last")
ggsave(file="./Clinical_Mutation_figure3.pdf", plot=figure_1,bg = 'white', width = 22.85, height = 13.5, units = 'cm', dpi = 600)





