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


GenomicAlter <- cbind(Genome,Mut,Gain,Loss)
GenomicAlter$chrX.alter <- 0
GenomicAlter$chrX.alter[which(GenomicAlter$chrXq.Gain == 1 | GenomicAlter$chrXq.Loss == 1)] <- 1
#step2: conduct chisq.test for every-pair features and summarize the pvalue results into a table

GenomicAlter <- GenomicAlter[,c("AID","KIT","KRAS","RRAS2","NRAS","CBL","PTEN","MTOR","MAPK","PI3K","chr1p36.Loss","chr1q21.Gain", 
                                "chr2p24.Gain","chr3p23.Gain","chr4q12.Gain","chr4q.Loss","chr6q25.Loss","chr8p.Loss","chr8p.Gain","chr8q13.Gain","chr11q24.Loss","chr12p.Gain",
                                "chr13q.Loss","chr16p.Loss","chr18q21.Loss","chr20q.Gain", "chr21q.Gain","chr22q11.Gain","chrXq.Gain","chrXq.Loss", "WGD","High_SNV","High_CNA","High_TE","High_SV","SBS18","SBSgct")]

ClinicalInfo <- Interaction[,c( "AID","Seq" ,"Histology","Karyotype","Age","Gender","Site","Race")]

Final <- merge(ClinicalInfo,GenomicAlter,1,1)
DNM <- read.table("../../DNM/FunctionalDNMs.txt",header = T)


Final <- merge(Final,DNM,1,1,all.x = T,all.y = T)

Rare <- read.table("RareGermline.txt",header = T)

new <- merge(Rare,Final,all.x=T)


gene.table1<-read.table('../../Survival/Xinhua.v2.txt',header=T)
gene.table1<-gene.table1[which(!is.na(gene.table1$sx_date) & !is.na(gene.table1$last_fup_date)  ),]

gene.table1 <- gene.table1 %>% 
  mutate(
    sx_date = as.Date(sx_date, format = "%d/%m/%Y"), 
    last_fup_date = as.Date(last_fup_date, format = "%d/%m/%Y") 
  )

gene.table1 <- gene.table1 %>% 
  mutate(
    time = as.numeric( difftime(last_fup_date, sx_date, units = "days") ) / 30.4375
  )



sur <- data.frame(gene.table1$CaseID,gene.table1$Tumor,gene.table1$time,gene.table1$status)
colnames(sur) <- c("AID","TumorStage","SurvivalTime","SurvivalStatus")

new2 <- merge(new,sur,all.x=T)
write.table(new2,"Summary_all_GCT.txt",row.names = F,quote = F,sep="\t")


