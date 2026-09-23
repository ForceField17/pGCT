# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
library(ggplot2)
library(gridExtra)
library(grid)
library("dplyr") 
library("plyr")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
table1 <- read.table('./Final_Detail.txt',header = T,sep = "\t")
table1$UKBBhealth <- table1$UKBB - table1$UKBBneoplasm
########
F1 <- rep(NA,nrow(table1))
F2 <- rep(NA,nrow(table1))
F3 <- rep(NA,nrow(table1))
P1 <- rep(NA,nrow(table1))
P2 <- rep(NA,nrow(table1))
P3 <- rep(NA,nrow(table1))

for(i in c(1:nrow(table1))){
  F1[i] <- table1$GCT[i] / 229
  F2[i] <- table1$UKBBhealth[i] / 350703
  F3[i] <- table1$UKBBneoplasm[i] / 151811
  tmp1 <- data.frame(c(table1$GCT[i],229-table1$GCT[i]),c(table1$UKBBhealth[i],350703-table1$UKBBhealth[i]))
  tmp2 <- data.frame(c(table1$GCT[i],229-table1$GCT[i]),c(table1$UKBBneoplasm[i],151811-table1$UKBBneoplasm[i]))
  tmp3 <- data.frame(c(table1$UKBBneoplasm[i],151811-table1$UKBBneoplasm[i]),c(table1$UKBBhealth[i],350703-table1$UKBBhealth[i]))
  P1[i] <- fisher.test(tmp1)$p.value
  P2[i] <- fisher.test(tmp2)$p.value
  P3[i] <- fisher.test(tmp3)$p.value
}

types <- unique(table1$AID)
FDR1 <- p.adjust(P1,method = "fdr")   #,n = length(types))
FDR2 <- p.adjust(P2,method = "fdr")   #,n = length(types))
FDR3 <- p.adjust(P3,method = "fdr")   #,n = length(types))

#theResults <- data.frame(table1$AID,table1$CaseID,table1$GCT,F1,table1$UKBBhealth,F2,table1$UKBBneoplasm,F3,)
table1$FreqGCT <- F1
table1$FreqUKBBhealth <- F2
table1$FreqUKBBneoplasm <- F3
table1$P_GCT_vs_UKBBhealth <- P1
#table1$FDR_GCT_vs_UKBBhealth <- FDR1
table1$P_GCT_vs_UKBBneoplasm <- P2
#table1$FDR_GCT_vs_UKBBneoplasm <- FDR2
table1$P_UKBBneoplasm_vs_UKBBhealth <- P3
#table1$FDR_UKBBneoplasm_vs_UKBBhealth <- FDR3
table1$Diff_GCT_vs_UKBBhealth <- table1$FreqGCT - table1$FreqUKBBhealth
table1$Diff_GCT_vs_UKBBneoplasm <- table1$FreqGCT - table1$FreqUKBBneoplasm
table1$Diff_UKBBneoplasm_vs_UKBBhealth <- table1$FreqUKBBneoplasm - table1$FreqUKBBhealth
  
write.table(table1,file = './Results_FisherTest.txt',row.names = F,quote = F,sep = "\t")







