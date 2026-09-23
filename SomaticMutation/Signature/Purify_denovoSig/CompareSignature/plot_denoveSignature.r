# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library("gridExtra")
library(ggplot2)
library(sigfit)
library(foreach)
library(doParallel)
library(reshape2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}


cosmic.sigs <- read.table("../../COSMIC_v3.3.1_SBS_GRCh38.txt",header = T,row.names = "MutationType")

mutName <- read.table("../../StandardName.txt",header = T)

################################################################################
plotspectrum = function(freqs, samplename=NULL){
  if (is.null(samplename)){
    samplename = ""
  }
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  colvec = rep(c("dodgerblue","black","red","grey60","olivedrab3","#f1b6da"),each=16)
  #freqs_full = freqs[full_vec]
  freqs_full = freqs
  freqs_full[is.na(freqs_full)] = 0
  names(freqs_full) = full_vec
  y = freqs_full; maxy = max(y)
  #h=barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=0.5, cex.names=0.5, names.arg=xstr, ylab="Relative contribution")
  #axis(side=1, labels = FALSE, at=h, cex.names=0.5)
  #mtext(side=4, text=samplename)
  Table <- data.frame(full_vec,xstr,freqs_full)
  myplot <- ggplot(data=Table,aes(x=full_vec,y=freqs_full))+theme_classic()
  myplot <- myplot + geom_bar(data=Table,aes(x=full_vec,y=freqs_full),fill=colvec,alpha=1,width=0.6,stat='identity',position=position_stack())
  
  myplot <- myplot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),axis.line.y = element_blank(),
                            text=element_text(size=12,face='plain'),legend.position="none",axis.text.x=element_text(size=8,angle=90,vjust=0.5,hjust=1,face="plain",color='black'),
                            axis.text.y=element_text(size=12,face='plain',color='black'),axis.title.x=element_blank(),axis.title.y=element_text(size=12,face='plain',color='black'))
  myplot <- myplot+ggtitle(NULL)+xlab("Context")+ylab('Relative contribution')
  myplot <- myplot+scale_y_continuous(expand=c(0,0),limits = c(0,maxy*1.4),breaks = seq(0,1,0.02))+scale_x_discrete(labels=xstr)#
  
  for (j in 1:6) {
    xpos = c((j-1)*16+1,j*16)
    myplot <- myplot + geom_segment(x = xpos[1]-0.3, y = maxy*1.2, xend = xpos[2]+0.3 ,yend = maxy*1.2,color = colvec[j*16],linetype=1,size=5)
    myplot <- myplot+geom_text(x=mean(xpos), y=maxy*1.33, size=4, label=sub_vec[j])
  } 
  return(ggplotGrob(myplot))
}
################################################################################


### load in sigprofiler results
YSTsig = read.table("../../WGS/April6/SBS96/All_Solutions/SBS96_7_Signatures/Signatures/SBS96_S7_Signatures.txt", header=T )
rownames(YSTsig)<- YSTsig$MutationType
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigA = YSTsig[Order$Type1,]
figure <- plotspectrum(YSTsigA$SBS96G)
ggsave(file="./SigProfile_SBSyst.pdf", plot=figure,bg = 'white', width = 28, height = 7, units = 'cm', dpi = 600)


### load in sigfit results
YSTsig = read.table("../../WGS/SigFit/original_denovo.txt", header=T )
rownames(YSTsig)<- YSTsig$MutationType
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigB = YSTsig[Order$Type2,]
figure <- plotspectrum(YSTsigB$SignatureA)
ggsave(file="./SigFit_SBSyst.pdf", plot=figure,bg = 'white', width = 28, height = 7, units = 'cm', dpi = 600)


### load in Palimpsest results
YSTsig = read.table("../../WGS/DenoveSignaturePalimpsest.txt", header=T )
rownames(YSTsig)<- YSTsig$MutationType
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigC = YSTsig[Order$Type1,]
figure <- plotspectrum(YSTsigC$SBS96F)
ggsave(file="./Palimpsest_SBSyst.pdf", plot=figure,bg = 'white', width = 28, height = 7, units = 'cm', dpi = 600)



################################################################################
#Comparison
cosine_sim(YSTsigA$SBS96G,YSTsigB$SignatureA)  #  sigprofiler vs sigfit
cosine_sim(YSTsigA$SBS96G,YSTsigC$SBS96F)      #  sigprofiler vs Palimpsest
cosine_sim(YSTsigB$SignatureA,YSTsigC$SBS96F)  #  Palimpsest  vs sigfit



################################################################################
#Purified SBSyst
SBSyst = read.table("../SigFit/FitExt_denove.txt", header=T, row.names="MutationType" )
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigD = SBSyst[Order$Type2,]
figure <- plotspectrum(YSTsigD$SignatureB)
ggsave(file="./Purified_SBSyst_SigFit.pdf", plot=figure,bg = 'white', width = 28, height = 7, units = 'cm', dpi = 600)

SBSyst = read.table("../Palimpsest/FitExt_denove.txt", header=T, row.names="MutationType" )
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigE = SBSyst[Order$Type2,]
figure <- plotspectrum(YSTsigE$SignatureB)
ggsave(file="./Purified_SBSyst_Palimpsest.pdf", plot=figure,bg = 'white', width = 28, height = 7, units = 'cm', dpi = 600)

SBSyst = read.table("../SigProfiler/FitExt_denove.txt", header=T, row.names="MutationType" )
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigF = SBSyst[Order$Type2,]
figure <- plotspectrum(YSTsigF$SignatureD)
ggsave(file="./Purified_SBSyst_SigProfiler.pdf", plot=figure,bg = 'white', width = 28, height = 7, units = 'cm', dpi = 600)


#Comparison
cosine_sim(YSTsigD$SignatureB,YSTsigE$SignatureB)  #  Palimpsest  vs sigfit
cosine_sim(YSTsigD$SignatureB,YSTsigF$SignatureD)  #  sigprofiler vs sigfit
cosine_sim(YSTsigE$SignatureB,YSTsigF$SignatureD)      #  sigprofiler vs Palimpsest


################################## Compare with Cosmic
COSMIC <- cosmic.sigs[Order$Type1,]
SimilaritySigFit      <- rep(0,ncol(COSMIC))
SimilarityPalimpsest  <- rep(0,ncol(COSMIC))
SimilaritySigProfiler <- rep(0,ncol(COSMIC))
for(i in 1:ncol(COSMIC)){
  SimilaritySigFit[i] <- cosine_sim(YSTsigD$SignatureB,COSMIC[,i]) 
  SimilarityPalimpsest[i] <- cosine_sim(YSTsigE$SignatureB,COSMIC[,i]) 
  SimilaritySigProfiler[i] <- cosine_sim(YSTsigF$SignatureD,COSMIC[,i]) 
}

SimTable <- data.frame(colnames(COSMIC),SimilaritySigFit,SimilarityPalimpsest,SimilaritySigProfiler)

###################################### save final selected signatures
COSMIC$SBSyst <- YSTsigD$SignatureB
Order <- mutName[order(mutName$Type1),] # put in standard order
FinalSpectrum = COSMIC[Order$Type1,]
write.table(FinalSpectrum,file = "../mySpectrum.txt",sep = "\t",quote = F,row.names = T)

FinalSpectrum <- FinalSpectrum[,c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS31","SBS39","SBSyst")]
write.table(FinalSpectrum,file = "../mySpectrum_finalized.txt",sep = "\t",quote = F,row.names = T)


#SBSyst = read.table("../../WGS/SigFit/original_denovo.txt", header=T, row.names="MutationType" )
#Order <- mutName[order(mutName$Type1),] # put in standard order
#YSTsigF = SBSyst[Order$Type2,]
#rownames(YSTsigF) <- Order$Type1
#write.table(YSTsigF,file = "../../WGS/DenoveSignatureSigFit.txt",sep = "\t",quote = F,row.names = T)

