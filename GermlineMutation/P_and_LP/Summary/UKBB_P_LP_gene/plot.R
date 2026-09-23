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
table2 <- read.table('./Total_case_num.txt',header = T,sep = "\t")

table2$Fraction <- table2$Num/table2$Total*100

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=label,color=label),width=0.4,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,3.4),breaks = seq(0,3,1))+
  scale_x_discrete(labels = table2$label)
  
#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,5),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Frequency (%)")
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,GCT="#998EE0",UKBB_neoplasm="#bf812d",UKBB_healthy="#656665",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,GCT="#998EE0",UKBB_neoplasm="#bf812d",UKBB_healthy="#656665",DDD="#bebada",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
#revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,AAA="#fb8072",BBB="#78c679",CCC="#bf812d",DDD="#bebada",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
#revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,AAA="#fb8072",BBB="#78c679",CCC="#bf812d",DDD="#bebada",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)


revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="All_AF_3Gene.pdf", plot=figure_2,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)



myTable <- table2
table2 <- myTable[which(myTable$Name=="CHEK2"),]
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=label,color=label),width=0.4,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,3.4),breaks = seq(0,3,1))+
  scale_x_discrete(labels = table2$label)
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,0,2,2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Frequency (%)")
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,GCT="#998EE0",UKBB_neoplasm="#bf812d",UKBB_healthy="#656665",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,GCT="#998EE0",UKBB_neoplasm="#bf812d",UKBB_healthy="#656665",DDD="#bebada",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
p1 <- revp_plot
p1


table2 <- myTable[which(myTable$Name=="PTEN"),]
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=label,color=label),width=0.4,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,3.4),breaks = seq(0,3,1))+
  scale_x_discrete(labels = table2$label)
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,0,2,0.2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Frequency (%)")
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,GCT="#998EE0",UKBB_neoplasm="#bf812d",UKBB_healthy="#656665",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,GCT="#998EE0",UKBB_neoplasm="#bf812d",UKBB_healthy="#656665",DDD="#bebada",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
p2 <- revp_plot
p2


table2 <- myTable[which(myTable$Name=="BLM"),]
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=label,color=label),width=0.4,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,3.4),breaks = seq(0,3,1))+
  scale_x_discrete(labels = table2$label)
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,0.2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Frequency (%)")
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,GCT="#998EE0",UKBB_neoplasm="#bf812d",UKBB_healthy="#656665",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,GCT="#998EE0",UKBB_neoplasm="#bf812d",UKBB_healthy="#656665",DDD="#bebada",EEE="#ffd92f",FFF="#386cb0",GGG="#e78ac3"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
p3 <- revp_plot
p3


figure_2<-cbind(ggplotGrob(p1),ggplotGrob(p2),ggplotGrob(p3),size="first")
ggsave(file="All_AF_3Gene.pdf", plot=figure_2,bg = 'white', width = 14, height = 11, units = 'cm', dpi = 600)







########
Gene <- c("CHEK2","PTEN","BLM")
table2 <- myTable[which(myTable$Name %in% Gene),]
i<-1
zzz <- table2[which(table2$Name==Gene[i]),]
A1 <- rep(Gene[i],3)
A2 <- rep(NA,3)
C1 <- rep(NA,3)
C2 <- rep(NA,3)
C3 <- rep(NA,3)
C4 <- rep(NA,3)
A2[1] <- paste0(zzz$label[1],"-",zzz$label[2])
C1[1] <- zzz$Num[1];
C2[1] <- zzz$Total[1] - zzz$Num[1]
C3[1] <- zzz$Num[2];
C4[1] <- zzz$Total[2] - zzz$Num[2]

A2[2] <- paste0(zzz$label[1],"-",zzz$label[3])
C1[2] <- zzz$Num[1];
C2[2] <- zzz$Total[1] - zzz$Num[1]
C3[2] <- zzz$Num[3];
C4[2] <- zzz$Total[3] - zzz$Num[3]

A2[3] <- paste0(zzz$label[2],"-",zzz$label[3])
C1[3] <- zzz$Num[2];
C2[3] <- zzz$Total[2] - zzz$Num[2]
C3[3] <- zzz$Num[3];
C4[3] <- zzz$Total[3] - zzz$Num[3]
table <- data.frame(A1,A2,C1,C2,C3,C4)


for(i in c(2:3)){
zzz <- table2[which(table2$Name==Gene[i]),]
A1 <- rep(Gene[i],3)
A2 <- rep(NA,3)
C1 <- rep(NA,3)
C2 <- rep(NA,3)
C3 <- rep(NA,3)
C4 <- rep(NA,3)
A2[1] <- paste0(zzz$label[1],"-",zzz$label[2])
C1[1] <- zzz$Num[1];
C2[1] <- zzz$Total[1] - zzz$Num[1]
C3[1] <- zzz$Num[2];
C4[1] <- zzz$Total[2] - zzz$Num[2]

A2[2] <- paste0(zzz$label[1],"-",zzz$label[3])
C1[2] <- zzz$Num[1];
C2[2] <- zzz$Total[1] - zzz$Num[1]
C3[2] <- zzz$Num[3];
C4[2] <- zzz$Total[3] - zzz$Num[3]

A2[3] <- paste0(zzz$label[2],"-",zzz$label[3])
C1[3] <- zzz$Num[2];
C2[3] <- zzz$Total[2] - zzz$Num[2]
C3[3] <- zzz$Num[3];
C4[3] <- zzz$Total[3] - zzz$Num[3]
tmp <- data.frame(A1,A2,C1,C2,C3,C4)
table <- rbind(table,tmp)
}

write.table(table,file = './InputForFisherTest.txt',row.names = F,quote = F,sep = "\t")

ID <- rep(NA,nrow(table))
OR <- rep(NA,nrow(table))
pValue <- rep(NA,nrow(table))
for(i in c(1:nrow(table))){
  ID[i] <- paste0(table$A1[i],"---",table$A2[i])
  tmp <- data.frame(c(table$C1[i],table$C2[i]),c(table$C3[i],table$C4[i]))
  pValue[i] <- fisher.test(tmp,alternative = "two.sided")$p.value
  OR[i] <- fisher.test(tmp)$estimate
}

FinalTest <- data.frame(ID,OR,pValue)
FinalTest$FDR <- p.adjust(FinalTest$pValue,method = "fdr")
FinalTest

write.table(FinalTest,file = "Fisher.txt",row.names = F,quote = F,sep = "\t")









