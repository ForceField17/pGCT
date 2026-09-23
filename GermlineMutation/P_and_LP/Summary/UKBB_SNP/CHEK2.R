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
table2 <- read.table('SNP1_CHEK2.txt',header = T,sep = "\t")

table2$Fraction <- table2$Num/table2$Total*100
table2 <- table2[which(table2$Mut=="AAA" & table2$Name!="XXX"),]

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=Mut,color=Mut),width=0.25,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.42),breaks = seq(0,2,0.4))+
  scale_x_discrete(labels = table2$label)
  
#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,5),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("% cases harboring\n22:28695868 AG>A")
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,AAA="#bebada",BBB="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,AAA="#bebada",BBB="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)

revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="CHEK2_22_28695868_AG_A.pdf", plot=figure_2,bg = 'white', width = 16, height = 12, units = 'cm', dpi = 600)


#########
table2 <- read.table('SNP1_CHEK2.txt',header = T,sep = "\t")

table2$Fraction <- table2$Num/table2$Total*100
table2 <- table2[which(table2$Mut=="BBB" & table2$Name!="XXX"),]

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=Mut,color=Mut),width=0.25,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.48),breaks = seq(0,2,0.1))+
  scale_x_discrete(labels = table2$label)

#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,5),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("% cases harboring\n22:28734532 C>T")
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,AAA="#bebada",BBB="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,AAA="#bebada",BBB="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)

revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="CHEK2_22_28734532_C_T.pdf", plot=figure_2,bg = 'white', width = 16, height = 12, units = 'cm', dpi = 600)


#########
table2 <- read.table('SNP1_CHEK2.txt',header = T,sep = "\t")

table2$Fraction <- table2$Num/table2$Total*100
table2 <- table2[which(table2$Mut=="CCC" & table2$Name!="XXX"),]

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=Mut,color=Mut),width=0.25,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.48),breaks = seq(0,2,0.1))+
  scale_x_discrete(labels = table2$label)

#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,5),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("% cases harboring\n22:28695219 G>A")
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,AAA="#bebada",BBB="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,AAA="#bebada",BBB="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)

revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="CHEK2_22_28695219_G_A.pdf", plot=figure_2,bg = 'white', width = 16, height = 12, units = 'cm', dpi = 600)






#########
#########
table2 <- read.table('SNP1_CHEK2.txt',header = T,sep = "\t")

table2$Fraction <- table2$Num/table2$Total*100
table2 <- table2[which(table2$Mut=="EEE" & table2$Name!="XXX"),]
table2
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=Mut,color=Mut),width=0.25,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.48),breaks = seq(0,2,0.1))+
  scale_x_discrete(labels = table2$label)

#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,5),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("% cases harboring\nATM mut")
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,EEE="#bebada",BBB="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,EEE="#bebada",BBB="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)

revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="ATM.pdf", plot=figure_2,bg = 'white', width = 16, height = 12, units = 'cm', dpi = 600)

#########

#########
table2 <- read.table('SNP1_CHEK2.txt',header = T,sep = "\t")

table2$Fraction <- table2$Num/table2$Total*100
table2 <- table2[which(table2$Mut=="GGG" & table2$Name!="XXX"),]
table2
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=Mut,color=Mut),width=0.25,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.96),breaks = seq(0,2,0.3))+
  scale_x_discrete(labels = table2$label)

#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,5),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("% cases harboring\n15:90790746 AC>A")
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=NA,GGG="#bebada",DDD="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=NA,GGG="#bebada",DDD="#78c679",CCC="#fc9272",zFragile="#F8766D"),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)

revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="BLM_15_90790746_AC_A.pdf", plot=figure_2,bg = 'white', width = 16, height = 12, units = 'cm', dpi = 600)

