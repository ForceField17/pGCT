# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library("gridExtra")
library(ggplot2)
library("reshape2")
library(ggbeeswarm)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l =80, c = 90)[1:n]
}


data_summary <- function(x) {
  m <- mean(x)
  ymin <- max(m-sd(x),0)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

Clinical <- read.table("../../lib/Clinical.txt",sep = "\t",header = T)
PointMut <- read.table("../Landscape/results/Final_PointMut.txt",header=T)
RAS <- data.frame(PointMut$AID,rep("WT",nrow(PointMut)))
colnames(RAS) <- c("AID","RASmut")
RAS$RASmut[which(PointMut$RRAS2 != "WT")] <- "RRAS2"
RAS$RASmut[which(PointMut$KRAS != "WT")] <- "KRAS"
RAS$RASmut[which(PointMut$NRAS != "WT")] <- "NRAS"

qq <- merge(RAS,Clinical)
qq <- qq[which(qq$RASmut!="WT"),]

NAN_plot <- ggplot(data=qq,aes(c)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(data=qq,aes(x=RASmut,y=Age),width=0.3,size=0.5,alpha=0,outlier.shape = NA)+
  geom_quasirandom(data=qq,aes(x=RASmut,y=Age,color=Histology,fill=Histology),shape=21,width = 0.4,size=2,alpha=0.7,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
NAN_plot <- NAN_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot+ggtitle(NULL)+ylab("Onset age of patient")+xlab(NULL)
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(5,21),breaks = seq(0,200,4) )
NAN_plot <- NAN_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
NAN_plot <- NAN_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))

plotxxx<-cbind(ggplotGrob(NAN_plot ),size="first")
ggsave(file=paste0("./Age_RAS_mut.pdf"), plot=plotxxx,bg = 'white', width = 7, height = 9, units = 'cm', dpi = 600)

the <- compare_means(Age ~ RASmut,data = qq,p.adjust.method ="BH" )
the 


qq$xx <- qq$RASmut
qq$xx[which(qq$RASmut!="KRAS")] <- "NRAS/RRAS2"
NAN_plot <- ggplot(data=qq,aes(c)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(data=qq,aes(x=xx,y=Age),width=0.3,size=0.5,alpha=0,outlier.shape = NA)+
  geom_quasirandom(method="quasirandom",data=qq,aes(x=xx,y=Age,color=Histology,fill=Histology),shape=21,width = 0.4,size=2,alpha=0.7,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
NAN_plot <- NAN_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot+ggtitle(NULL)+ylab("Onset age of patient")+xlab(NULL)
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(5,21),breaks = seq(0,200,4) )
NAN_plot <- NAN_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
NAN_plot <- NAN_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))

plotxxx<-cbind(ggplotGrob(NAN_plot ),size="first")
ggsave(file=paste0("./Age_RAS_mut_2.pdf"), plot=plotxxx,bg = 'white', width = 7, height = 9, units = 'cm', dpi = 600)

the <- compare_means(Age ~ xx,data = qq,p.adjust.method ="BH" )
the 
