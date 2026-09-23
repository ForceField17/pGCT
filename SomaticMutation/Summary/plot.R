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
library(ggbeeswarm)
library(ggpubr)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
sample <- read.table('./results/Summary.txt',header = T)
#clinical <- read.table('../../lib/Clinical.txt',header=T)


sample$xx <- sample$His
sample$xx[which(sample$His=="Teratoma")] <- "b.Teratoma"
sample$xx[which(sample$His=="IMT")] <- "a.IMT"
sample$xx[which(sample$His=="Germinoma")] <- "c.Germinoma"
sample$xx[which(sample$His=="Mix")] <- "c.Mix"

sample$TE <- sample$LINE1 + sample$SVA + sample$Alu #+ sample$HERV
sample$Indel <- sample$PointMut - sample$SNV 

sample <- sample[which(sample$SNV>=0 & sample$AID!="T135"),]

#############


F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_boxplot(data=sample,aes(x=xx,y=CNA),width=0.5,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=sample,aes(x=xx,y=CNA,color=Seq),width = 0.3,size=1.5,alpha=0.75,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,2.5,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("#CNA"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,66),breaks = seq(0,160,20)) 

F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial='#addd8e',No="grey",WGS="#9970ab"))
PartA <- F1a.plot  
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./CNA_boxplot.pdf"), plot=plotxxx,bg = 'white', width = 8, height = 7, units = 'cm', dpi = 600)

the <- compare_means(CNA ~ xx,  data = sample,p.adjust.method = "fdr")
the 

#############
F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_boxplot(data=sample,aes(x=xx,y=TotalSV),width=0.5,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=sample,aes(x=xx,y=TotalSV,color=Seq),width = 0.3,size=1.5,alpha=0.75,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,2.5,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("#SV"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,130),breaks = seq(0,120,40)) 

F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial='#addd8e',No="grey",WGS="#9970ab"))
PartB<- F1a.plot  
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./TotalSV_boxplot.pdf"), plot=plotxxx,bg = 'white', width = 8, height = 7, units = 'cm', dpi = 600)

the <- compare_means(TotalSV ~ xx,  data = sample,p.adjust.method = "fdr")
the 


#############
F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_boxplot(data=sample,aes(x=xx,y=TE),width=0.5,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=sample,aes(x=xx,y=TE,color=Seq),width = 0.3,size=1.5,alpha=0.75,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("#TE"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,49),breaks = seq(0,120,15)) 

F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial='#addd8e',No="grey",WGS="#9970ab"))
PartC <- F1a.plot  
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./TE_boxplot.pdf"), plot=plotxxx,bg = 'white', width = 8, height = 7, units = 'cm', dpi = 600)

the <- compare_means(TE ~ xx,  data = sample,p.adjust.method = "fdr")
the 



plotxxx<-rbind(ggplotGrob(PartA),ggplotGrob(PartB),ggplotGrob(PartC),size="first")
ggsave(file=paste0("./figure3f.pdf"), plot=plotxxx,bg = 'white', width = 7, height = 15, units = 'cm', dpi = 600)

