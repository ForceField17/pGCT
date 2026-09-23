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
F1a.plot<-F1a.plot+ geom_boxplot(data=sample,aes(x=xx,y=LINE1),width=0.5,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=sample,aes(x=xx,y=LINE1,color=Site,fill=Site),width = 0.3,size=1.5,alpha=0.75,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("#LINE1"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,49),breaks = seq(0,45,15)) 

F1a.plot <- F1a.plot + scale_fill_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial="#9970ab",No="grey"))
F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial="#9970ab",No="grey"))
PartA <- F1a.plot  
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./LINE1_boxplot.pdf"), plot=plotxxx,bg = 'white', width = 6, height = 7, units = 'cm', dpi = 600)




#############
#############
F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_boxplot(data=sample,aes(x=xx,y=HERV),width=0.5,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=sample,aes(x=xx,y=HERV,color=Site,fill=Site),width = 0.3,size=1.5,alpha=0.75,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(0,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("#HERV"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,49),breaks = seq(0,45,15)) 

F1a.plot <- F1a.plot + scale_fill_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial="#9970ab",No="grey"))
F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial="#9970ab",No="grey"))
PartB <- F1a.plot  
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./HERV_boxplot.pdf"), plot=plotxxx,bg = 'white', width = 6, height = 6.5, units = 'cm', dpi = 600)

#############
F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_boxplot(data=sample,aes(x=xx,y=Alu),width=0.5,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=sample,aes(x=xx,y=Alu,color=Site,fill=Site),width = 0.3,size=1.5,alpha=0.75,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(0,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("#Alu"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,49),breaks = seq(0,45,15)) 

F1a.plot <- F1a.plot + scale_fill_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial="#9970ab",No="grey"))
F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial="#9970ab",No="grey"))
PartB2 <- F1a.plot  
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./Alu_boxplot.pdf"), plot=plotxxx,bg = 'white', width = 6, height = 6.5, units = 'cm', dpi = 600)



#############
F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_boxplot(data=sample,aes(x=xx,y=SVA),width=0.5,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=sample,aes(x=xx,y=SVA,color=Seq),width = 0.3,size=1.5,alpha=0.75,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(0,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("#SVA"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,49),breaks = seq(0,45,15)) 

F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Intracranial="#9970ab",Extracranial="#9970ab",No="grey",WGS="#9970ab"))
PartC <- F1a.plot  
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./SVA_boxplot.pdf"), plot=plotxxx,bg = 'white', width = 8, height = 7, units = 'cm', dpi = 600)



plotxxx<-rbind(ggplotGrob(PartA),ggplotGrob(PartB2),ggplotGrob(PartC),size="first")
ggsave(file=paste0("./figure_RT.pdf"), plot=plotxxx,bg = 'white', width = 7, height = 11, units = 'cm', dpi = 600)

the <- compare_means(HERV ~ xx,  data = sample,p.adjust.method = "fdr")
the 
the <- compare_means(LINE1 ~ xx,  data = sample,p.adjust.method = "fdr")
the 
the <- compare_means(Alu ~ xx,  data = sample,p.adjust.method = "fdr")
the 
the <- compare_means(SVA ~ xx,  data = sample,p.adjust.method = "fdr")
the 

