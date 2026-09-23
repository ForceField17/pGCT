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
library(ggpubr)
library(ggbeeswarm)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
sample <- read.table('all_eroll_patients_2025.txt',header = T)

ID<-as.character(sample[,1])
case <- unique(sample$ID)


##########################################################################################################################

sample$His    <- sample$Histology
sample$His[which(sample$Histology=="NGGCT")] <- "*Notspecified"
sample$His[which(sample$Histology=="GCT")] <- "*Notspecified"


#######

#######




##################
gene.table <- sample[which(sample$Age != "No"),]
gene.table$age <-  as.numeric(gene.table$Age)
gene.table <- gene.table[which(gene.table$age<21),]


#gene.table <- gene.table[which(gene.table$Histology!="EmbryonalC" & gene.table$Histology!="GCT"  & gene.table$Histology!="NGGCT" & gene.table$Histology!="ChorioC" & gene.table$Cohort!="Japan"),]
#gene.table2 <- gene.table[which(gene.table$PrimarySite=="Intracranial" & gene.table$Histology!="GCT" ),]
gene.table2 <- gene.table[which(gene.table$PrimarySite=="Intracranial"),]
gene.table2$x <- round(gene.table2$age)
gene.table2$y <- 1

use_pathway <- group_by(gene.table2, x) %>%
  group_by(Histology) %>%
  ungroup() %>%
  mutate(Histology = factor(Histology,levels=c("YST", "IMT", "Teratoma","Germinoma","Mix"))) %>%
  dplyr::arrange( x,Histology) 

y <- c()
for(i in unique(use_pathway$x)){
  tmp <- c(1:length(which(use_pathway$x==i)))
  y <- c(y,tmp)
}
use_pathway$y <- y-0.4

gene.table2 <- use_pathway
F1a.plot<-ggplot(data=gene.table2,aes(x=x,y=y))+theme_classic()
F1a.plot<-F1a.plot+ geom_point(data=gene.table2,aes(x=x,y=y,fill=His),shape=21,color="grey20",size=5,alpha=0.6)
F1a.plot<-F1a.plot+geom_density(linetype=1,aes(y = after_stat(density*90),color=His),size=0.7,alpha=0.9,type=1)
#F1a.plot<-F1a.plot+ geom_density(data=gene.table2 ,aes(x=x,weight=100,color=His),size=0.75,alpha=0.9)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="right",legend.text=element_text(size=10,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab("Age of disease diagnosis")+ylab('#patients')+scale_fill_manual(name=NULL,values=c(gg_color_hue(5) ))
cdat <- ddply(gene.table2, "His", summarise, AF.mean=median(age))
cdat
#F1a.plot<-F1a.plot+geom_vline(data=cdat, aes(xintercept= AF.mean,colour=His),linetype=2, size=0.5) 
F1a.plot<-F1a.plot+scale_fill_manual(name=NULL,values=c("grey30",gg_color_hue(5)[1],gg_color_hue(5)[2],gg_color_hue(5)[3],"#1aa3d9",gg_color_hue(5)[5]))
F1a.plot<-F1a.plot+scale_color_manual(name=NULL,values=c("grey30",gg_color_hue(5)[1],gg_color_hue(5)[2],gg_color_hue(5)[3],"#1aa3d9",gg_color_hue(5)[5]))

F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(0,13.1),breaks = seq(0,12,2))+scale_x_continuous(expand=c(0,0),limits=c(-0.5,20.5),breaks=seq(0,20,5))#

F1a.plot
figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="./Fig_Age/test_Age_vs_Histology_Intracranial.pdf", plot=figure_1,bg = 'white', width = 18, height = 9, units = 'cm', dpi = 600)

the2 <- compare_means(age ~ Gender, data = gene.table2)
the2

the2 <- compare_means(age ~ His, data = gene.table2, method = "wilcox")
the2


#testin matrix




##################
gene.table <- sample[which(sample$Age != "No"),]
gene.table$age <-  as.numeric(gene.table$Age)
gene.table <- gene.table[which(gene.table$age<21),]
gene.table <- gene.table[which(gene.table$Histology!="EmbryonalC" & gene.table$Histology!="GCT"  & gene.table$Histology!="NGGCT" & gene.table$Histology!="ChorioC" & gene.table$Cohort!="Japan"),]
gene.table2 <- gene.table[which(gene.table$PrimarySite=="Extracranial"),]

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_histogram(data=gene.table2,aes(x=age,y=(..density..)/3,fill=His),binwidth=1,color="grey20",size=0.2,alpha=0.5,position='stack')
F1a.plot<-F1a.plot+ geom_density(data=gene.table2 ,aes(x=age,color=His),size=1,alpha=1)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="right",legend.text=element_text(size=10,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab("Age of disease diagnosis")+ylab('Density')+scale_fill_manual(name=NULL,values=c(gg_color_hue(5) ))
cdat <- ddply(gene.table2, "His", summarise, AF.mean=median(age))
cdat
F1a.plot<-F1a.plot+geom_vline(data=cdat, aes(xintercept= AF.mean,colour=His),linetype=2, size=0.5) 
F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(0,0.25))+scale_x_continuous(expand=c(0,0),limits=c(-0.5,20.5),breaks=seq(0,20,5))#
F1a.plot<-F1a.plot+scale_fill_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[3],gg_color_hue(5)[4],gg_color_hue(5)[5]))
F1a.plot<-F1a.plot+scale_color_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[3],gg_color_hue(5)[4],gg_color_hue(5)[5]))

F1a.plot
figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="./Fig_Age/Age_vs_Histology_Extracranial.pdf", plot=figure_1,bg = 'white', width = 14, height = 8, units = 'cm', dpi = 600)

the2 <- compare_means(age ~ Gender, data = gene.table2)
the2

the2 <- compare_means(age ~ His, data = gene.table2)
the2


library(ggpubr)


##################
gene.table <- sample[which(sample$Age != "No" & sample$PrimarySite!="No" & sample$His!="IMT" & sample$His!="Teratoma"),]
gene.table$age <-  as.numeric(gene.table$Age)
gene.table <- gene.table[which(gene.table$age<21),]
gene.table2 <- gene.table[which( gene.table$Histology!="GCT"  & gene.table$Histology!="NGGCT"  & gene.table$Cohort!="Japan"),]
gene.table2$class <- paste0(gene.table2$His,gene.table2$PrimarySite)
gene.table2$Order <- 1
gene.table2$Order[which(gene.table2$class=="GerminomaIntracranial")] <- 2
gene.table2$Order[which(gene.table2$class=="YSTExtracranial")] <- 3
gene.table2$Order[which(gene.table2$class=="YSTIntracranial")] <- 4
gene.table2$Order[which(gene.table2$class=="MixExtracranial")] <- 5
gene.table2$Order[which(gene.table2$class=="MixIntracranial")] <- 6

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_boxplot(data=gene.table2,aes(x=reorder(class,Order),y=age),width=0.3,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=gene.table2,aes(x=reorder(class,Order),y=age,color=His,fill=His),width = 0.25,size=1.7,alpha=0.5,stroke=0.8, varwidth = T)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("Age"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,21),breaks = seq(0,21,5)) 
the <- compare_means(age ~ class,  data = gene.table2)
the 

F1a.plot <- F1a.plot + scale_fill_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[3],gg_color_hue(5)[5]))
F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[3],gg_color_hue(5)[5]))
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./Fig_Age/All_His_Site.pdf"), plot=plotxxx,bg = 'white', width = 11, height = 10, units = 'cm', dpi = 600)

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_quasirandom(data=gene.table2,aes(x=reorder(class,Order),y=age,color=PrimarySite,fill=PrimarySite),width = 0.25,size=1.7,alpha=0.75,stroke=0.8, varwidth = T)
F1a.plot<-F1a.plot+ geom_boxplot(data=gene.table2,aes(x=reorder(class,Order),y=age),width=0.3,size=0.5,alpha=0,outlier.shape = NA)

F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("Age"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,21),breaks = seq(0,21,5)) 
the <- compare_means(age ~ class,  data = gene.table2)
the 

F1a.plot <- F1a.plot + scale_fill_manual(name=NULL,values=c("#a6dba0","#9970ab"))
F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c("#a6dba0","#9970ab"))
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./Fig_Age/All_His_Site2.pdf"), plot=plotxxx,bg = 'white', width = 11, height = 10, units = 'cm', dpi = 600)

the <- compare_means(age ~ Histology,  data = gene.table2)
the 
