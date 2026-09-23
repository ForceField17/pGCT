# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(DESeq2)
library(stringr)
library("pheatmap")
library("RColorBrewer")
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(preprocessCore)
library(clusterProfiler)
library(enrichplot)
library(ggbeeswarm)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#raw data preprocessing
Sample_features <- read.table("../../../TGCT_subtyping/Info_temp.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- gsub("-",".",samples$AID)
samples$AID <- rownames(samples) 


##########################################
myTable <- samples
table1 <- myTable[which(!is.na(myTable$RT) & !is.na(myTable$chr1p36) & !is.na(myTable$methy_L1HSpro) ),]

table1$His <- "aOthers"
table1$His[which(table1$Histology == "YST")] <- "YST"
table1$His[which(table1$Histology == "Germinoma")] <- "Germinoma"
#table1$Tag <- paste0(table1$Tag,".",table1$chr1p36)

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_histogram(data=table1,aes(x=methy_L1HSpro,y=..density..),fill="#c69c72",binwidth=0.025,alpha=0.3,position='identity')
F1a.plot<-F1a.plot+ geom_density(data=table1 ,aes(x=methy_L1HSpro),color="#c69c72",size=1,alpha=1)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="bottom",legend.text=element_text(size=10,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab("L1HSpro methylation")+ylab('Density')
F1a.plot<-F1a.plot+geom_vline( aes(xintercept=median(table1$methy_L1HSpro)),colour="#f25c54",linetype=2, size=0.5) + 
  scale_color_manual(name=NULL,values=c("#f25c54","#01baef"))
F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits =c(0,9))+scale_x_continuous(expand=c(0,0),limits=c(0.5,1),breaks=seq(-1,1,0.1))#
F1a.plot
figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="distri_methy.pdf", plot=figure_1,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)
cor.test(table1$methy_L1HSpro ,table1$RT,method="spearman")

median(table1$methy_L1HSpro)
summary(table1$methy_L1HSpro)




#########

table1$group <- "Methy_high"
table1$group[which(table1$methy_L1HSpro < summary(table1$methy_L1HSpro)[3] )] <- "aMethy_low"
table(table1$group)
table1$Tag <- paste0(table1$group," + ",table1$chr1p36)
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=His),shape=21,width = 0.16,size=2,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits = c(-0.5,23),breaks = c(0,5,10,15,20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =c(Germinoma="#FFA097",YST="#FF99FF",aOthers="grey"))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("L1") +xlab(NULL)

the <- compare_means(RT ~ Tag,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the
figure_2<-cbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./Tag2_L1_methy_groups.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)
####################################


table(table1$Tag)
res.aov2 <- aov(RT ~  group * chr1p36, data = table1)
summary(res.aov2)
eta_squared(res.aov2 , partial = T)





#######################################33
F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_histogram(data=table1,aes(x=exp_L1HS,y=..density..),fill="#c69c72",binwidth=0.5,alpha=0.3,position='identity')
F1a.plot<-F1a.plot+ geom_density(data=table1 ,aes(x=exp_L1HS),color="#c69c72",size=1,alpha=1)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="bottom",legend.text=element_text(size=10,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab("L1HSpro methylation")+ylab('Density')
F1a.plot<-F1a.plot+geom_vline( aes(xintercept=median(table1$exp_L1HS)),colour="#f25c54",linetype=2, size=0.5) + 
  scale_color_manual(name=NULL,values=c("#f25c54","#01baef"))
F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits =c(0,0.5))+scale_x_continuous(expand=c(0,0.5))#
F1a.plot
figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="distri_exp_L1HS.pdf", plot=figure_1,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)
cor.test(table1$exp_L1HS ,table1$RT,method="spearman")

median(table1$exp_L1HS)
summary(table1$exp_L1HS)




#########

table1$group <- "Exp_low"
table1$group[which(table1$exp_L1HS > summary(table1$exp_L1HS)[3] )] <- "aExp_high"
table(table1$group)
table1$Tag <- paste0(table1$group," + ",table1$chr1p36)
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=His),shape=21,width = 0.2,size=2,alpha=0.75,stroke=0.5, varwidth = F,dodge.width = 0.1)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits = c(-0.5,23),breaks = c(0,5,10,15,20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =c(Germinoma="#FFA097",YST="#FF99FF",aOthers="grey"))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("L1") +xlab(NULL)

the <- compare_means(RT ~ Tag,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the
figure_2<-cbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./Tag3_exp_L1HS_groups.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)
####################################


table(table1$Tag)
res.aov2 <- aov(RT ~  group * chr1p36, data = table1)
summary(res.aov2)
eta_squared(res.aov2 , partial = T)