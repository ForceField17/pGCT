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
sample <- read.table('Depth_results.txt',header = T)


gene.table <- sample[which(sample$Seq=="WGS"),]

Tumor <- gene.table[which(gene.table$Tissue=="Tumor"),]
Blood <- gene.table[which(gene.table$Tissue=="Blood"),]
Parent <- gene.table[which(gene.table$Tissue=="aParents"),]

#gene.table2 <- gene.table[which((gene.table$Case %in% Tumor & gene.table$Case %in% Blood) | gene.table$Cohort=="Xinhua"),]
gene.table2 <- gene.table[which(gene.table$Case %in% c(Tumor$Case,Blood$Case,Parent$Case)),]

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_quasirandom(data=gene.table2,aes(x=Tissue,y=MedianDepth,color=Cohort,fill=Cohort),width = 0.4,size=1.5,alpha=0.5,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ geom_boxplot(data=gene.table2,aes(x=Tissue,y=MedianDepth),width=0.3,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="top",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("Median reads coverage"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks = seq(0,100,20)) 

F1a.plot <- F1a.plot + scale_fill_manual(name=NULL,values=c(Xinhua='#e31a1c',CBTTC='grey30',StJude='grey30',WANG='grey30',BASIC3='grey30'))
F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Xinhua='#e31a1c',CBTTC='grey30',StJude='grey30',WANG='grey30',BASIC3='grey30'))
F1b.plot <- F1a.plot  
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./SequencingDepthWGS.pdf"), plot=plotxxx,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)

the <- compare_means(MedianDepth ~ Tissue,  data = gene.table2)
the 

########################################3
gene.table <- sample[which(sample$Seq=="WES" & !(sample$Case %in% gene.table2$Case[which(gene.table2$Tissue=="Tumor")])),]

gene.tableA <- sample[which(sample$Seq=="WES" & sample$Tissue=="Tumor" & !(sample$Case %in% gene.table2$Case[which(gene.table2$Tissue=="Tumor")])),]
gene.tableB <- sample[which(sample$Seq=="WES" & sample$Tissue=="Blood" & !(sample$Case %in% gene.table2$Case[which(gene.table2$Tissue=="Blood")])),]


Tumor <- unique(gene.table$Case[which(gene.table$Tissue=="Tumor")])
Blood <- unique(gene.table$Case[which(gene.table$Tissue=="Blood")])

gene.table3 <- gene.table[which(gene.table$Case %in% Tumor & gene.table$Case %in% Blood),]

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_quasirandom(data=gene.table3,aes(x=Tissue,y=MedianDepth,color=Cohort,fill=Cohort),width = 0.4,size=1.5,alpha=0.5,stroke=0.5, varwidth = T,stat = "identity",dodge.width = 0)
F1a.plot<-F1a.plot+ geom_boxplot(data=gene.table3,aes(x=Tissue,y=MedianDepth),width=0.3,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="top",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("Median reads coverage"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,500),breaks = seq(0,500,100)) 
F1a.plot <- F1a.plot + scale_fill_manual(name=NULL,values=c(Xinhua='red',CBTTC='grey30',StJude='grey30',WANG='grey30',BASIC3='grey30'))
F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(Xinhua='red',CBTTC='grey30',StJude='grey30',WANG='grey30',BASIC3='grey30'))

plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./SequencingDepthWES.pdf"), plot=plotxxx,bg = 'white', width = 8, height = 10, units = 'cm', dpi = 600)

the2 <- compare_means(MedianDepth ~ Tissue,  data = gene.table3)
the2 







figure_1<-cbind(ggplotGrob(F1b.plot),ggplotGrob(F1a.plot),size="last")

panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]

figure_1$widths[7]  <- unit(3/5,'null')
figure_1$widths[20] <- unit(2/5,'null')
ggsave(file=paste0("./SequencingDepth.pdf"), plot=figure_1,bg = 'white', width = 22, height = 10, units = 'cm', dpi = 600)

nrow(gene.table2[which(gene.table2$Tissue=="aParents"),])
nrow(gene.table2[which(gene.table2$Tissue=="Blood"),])
nrow(gene.table2[which(gene.table2$Tissue=="Tumor"),])

median(gene.table2$MedianDepth[which(gene.table2$Tissue=="aParents")])
median(gene.table2$MedianDepth[which(gene.table2$Tissue=="Blood")])
median(gene.table2$MedianDepth[which(gene.table2$Tissue=="Tumor")])


nrow(gene.table3[which(gene.table3$Tissue=="Blood"),])
nrow(gene.table3[which(gene.table3$Tissue=="Tumor"),])

median(gene.table3$MedianDepth[which(gene.table3$Tissue=="Blood")])
median(gene.table3$MedianDepth[which(gene.table3$Tissue=="Tumor")])

mean(gene.table2$MedianDepth[which(gene.table2$Tissue=="aParents")])
mean(gene.table2$MedianDepth[which(gene.table2$Tissue=="Blood")])
mean(gene.table2$MedianDepth[which(gene.table2$Tissue=="Tumor")])
mean(gene.table3$MedianDepth[which(gene.table3$Tissue=="Blood")])
mean(gene.table3$MedianDepth[which(gene.table3$Tissue=="Tumor")])

