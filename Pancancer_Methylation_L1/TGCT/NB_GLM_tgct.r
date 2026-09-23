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
Sample_features <- read.table("../../TGCT_subtyping/Info_temp.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- gsub("-",".",samples$AID)
samples$AID <- rownames(samples) 


##########################################
myTable <- samples
table1 <- myTable[which(!is.na(myTable$RT) & !is.na(myTable$chr1p36) ),]
table1$His <- table1$Histology
table1$His[which(is.na(table1$Histology) ) ] <- table1$Original_report[which(is.na(table1$Histology) ) ]

table1$Tag <- "nonYST"
table1$Tag[which(table1$His == "YST")] <- "aYST"
#table1$Tag <- paste0(table1$Tag,".",table1$chr1p36)
table(table1$Tag)
table1$Histology <- table1$Tag
res.aov2 <- aov(RT ~  Histology * chr1p36, data = table1)
summary(res.aov2)
eta_squared(res.aov2 , partial = T)

library(MASS)
m_add <- glm.nb(RT ~ Histology + chr1p36, data = table1)
m_int <- glm.nb(RT ~ Histology * chr1p36, data = table1)
anova(m_add, m_int)          # 交互项的似然比检验
exp(coef(m_int)); exp(confint(m_int))

library(ARTool)
anova(art(RT ~ as.factor(Histology) * as.factor(chr1p36), data = table1))

#"YST","YST_domin","Yes","ReportYST"
table1$his <- "nonYST"
table1$his[which(table1$YST_1st_his %in%  c("YST") )] <- "aYST"
table1$Tag <- paste0(table1$his," + ",table1$chr1p36)
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=chr1p36),shape=21,width = 0.16,size=2,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0.8)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits = c(-0.5,23),breaks = c(0,5,10,15,20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =c("#2166ac","grey90"))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("L1") +xlab(NULL)

the <- compare_means(RT ~ Tag,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the
figure_2<-cbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./figs_GLM/Tag_L1.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)
####################################
table(table1$Tag)
anova(art(RT ~ as.factor(his) * as.factor(chr1p36), data = table1))


#########
myTable <- samples
table1 <- myTable[which(!is.na(myTable$RT) & !is.na(myTable$chr1p36) & !(is.na(myTable$methy_L1HSpro)) ),] 
nrow(table1)
summary(table1$methy_L1HSpro)
table1$his <- "Methy"
table1$his[which(table1$methy_L1HSpro < summary(table1$methy_L1HSpro)[3] )] <- "hypoMethy"
table(table1$his)
table1$Tag <- paste0(table1$his," + ",table1$chr1p36)
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=chr1p36),shape=21,width = 0.16,size=2,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0.8)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits = c(-0.5,23),breaks = c(0,5,10,15,20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =c("#2166ac","grey90"))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("L1") +xlab(NULL)

the <- compare_means(RT ~ Tag,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the
figure_2<-cbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./figs_GLM/Tag2_L1.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)
####################################


#########
myTable <- samples
table1 <- myTable[which(!is.na(myTable$RT) & !is.na(myTable$chr1p36) & !(is.na(myTable$exp_L1HS)) ),] 
nrow(table1)
summary(table1$exp_L1HS)
table1$his <- "noOE_L1HS"
table1$his[which(table1$exp_L1HS > summary(table1$exp_L1HS)[3] )] <- "aOE_L1HS"
table(table1$his)
table1$Tag <- paste0(table1$his," + ",table1$chr1p36)
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=chr1p36),shape=21,width = 0.16,size=2,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0.8)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits = c(-0.5,23),breaks = c(0,5,10,15,20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =c("#2166ac","grey90"))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("L1") +xlab(NULL)

the <- compare_means(RT ~ Tag,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the
figure_2<-cbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./figs_GLM/Tag3_L1.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)
####################################


###############3
table1 <- samples 
table1$group <- "aOther"
table1$group[which(table1$YST_1st_his %in%  c("YST") )] <- "YST"
table1$group[which(table1$Histology %in%  c("Germinoma") )] <- "Germinoma"
myTable <- table1[which(!(is.na(table1$exp_L1HS))),]
myTable$rank <- 1
myTable$rank[which(myTable$group=="YST")] <- 3
myTable$rank[which(myTable$group=="Germinoma")] <- 2
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(data=myTable ,aes(y=exp_L1HS ,x=reorder(group,rank)),width=0.8,size=0.5,alpha=0,outlier.shape = NA)
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable ,aes(y=exp_L1HS ,x=reorder(group,rank),fill=group),shape=21,width = 0.3,size=1.6,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0)
#geom_beeswarm(dodge.width=0.8,cex=2,data=table5,aes(y=Freq ,x=Type,fill=chr1p36),size=2.6,shape=21,alpha=0.7,stroke=0.5, stat = "identity")
NAN_plot <- NAN_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                            axis.text.x=element_text(size=11,angle=45,vjust = 1,hjust = 1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot+ggtitle(NULL)+ylab("L1HS expression")+xlab(NULL)
NAN_plot <- NAN_plot +  scale_y_continuous(expand=c(0,0),limits=c(9.6,15.4),breaks = seq(0,100,1))
NAN_plot <- NAN_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA097",YST="#FF99FF"))
plotxxx<-cbind(ggplotGrob(NAN_plot ),size="first")
ggsave(file=paste0("figs_GLM/exp_L1HS.pdf"), plot=plotxxx,bg = 'white', width = 6, height = 9, units = 'cm', dpi = 600)
the <- compare_means(exp_L1HS ~ group,  data = myTable,method =  "wilcox")
the 

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(data=myTable ,aes(y=methy_L1HSpro ,x=reorder(group,rank)),width=0.8,size=0.5,alpha=0,outlier.shape = NA)
NAN_plot <- NAN_plot + geom_quasirandom(data=myTable ,aes(y=methy_L1HSpro ,x=reorder(group,rank),fill=group),shape=21,width = 0.3,size=1.6,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0)
#geom_beeswarm(dodge.width=0.8,cex=2,data=table5,aes(y=Freq ,x=Type,fill=chr1p36),size=2.6,shape=21,alpha=0.7,stroke=0.5, stat = "identity")
NAN_plot <- NAN_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                            axis.text.x=element_text(size=11,angle=45,vjust = 1,hjust = 1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot+ggtitle(NULL)+ylab("L1HS expression")+xlab(NULL)
NAN_plot <- NAN_plot +  scale_y_continuous(expand=c(0,0),limits=c(0.5,0.9),breaks = seq(0,100,0.1))
NAN_plot <- NAN_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA097",YST="#FF99FF"))
plotxxx<-cbind(ggplotGrob(NAN_plot ),size="first")
ggsave(file=paste0("figs_GLM/methy_L1HS.pdf"), plot=plotxxx,bg = 'white', width = 6, height = 9, units = 'cm', dpi = 600)
the <- compare_means(methy_L1HSpro ~ group,  data = myTable,method =  "wilcox")
the 







NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=His),shape=21,width = 0.26,size=2,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits = c(-0.5,23),breaks = c(0,5,10,15,20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ Tag,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/RT_Histology.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)


##########################################
myTable <- samples
table1 <- myTable[which(!is.na(myTable$RT) & !is.na(myTable$chr1p36) ),]
table1$His <- table1$Histology
table1$His[which(is.na(table1$Histology) ) ] <- table1$Original_report[which(is.na(table1$Histology) ) ]

table1$Tag <- table1$chr1p36
table(table1$chr1p36)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=His),shape=21,width = 0.26,size=2,alpha=0.75,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits = c(-0.5,23),breaks = c(0,5,10,15,20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ Tag,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/RT_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)




gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################
TCGA <- read.table("../../Pancancer_L1/Final_CNV_RT.txt",header = T)
TCGA$RT <- TCGA$RT_burden
check <- TCGA[which(TCGA$chr1p36=="Loss"),]
temp <- as.data.frame(table(check$subtype))
TCGA <- TCGA[which(TCGA$subtype %in% temp$Var1[which(temp$Freq>2)]),]
cancer <- unique(TCGA$subtype)
p <- rep(NA,length(cancer))
rho <- rep(NA,length(cancer))
diff <- rep(NA,length(cancer))
for(i in c(1:length(cancer))){
  xx <- TCGA[which(TCGA$subtype==cancer[i]),]
  #tmp <- cor.test(formula = ~ active_loci_TPM + RT_burden, data = xx,method = "spearman")
  loss <- xx$RT[which(xx$chr1p36=="Loss")]
  WT <- xx$RT[which(xx$chr1p36=="WT")]
  p[i] <- wilcox.test(loss,WT)$p.value
  diff[i] <- mean(loss) - mean(WT)
}
final <- data.frame(cancer,p,diff)
final$fdr <- -log10(p.adjust(final$p,method = "fdr"))

name <- final[which(final$fdr> -log10(0.05)),]
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_point(data = final , aes(x = diff, y = fdr),fill="grey90",color="grey10",alpha=0.7,size=2, shape = 21,stroke=0.5) 
NAN_plot <- NAN_plot + geom_hline(aes(yintercept=-log10(0.05)),color='black',size=.5,linetype='dashed') 
NAN_plot <- NAN_plot + geom_vline(aes(xintercept=0),color='black',size=.5,linetype=2)
#NAN_plot <- NAN_plot + geom_segment(aes(x =-50, y =-0.5, xend = -50, yend =30),size=1,linetype=1 )
NAN_plot <- NAN_plot + geom_point(data = name, aes(x = diff, y = fdr),fill="#f46d43",size=2,shape=21,alpha=0.7,stroke=0.5)#,fill="#009BFF") 

#NAN_plot <- NAN_plot + geom_text_repel(data = final  , aes(x = diff, y = fdr,label=cancer),
#                                       segment.color = "black",size =2,color ='black',min.segment.length = 0.4 , box.padding = 0.2)

NAN_plot <- NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,2,1,2),'lines'),
                             plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                             legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=12,face='plain'),legend.title=element_text(size=12,face='plain'),
                             axis.text.x=element_text(size=14,face='plain',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),
                             axis.title.x=element_text(size=16,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-100,100),breaks=seq(-80,80,40))+ scale_y_continuous(expand=c(0,0),limits=c(-0.05,5.2),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('[Diff RT 1p36 Loss-WT]')+ylab('FDR')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_TCGA_1p36_ppt.pdf", plot=AllTumorure_2,bg = 'white', width = 11, height = 8, units = 'cm', dpi = 600)


L1HS <- read.table("../Pancancer/L1HSpro_methylation_TCGA.txt",header = T)
myTest <- merge(L1HS,final)
cor.test(formula = ~ L1HSmethy + diff, data = myTest,method = "spearman")

new.plot<-ggplot() +theme_classic()
new.plot <- new.plot+  geom_smooth(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff)),method = "lm", se = F,formula= y ~ x,color="grey30",size=1)
new.plot <- new.plot+  geom_point(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff),fill=cancer),shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_text_repel(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff),label=cancer),size=1.8) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='plain'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=12,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
new.plot<- new.plot+ggtitle(NULL)+ylab(paste0("#LINE1 insertion diff\n(1p36 loss - no loss)"))+xlab(paste0("L1HS methylation"))
new.plot <- new.plot +scale_y_continuous(expand=c(0,0),limits = c(-5,100),breaks=seq(0,400,20)) +scale_x_continuous(expand=c(0,0),limits = c(0.65,0.95),breaks = seq(0.5,1,0.1)) +scale_fill_manual(name=NULL,values =c(Germinoma="#FFA097",IMT="#CBCE00",Mix="#00E6AA",Loss="#2166ac",YST="#FF99FF",WT="grey90"))
plot7 <- cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./L1HSmethy_L1_diff_pancancer.pdf"), plot=plot7,bg = 'white', width = 6.9, height = 7, units = 'cm', dpi = 600)
###############
cor.test(formula = ~ L1HSmethy + diff, data = myTest,method = "spearman")



table1 <- TCGA
table1$Tag <- table1$chr1p36
table(table1$chr1p36)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + stat_summary(fun = mean, geom = "point",   shape = 18, size = 3, color = "red")
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
#NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=TCGA$subtype),shape=21,width = 0.3,size=1,alpha=0.25,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=9,face='plain'),
                            legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.text=element_text(size=9,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + coord_cartesian(ylim = c(0, 20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = gg_color_hue(31))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/TCGA_RT_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)



##########
table1 <- TCGA[which(TCGA$subtype=="UCEC"),]
table1$Tag <- table1$chr1p36
table(table1$chr1p36)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + stat_summary(fun = mean, geom = "point",   shape = 18, size = 3, color = "red")
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
#NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=TCGA$subtype),shape=21,width = 0.3,size=1,alpha=0.25,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=9,face='plain'),
                            legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.text=element_text(size=9,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + coord_cartesian(ylim = c(0, 20)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = gg_color_hue(31))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/UCEC_RT_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)



##########
table1 <- TCGA[which(TCGA$subtype=="LUSC"),]
table1$Tag <- table1$chr1p36
table(table1$chr1p36)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + stat_summary(fun = mean, geom = "point",   shape = 18, size = 3, color = "red")
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
#NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=TCGA$subtype),shape=21,width = 0.3,size=1,alpha=0.25,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=9,face='plain'),
                            legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.text=element_text(size=9,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + coord_cartesian(ylim = c(0, 200)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = gg_color_hue(31))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/LUSC_RT_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)




##########
table1 <- TCGA[which(TCGA$subtype=="HNSC"),]
table1$Tag <- table1$chr1p36
table(table1$chr1p36)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + stat_summary(fun = mean, geom = "point",   shape = 18, size = 3, color = "red")
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
#NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=TCGA$subtype),shape=21,width = 0.3,size=1,alpha=0.25,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=9,face='plain'),
                            legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.text=element_text(size=9,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + coord_cartesian(ylim = c(0, 300)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = gg_color_hue(31))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/HNSC_RT_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)




##########
################
TCGA <- read.table("../../Pancancer_L1/Final_CNV_RT_with_RNA.txt",header = T)
TCGA$RT <- TCGA$RT_burden
check <- TCGA[which(TCGA$chr1p36=="Loss"),]
temp <- as.data.frame(table(check$subtype))
TCGA <- TCGA[which(TCGA$subtype %in% temp$Var1[which(temp$Freq>2)]),]
cancer <- unique(TCGA$subtype)
p <- rep(NA,length(cancer))
rho <- rep(NA,length(cancer))
diff <- rep(NA,length(cancer))
for(i in c(1:length(cancer))){
  xx <- TCGA[which(TCGA$subtype==cancer[i]),]
  #tmp <- cor.test(formula = ~ active_loci_TPM + RT_burden, data = xx,method = "spearman")
  loss <- xx$RT[which(xx$chr1p36=="Loss")]
  WT <- xx$RT[which(xx$chr1p36=="WT")]
  p[i] <- wilcox.test(loss,WT)$p.value
  diff[i] <- mean(loss) - mean(WT)
}
final <- data.frame(cancer,p,diff)
final$fdr <- -log10(p.adjust(final$p,method = "fdr"))

name <- final[which(final$fdr> -log10(0.05)),]
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_point(data = final , aes(x = diff, y = fdr),fill="grey90",color="grey10",alpha=0.7,size=2, shape = 21,stroke=0.5) 
NAN_plot <- NAN_plot + geom_hline(aes(yintercept=-log10(0.05)),color='black',size=.5,linetype='dashed') 
NAN_plot <- NAN_plot + geom_vline(aes(xintercept=0),color='black',size=.5,linetype=2)
#NAN_plot <- NAN_plot + geom_segment(aes(x =-50, y =-0.5, xend = -50, yend =30),size=1,linetype=1 )
NAN_plot <- NAN_plot + geom_point(data = name, aes(x = diff, y = fdr),fill="#f46d43",size=2,shape=21,alpha=0.7,stroke=0.5)#,fill="#009BFF") 

NAN_plot <- NAN_plot + geom_text_repel(data = final  , aes(x = diff, y = fdr,label=cancer),
                                       segment.color = "black",size =2,color ='black',min.segment.length = 0.4 , box.padding = 0.2)

NAN_plot <- NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,2,1,2),'lines'),
                             plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                             legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=12,face='plain'),legend.title=element_text(size=12,face='plain'),
                             axis.text.x=element_text(size=14,face='plain',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),
                             axis.title.x=element_text(size=16,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-80,80),breaks=seq(-80,80,40))+ scale_y_continuous(expand=c(0,0),limits=c(-0.05,6.7),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('[Diff RT 1p36 Loss-WT]')+ylab('FDR')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_TCGA_L1exp_1p36_ppt.pdf", plot=AllTumorure_2,bg = 'white', width = 11, height = 8, units = 'cm', dpi = 600)



################
TCGA <- read.table("../../Pancancer_L1/Final_CNV_RT_with_RNA.txt",header = T)
table(TCGA$subtype)
cancer <- unique(TCGA$subtype)
p <- rep(NA,length(cancer))
rho <- rep(NA,length(cancer))
for(i in c(1:length(cancer))){
  xx <- TCGA[which(TCGA$subtype==cancer[i]),]
  tmp <- cor.test(formula = ~ active_loci_TPM + RT_burden, data = xx,method = "spearman")
  p[i] <- tmp$p.value
  rho[i] <- tmp$estimate
}
final <- data.frame(cancer,p,rho)



TCGA <- read.table("../../Pancancer_L1/Final_CNV_RT_with_RNA.txt",header = T)
TCGA <- TCGA[which(TCGA$chr1p36!="WT" & TCGA$TP53=="noalteration"),]
temp <- as.data.frame(table(TCGA$subtype))
TCGA <- TCGA[which(TCGA$subtype %in% temp$Var1[which(temp$Freq>1)]),]
cancer <- unique(TCGA$subtype)
p <- rep(NA,length(cancer))
rho <- rep(NA,length(cancer))
for(i in c(1:length(cancer))){
  xx <- TCGA[which(TCGA$subtype==cancer[i]),]
  tmp <- cor.test(formula = ~ active_loci_TPM + RT_burden, data = xx,method = "spearman")
  p[i] <- tmp$p.value
  rho[i] <- tmp$estimate
}
final <- data.frame(cancer,p,rho)



TCGA$RT <- TCGA$active_loci_TPM
TCGA <- TCGA[which(!is.na(TCGA$TP53)),]
TCGA$p53 <- "Mut"
TCGA$p53[which(TCGA$TP53=="noalteration")] <- "WT"
table(TCGA$p53)

compare_means(RT ~ p53,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
compare_means(RT_burden ~ p53,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )


table1 <- TCGA[which(TCGA$TP53=="noalteration"),]
table1$Tag <- table1$chr1p36
table(table1$chr1p36)
table(table1$p53)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + stat_summary(fun = mean, geom = "point",   shape = 18, size = 3, color = "red")
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
#NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=TCGA$subtype),shape=21,width = 0.3,size=1,alpha=0.25,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=9,face='plain'),
                            legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.text=element_text(size=9,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + coord_cartesian(ylim = c(0, 60)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = gg_color_hue(31))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/TCGA_TPM_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)


table1 <- TCGA[which(TCGA$subtype=="HNSC"),]
table1$Tag <- table1$chr1p36
table(table1$chr1p36)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + stat_summary(fun = mean, geom = "point",   shape = 18, size = 3, color = "red")
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
#NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=TCGA$subtype),shape=21,width = 0.3,size=1,alpha=0.25,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=9,face='plain'),
                            legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.text=element_text(size=9,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + coord_cartesian(ylim = c(0, 200)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = gg_color_hue(31))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/HNSC_TPM_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)




###########

table1 <- TCGA[which(TCGA$subtype=="LUSC"),]
table1$Tag <- table1$chr1p36
table(table1$chr1p36)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=RT)) + theme_classic() 
NAN_plot <- NAN_plot + stat_summary(fun = mean, geom = "point",   shape = 18, size = 3, color = "red")
NAN_plot <- NAN_plot + geom_boxplot(color="grey10",fill="transparent",width=0.5,size=0.5,outlier.shape = NA) 
#NAN_plot <- NAN_plot + geom_quasirandom(data=table1,aes(fill=TCGA$subtype),shape=21,width = 0.3,size=1,alpha=0.25,stroke=0.5, varwidth = T,dodge.width = 0)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.5,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=9,face='plain'),
                            legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.text=element_text(size=9,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + coord_cartesian(ylim = c(0, 200)) 
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = gg_color_hue(31))+scale_x_discrete()
NAN_plot<- NAN_plot +ylab("LINE1 insertion") +xlab(NULL)

the <- compare_means(RT ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="test/LUSC_TPM_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)



###########
TCGA<- TCGA[which(TCGA$TP53=="noalteration"),]
xx <- TCGA[which(TCGA$chr1p36=="Loss" & TCGA$subtype=="LUSC"),]

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = xx, aes(x = active_loci_TPM, y = RT_burden),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = xx,  aes(x = active_loci_TPM, y = RT_burden,fill = subtype),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1HS RT")+xlab("RNA expresssion of L1HS")
new.plot <- new.plot + scale_x_continuous(expand=c(0,1))+scale_y_continuous(expand=c(0,1))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=gg_color_hue(31))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="test/linear_TCGA_exp_vs_RT_LOSS.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ active_loci_TPM + RT_burden, data = xx,method = "spearman")


xx <- TCGA[which(TCGA$chr1p36=="WT" & TCGA$subtype=="LUSC"),]

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = xx, aes(x = active_loci_TPM, y = RT_burden),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = xx,  aes(x = active_loci_TPM, y = RT_burden,fill = subtype),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("L1HS RT")+xlab("RNA expresssion of L1HS")
new.plot <- new.plot + scale_x_continuous(expand=c(0,1))+scale_y_continuous(expand=c(0,1))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=gg_color_hue(31))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="test/linear_TCGA_exp_vs_RT_WT.pdf", plot=plot7,bg = 'white', width =6.9, height = 7, units = 'cm', dpi = 600)
cor.test(formula = ~ active_loci_TPM + RT_burden, data = xx,method = "spearman")



