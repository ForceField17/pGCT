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
  hcl(h = hues, l = 65, c = 100)[1:n]
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



################
TCGA <- read.table("../../Pancancer_L1/Final_CNV_RT.txt",header = T)
TCGA$RT <- TCGA$RT_burden
Methy <- read.table("../TCGA_samples_with_methy_and_CNV.txt",header = T)
TCGA <- merge(TCGA,Methy)
check_N_loss <- TCGA[which(TCGA$chr1p36=="Loss"),]
check_N_no_loss <- TCGA[which(TCGA$chr1p36!="Loss"),]
temp_N_loss <- as.data.frame(table(check_N_loss$subtype));colnames(temp_N_loss) <- c("cancerType","N_loss")
temp_N_no_loss <- as.data.frame(table(check_N_no_loss$subtype));colnames(temp_N_no_loss) <- c("cancerType","N_no_loss")
temp_N <- as.data.frame(table(TCGA$subtype));colnames(temp_N) <- c("cancerType","N")
Ntable <- merge(temp_N,temp_N_loss)
Ntable <- merge(Ntable,temp_N_no_loss)

TCGA <- TCGA[which(TCGA$subtype %in% Ntable$cancerType[which(Ntable$N>=10 & Ntable$N_loss>=10 & Ntable$N_no_loss>=10)]),]
table(TCGA$subtype,TCGA$chr1p36)

cancer <- unique(TCGA$subtype)
p <- rep(NA,length(cancer))
rho <- rep(NA,length(cancer))
diff <- rep(NA,length(cancer))
FC <- rep(NA,length(cancer))
for(i in c(1:length(cancer))){
  xx <- TCGA[which(TCGA$subtype==cancer[i]),]
  #tmp <- cor.test(formula = ~ active_loci_TPM + RT_burden, data = xx,method = "spearman")
  loss <- xx$RT[which(xx$chr1p36=="Loss")]
  WT <- xx$RT[which(xx$chr1p36=="WT")]
  p[i] <- wilcox.test(loss,WT)$p.value
  diff[i] <- mean(loss) - mean(WT)
  FC[i] <- log2( (mean(loss)+1) / (mean(WT)+1) )
}
final <- data.frame(cancer,p,diff,FC)
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

NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-100,100),breaks=seq(-80,80,40))+ scale_y_continuous(expand=c(0,0),limits=c(-0.05,5.3),breaks=seq(-100,200,1))
NAN_plot <- NAN_plot + xlab('[Diff RT 1p36 Loss-WT]')+ylab('FDR')+ggtitle(NULL) 

AllTumorure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Volcano_TCGA_1p36_ppt.pdf", plot=AllTumorure_2,bg = 'white', width = 11, height = 8, units = 'cm', dpi = 600)



meth <-  data.frame(aggregate(TCGA$MeanL1HSmethy, list(TCGA$subtype), FUN=mean))
colnames(meth) <- c("cancer","L1HSmethy")
myTest <- merge(meth,final)
myTest_diff <- myTest[which(myTest$cancer %in% "TGCT"),]

new.plot<-ggplot() +theme_classic()
new.plot <- new.plot+  geom_smooth(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff)),method = "lm", se = F,formula= y ~ x,color="grey30",size=1)
new.plot <- new.plot+  geom_point(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff)),fill="grey60",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_point(data = myTest_diff, aes(x = L1HSmethy, y = as.numeric(diff)),fill="#f46d43",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_text_repel(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff),label=cancer),size=1.8) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='plain'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=12,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
new.plot<- new.plot+ggtitle(NULL)+ylab(paste0("#LINE1 insertion \n(1p36 loss - no loss)"))+xlab(paste0("L1HS methylation"))
new.plot <- new.plot +scale_y_continuous(expand=c(0,1)) +scale_x_continuous(expand=c(0,0),limits = c(0.72,0.92),breaks = seq(0.5,1,0.05)) +scale_fill_manual(name=NULL,values =c(Germinoma="#FFA097",IMT="#CBCE00",Mix="#00E6AA",Loss="#2166ac",YST="#FF99FF",WT="grey90"))
plot7 <- cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./L1HSmethy_L1_diff_pancancer.pdf"), plot=plot7,bg = 'white', width = 7.7, height = 7, units = 'cm', dpi = 600)
###############
cor.test(formula = ~ L1HSmethy + diff, data = myTest,method = "spearman")
cor.test(formula = ~ L1HSmethy + diff, data = myTest,method = "pearson")

new.plot<-ggplot() +theme_classic()
new.plot <- new.plot+  geom_smooth(data = myTest, aes(x = L1HSmethy, y = as.numeric(FC)),method = "lm", se = F,formula= y ~ x,color="grey30",size=1)
new.plot <- new.plot+  geom_point(data = myTest, aes(x = L1HSmethy, y = as.numeric(FC)),fill="grey60",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_point(data = myTest_diff, aes(x = L1HSmethy, y = as.numeric(FC)),fill="#f46d43",shape=21,alpha=0.75,size=2.2) 
#new.plot <- new.plot+  geom_text_repel(data = myTest, aes(x = L1HSmethy, y = as.numeric(FC),label=cancer),size=1.8) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='plain'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=12,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
new.plot<- new.plot+ggtitle(NULL)+ylab(paste0("#LINE1 insertion log2FC\n(1p36 loss - no loss)"))+xlab(paste0("L1HS methylation"))
new.plot <- new.plot +scale_y_continuous(expand=c(0,0),limits = c(-0.8,2.2),breaks=seq(-1,2,0.5)) +scale_x_continuous(expand=c(0,0),limits = c(0.72,0.92),breaks = seq(0.5,1,0.05)) +scale_fill_manual(name=NULL,values =c(Germinoma="#FFA097",IMT="#CBCE00",Mix="#00E6AA",Loss="#2166ac",YST="#FF99FF",WT="grey90"))
plot7 <- cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./L1HSmethy_L1_FC_pancancer.pdf"), plot=plot7,bg = 'white', width = 7.7, height = 7, units = 'cm', dpi = 600)
###############
cor.test(formula = ~ L1HSmethy + FC, data = myTest,method = "spearman")
#cor.test(formula = ~ L1HSmethy + FC, data = myTest,method = "pearson")




meth <-  data.frame(aggregate(TCGA$MeanL1methy, list(TCGA$subtype), FUN=mean))
colnames(meth) <- c("cancer","L1methy")
myTest <- merge(meth,final)
myTest_diff <- myTest[which(myTest$cancer %in% name$cancer),]

new.plot<-ggplot() +theme_classic()
new.plot <- new.plot+  geom_smooth(data = myTest, aes(x = L1methy, y = as.numeric(FC)),method = "lm", se = F,formula= y ~ x,color="grey30",size=1)
new.plot <- new.plot+  geom_point(data = myTest, aes(x = L1methy, y = as.numeric(FC)),fill="grey60",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_point(data = myTest_diff, aes(x = L1methy, y = as.numeric(FC)),fill="#f46d43",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_text_repel(data = myTest, aes(x = L1methy, y = as.numeric(FC),label=cancer),size=1.8) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='plain'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=12,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
new.plot<- new.plot+ggtitle(NULL)+ylab(paste0("#LINE1 insertion log2FC\n(1p36 loss - no loss)"))+xlab(paste0("L1 methylation"))
new.plot <- new.plot +scale_y_continuous(expand=c(0,0),limits = c(-0.8,2.2),breaks=seq(-1,2,0.5)) +scale_x_continuous(expand=c(0,0),limits = c(0.62,0.83),breaks = seq(0.5,1,0.05)) +scale_fill_manual(name=NULL,values =c(Germinoma="#FFA097",IMT="#CBCE00",Mix="#00E6AA",Loss="#2166ac",YST="#FF99FF",WT="grey90"))
plot7 <- cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./L1methy_L1_FC_pancancer.pdf"), plot=plot7,bg = 'white', width = 7.7, height = 7, units = 'cm', dpi = 600)
###############
cor.test(formula = ~ L1methy + FC, data = myTest,method = "spearman")

new.plot<-ggplot() +theme_classic()
new.plot <- new.plot+  geom_smooth(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff)),method = "lm", se = F,formula= y ~ x,color="grey30",size=1)
new.plot <- new.plot+  geom_point(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff)),fill="grey60",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_point(data = myTest_diff, aes(x = L1HSmethy, y = as.numeric(diff)),fill="#f46d43",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_text_repel(data = myTest, aes(x = L1HSmethy, y = as.numeric(diff),label=cancer),size=1.8) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='plain'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=12,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
new.plot<- new.plot+ggtitle(NULL)+ylab(paste0("#LINE1 insertion diff\n(1p36 loss - no loss)"))+xlab(paste0("L1HS methylation"))
new.plot <- new.plot +scale_y_continuous(expand=c(0,0),limits = c(-5,100),breaks=seq(0,400,20)) +scale_x_continuous(expand=c(0,0),limits = c(0.65,0.95),breaks = seq(0.5,1,0.1)) +scale_fill_manual(name=NULL,values =c(Germinoma="#FFA097",IMT="#CBCE00",Mix="#00E6AA",Loss="#2166ac",YST="#FF99FF",WT="grey90"))
plot7 <- cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./L1HSmethy_L1_diff_pancancer.pdf"), plot=plot7,bg = 'white', width = 7.7, height = 7, units = 'cm', dpi = 600)
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

ggsave(file="./LUSC_RT_chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)



##########
table1 <- TCGA[which(TCGA$subtype=="TGCT"),]
table1 <- table1[which(table1$MeanL1HSmethy<=summary(table1$MeanL1HSmethy[which(table1$chr1p36=="Loss")])[5]),]
the <- compare_means(MeanL1HSmethy  ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the
the <- compare_means(RT  ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

table1$Tag <- table1$chr1p36
table(table1$chr1p36)

NAN_plot <- ggplot(data=table1,aes(x=reorder(Tag,-as.integer(as.factor(Tag))),y=MeanL1HSmethy )) + theme_classic() 
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

the <- compare_means(MeanL1HSmethy  ~ chr1p36,  data = table1,paired = F ,method = "wilcox.test",alternative = "two.sided" )
the

figure_2<-cbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="./LUSC_L1HSmethy _chr1p36.pdf", plot=figure_2,bg = 'white',  width =6, height = 7, units = 'cm', dpi = 600)


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



