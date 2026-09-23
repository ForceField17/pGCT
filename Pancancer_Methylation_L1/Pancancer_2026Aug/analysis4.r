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
library(dplyr)
library(ggpubr)
library(rstatix)
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


#CNV versus methylatio
#############################
RT <- read.table("../../Pancancer_L1/Final_CNV_RT.txt",header = T)
TCGA <- read.table("../../Pancancer_L1/Final_CNV_RT_with_RNA.txt",header = T)

myTable <- TCGA;myTable$cancerType <- myTable$subtype
# 过滤掉样本量过少的癌种（例如每组少于3个）
check_N_loss <- myTable[which(myTable$chr1p36=="Loss"),]
check_N_no_loss <- myTable[which(myTable$chr1p36!="Loss"),]
temp_N_loss <- as.data.frame(table(check_N_loss$cancerType));colnames(temp_N_loss) <- c("cancerType","N_loss")
temp_N_no_loss <- as.data.frame(table(check_N_no_loss$cancerType));colnames(temp_N_no_loss) <- c("cancerType","N_no_loss")
temp_N <- as.data.frame(table(myTable$cancerType));colnames(temp_N) <- c("cancerType","N")
Ntable <- merge(temp_N,temp_N_loss)
Ntable <- merge(Ntable,temp_N_no_loss)

myTable_filtered <- myTable[which(myTable$cancerType %in% Ntable$cancerType[which(Ntable$N_loss>=3 & Ntable$N_no_loss>=3)]),]
table(myTable_filtered$chr1p36,myTable_filtered$cancerType)

TCGA <- myTable_filtered
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

tmpTable <- data.frame(aggregate(myTable_filtered$active_loci_TPM, list(myTable_filtered$cancerType), FUN=mean))
colnames(tmpTable) <- c("cancer","L1exp")


myTest <- merge(tmpTable,final)
myTest_diff <- myTest[which(myTest$cancer %in% name$cancer),]

new.plot<-ggplot() +theme_classic()
new.plot <- new.plot+  geom_smooth(data = myTest, aes(x = L1exp, y = as.numeric(FC)),method = "lm", se = F,formula= y ~ x,color="grey30",size=1)
new.plot <- new.plot+  geom_point(data = myTest, aes(x = L1exp, y = as.numeric(FC)),fill="grey60",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_point(data = myTest_diff, aes(x = L1exp, y = as.numeric(FC)),fill="#f46d43",shape=21,alpha=0.75,size=2.2) 
new.plot <- new.plot+  geom_text_repel(data = myTest, aes(x = L1exp, y = as.numeric(FC),label=cancer),size=1.8) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='plain'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=12,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
new.plot<- new.plot+ggtitle(NULL)+ylab(paste0("#LINE1 insertion log2FC\n(1p36 loss - no loss)"))+xlab(paste0("L1 expression"))
new.plot <- new.plot +scale_y_continuous(expand=c(0,0),limits = c(-0.8,2.2),breaks=seq(-1,2,0.5)) +scale_x_continuous(expand=c(0,2)) +scale_fill_manual(name=NULL,values =c(Germinoma="#FFA097",IMT="#CBCE00",Mix="#00E6AA",Loss="#2166ac",YST="#FF99FF",WT="grey90"))
plot7 <- cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./L1HSexp_L1_FC_pancancer.pdf"), plot=plot7,bg = 'white', width = 7.7, height = 7, units = 'cm', dpi = 600)
###############
cor.test(formula = ~ L1exp  + FC, data = myTest,method = "spearman")









# 计算每个癌种的Wilcoxon检验P值
p_values <- myTable_filtered %>%
  group_by(cancerType) %>%
  wilcox_test(active_loci_TPM ~ chr1p36) %>%
  adjust_pvalue(method = "fdr") %>%  # 可选：FDR校正
  add_significance("p.adj") %>%
  mutate(
    p_label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "NS"
    ),
    # 用于在图中标注P值的位置（取该癌种Y轴最大值 + 5%余量）
    p_y_pos = max(myTable_filtered$active_loci_TPM[myTable_filtered$cancerType == cancerType]) * 1.08
  ) %>%
  select(cancerType, p_label, p_y_pos)

# 将P值信息合并回主数据
myTable_annotated <- myTable_filtered %>%
  left_join(p_values, by = "cancerType")

# 按癌种中位数排序（便于阅读）
cancer_order <- myTable_annotated %>%
  group_by(cancerType) %>%
  summarise(median_methy = median(active_loci_TPM, na.rm = TRUE)) %>%
  arrange(desc(median_methy)) %>%
  pull(cancerType)

myTable_annotated$cancerType <- factor(myTable_annotated$cancerType, levels = cancer_order)


# 定义颜色
loss_color <- "#2166ac"   # 红色 (loss)
  wt_color <- "grey90"     # 蓝色 (WT/no loss)
    
  # 绘制主图
  p <- ggplot(myTable_annotated, aes(x = cancerType, y = active_loci_TPM, color = chr1p36)) +
    # 抖动散点（半透明，避免重叠过多）
    geom_jitter(size = 0.4, alpha = 0.4, width = 0.2, aes(color = chr1p36)) +
    # 箱线图（只显示中位数和四分位，不显示离群点）
    geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.4, 
                 aes(fill = chr1p36), color = "black", linewidth = 0.3) +
    # 在癌种上方标注P值（用星号或 n.s.）
    geom_text(aes(y = p_y_pos, label = p_label), 
              size = 2, color = "black", vjust = 0.5, na.rm = TRUE) +
    # 颜色与填充
    scale_color_manual(values = c("Loss" = loss_color, "WT" = wt_color, "WT" = wt_color)) +
    scale_fill_manual(values = c("Loss" = loss_color, "WT" = wt_color, "WT" = wt_color)) +
    # 坐标轴标签
    labs(
      x = "Cancer Type",
      y = "L1HS exp",
      color = "1p36 Status",
      fill = "1p36 Status"
    ) +
    # 主题设置：X轴标签旋转90度，便于阅读
    theme_classic() +
    theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,2,1,2),'lines'),
          plot.title=element_text(size=12,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=9,face='plain'),legend.title=element_text(size=12,face='plain'),
          axis.text.x=element_text(size=9,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=9,vjust=0.5,hjust=1,face='plain',color='black'),
          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,hjust=0.5,vjust=2,face='plain',color='black')) +
    # 在癌种上方加一条细横线表示P值分组（可选）
    # 可以加入一条水平线标注整体中位数
    geom_hline(yintercept = median(myTable_annotated$active_loci_TPM, na.rm = TRUE), 
               linetype = "dashed", color = "gray30", linewidth = 0.3, alpha = 0.5)
  
  # 显示图形
  AllTumorure_2<-rbind(ggplotGrob(p),size="last")
  ggsave(file="Wilcon_TCGA_active_loci_TPM.pdf", plot=AllTumorure_2,bg = 'white', width = 24, height = 8, units = 'cm', dpi = 600)
  
  













TypeC <- "All TCGA 1p36 no loss"
test1 <- mergeTable[which(mergeTable$chr1p36!="Loss"),]
test1$RT <- log2(test1$RT_burden+1)
new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = test1, aes(x = MeanL1HSmethy, y = RT),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = test1,  aes(x = MeanL1HSmethy, y = RT,fill = chr1p36),alpha=0.75,size=1.5,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=11,vjust=0.5,hjust=0.5,face='plain',color='black'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='right',legend.text=element_text(size=11,face='plain'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(paste0(TypeC))+ylab("log2(#L1 ins. + 1)")+xlab("Methylation of L1HS promoter")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(0.45,1))+scale_y_continuous(expand=c(0,0.5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c("#2166ac","grey90"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./linear_methy_L1_",TypeC,".pdf"), plot=plot7,bg = 'white', width =9, height = 7.2, units = 'cm', dpi = 600)
cor.test(formula = ~ MeanL1HSmethy + RT, data = test1 ,method = "spearman")


TypeC <- "TCGA-LUSC"
test1 <- mergeTable[which(mergeTable$cancerType==TypeC),]
test1$RT <- log2(test1$RT_burden+1)
new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = test1, aes(x = MeanL1HSmethy, y = RT),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = test1,  aes(x = MeanL1HSmethy, y = RT,fill = chr1p36),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=11,vjust=0.5,hjust=0.5,face='plain',color='black'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='right',legend.text=element_text(size=11,face='plain'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(paste0(TypeC))+ylab("log2(#L1 ins. + 1)")+xlab("Methylation of L1HS promoter")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(0.45,1))+scale_y_continuous(expand=c(0,0.5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c("#2166ac","grey90"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./linear_methy_L1_",TypeC,".pdf"), plot=plot7,bg = 'white', width =9, height = 7.2, units = 'cm', dpi = 600)
cor.test(formula = ~ MeanL1HSmethy + RT, data = test1 ,method = "spearman")


TypeC <- "TCGA-TGCT"
test1 <- mergeTable[which(mergeTable$cancerType==TypeC & mergeTable$chr1p36=="Loss"),]
test1$RT <- log2(test1$RT_burden+1)
new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = test1, aes(x = MeanL1HSmethy, y = RT),color="grey30",method = "lm", se = F,formula= y ~ x,size=1)
new.plot <- new.plot + geom_point(data = test1,  aes(x = MeanL1HSmethy, y = RT,fill = chr1p36),alpha=0.75,size=2.2,shape=21)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = methy, y = L1HS,label = AID),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=11,vjust=0.5,hjust=0.5,face='plain',color='black'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='right',legend.text=element_text(size=11,face='plain'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(paste0(TypeC))+ylab("log2(#L1 ins. + 1)")+xlab("Methylation of L1HS promoter")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(0.45,1))+scale_y_continuous(expand=c(0,0.5))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c("#2166ac","grey90"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file=paste0("./linear_methy_L1_",TypeC,"_1p36loss.pdf"), plot=plot7,bg = 'white', width =9, height = 7.2, units = 'cm', dpi = 600)
cor.test(formula = ~ MeanL1HSmethy + RT, data = test1 ,method = "spearman")

