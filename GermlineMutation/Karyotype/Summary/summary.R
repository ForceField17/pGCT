# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#raw data preprocessing
chr21 <- read.delim('results_chr21_depth_log2.txt',header=F)
chrX <- read.delim('results_chrX_depth_log2.txt',header=F)

sample <- read.delim("Sample_feature.txt",sep = "\t",header = T)
colnames(sample)=c("ID","Histology","Karyotype","Age","Gender","Race","Site")


colnames(chr21) <- c("ID","chr21")
colnames(chrX) <- c("ID","chrX")

new1<-merge(chr21,chrX,1,1)
new2<-merge(new1,sample,1,1)

abnormal <- new2[which(new2$Site != "No" ),]
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_point(data = new2, aes(x = chr21, y = chrX,fill=Gender,shape=Gender),size=2.4,alpha=0.85,stroke=0.4)  
                      # geom_point(data = abnormal, aes(x = chr21, y = chrX, fill=Site,shape=Gender),size=3.7,alpha=0.8) + 
                      # geom_hline(aes(yintercept=0.5),color='black',size=.5,linetype='dashed')  +  
                      # geom_hline(aes(yintercept=-0.5),color='black',size=.5,linetype='dashed') +
                      ## geom_vline(aes(xintercept=0.3),color='black',size=.5,linetype='dashed')

NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(0.4,2,0.4,2),'lines'),
                            plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(1,'cm'),legend.position='right',legend.text=element_text(size=12,face='plain'),legend.title=element_text(size=12,face='plain'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold',color='black'),
                            axis.title.x=element_text(size=16,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

NAN_plot <-NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.15,1.15),breaks=seq(-1,1,0.5))+ scale_x_continuous(expand=c(0,0),limits=c(-0.25,0.85),breaks=seq(0,1.2,0.3))
NAN_plot <-NAN_plot + ylab('chrX coverage log2 ratio')+xlab('chr21 coverage log2 ratio')+ggtitle(NULL) 
#NAN_plot <-NAN_plot + scale_fill_manual(name="Dom His",values=c("grey",gg_color_hue(5)[2],gg_color_hue(5)[4],gg_color_hue(5)[5]),guide = guide_legend(override.aes=list(size=4,shape=21)))
NAN_plot <-NAN_plot + scale_fill_manual(name="Site",values=c(Male='#b3cde3',Female='#fccde5'),guide = guide_legend(override.aes=list(size=4,shape=21)))
NAN_plot <-NAN_plot + scale_shape_manual(name="Gender",values=c(21,23),guide = guide_legend(override.aes=list(size=4),nrow=2))


NAN_plot

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")


ggsave(file="figure2C_new.pdf", plot=figure_2,bg = 'white', width = 13.5, height = 8, units = 'cm', dpi = 600)
