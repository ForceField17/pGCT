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
library(DESeq2)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 80, c = 100)[1:n]
}

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}


Clinical <- read.table("./chr1p36_loss.txt",sep = "\t",header = T)
sortList <- Clinical[,c("ID","Order")]

###His
#####################################################################################################################################################################
table1 <- Clinical[which(Clinical$Cohort=="JP"),c("ID","Histology","Domin.His","Tumor_RNA")] 
rownames(table1) <- table1$ID;table1 <- table1[,-1]
Mut <- melt(as.matrix(table1))
colnames(Mut) <- c('ID','Feature','Value')
new_table1 <- merge(Mut,sortList,1,1);
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table1,aes(x=reorder(ID,Order),y=Feature,fill=Value),color='black',width=1,height=1,size=0.25,stat='identity')
F1A.pLot<-F1A.pLot+ scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF99FF",EmbryonalC="#BDBAD7",WES="#93C6E7",WGS="#54278f",TCS="#dfc27d",
                                                         Intracranial="#9970ab",Extracranial="#a6dba0",Male='#b3cde3',Female='#fccde5',Yes="#c51b7d",No="grey100",snRNA="white",Tiantna="#5aae61",RNAseq="#5aae61"))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_text(size=11,vjust=0.5,hjust=1,face='plain',color='black'),axis.line = element_blank(),
                         axis.text.x=element_text(size=10,angle = 90,vjust = 0.5,hjust = 1,face = "plain"),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
Part1 <- F1A.pLot
Part1

#####################################################################################################################################################################
table1 <- Clinical[which(Clinical$Cohort!="JP"),c("ID","Histology","Domin.His","Tumor_RNA")] 
rownames(table1) <- table1$ID;table1 <- table1[,-1]
Mut <- melt(as.matrix(table1))
colnames(Mut) <- c('ID','Feature','Value')
new_table1 <- merge(Mut,sortList,1,1);
F1A.pLot<-ggplot()+theme_classic()
F1A.pLot<-F1A.pLot+geom_tile(data = new_table1,aes(x=reorder(ID,Order),y=Feature,fill=Value),color='black',width=1,height=1,size=0.25,stat='identity')
F1A.pLot<-F1A.pLot+ scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF99FF",EmbryonalC="#BDBAD4",WES="#93C6E7",WGS="#54278f",TCS="#dfc27d",
                                                         Intracranial="#9970ab",Extracranial="#a6dba0",Male='#b3cde3',Female='#fccde5',Yes="#c51b7d",No="grey100",snRNA="grey90",Tiantna="#5aae61",RNAseq="#5aae61"))
F1A.pLot<-F1A.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0.1,0.1,0.1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='plain'),axis.text.y=element_text(size=11,vjust=0.5,hjust=1,face='plain',color='black'),axis.line = element_blank(),
                         axis.text.x=element_text(size=10,angle = 90,vjust = 0.5,hjust = 1,face = "plain"),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1A.pLot<-F1A.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
Part2 <- F1A.pLot
Part2

figure_5 <- rbind(ggplotGrob(Part1),ggplotGrob(Part2),size="last")

ggsave(file="heatmap_chr1p36.pdf", plot=figure_5,bg = 'white', width = 8, height = 5, units = 'cm', dpi = 600)
#####################################################################################################################################################################
