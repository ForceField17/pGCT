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

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
sample <- read.table('results/Summary_all_SomaticAlter.WGS.txt',header = F)
driver <- read.table("results/driver.WGS.txt",header = F)

pdf(file = paste0("./All_AF_WGS.pdf"),width=7,height=3.5, onefile = T)

for(i in 1:nrow(sample)){
  ID <- sample$V1[i]
  MutCount <- sample$V2[i]
  gene.table <- read.table(paste0("../WGS/GCT/",ID,".vcf"))
  a<-density(gene.table$V7,bw=1)
  Ceil <- max(a$y)*2
  F1a.plot<-ggplot()+theme_classic()
  F1a.plot<-F1a.plot+ geom_histogram(data=gene.table,aes(x=V7,y=(..density..)),fill="grey",binwidth=2,color="grey20",size=0.2,alpha=0.65,position='stack')
  F1a.plot<-F1a.plot+ geom_density(data=gene.table ,aes(x=V7),color="black",size=0.5,alpha=1)
  DriverGene <- driver[which(driver$V1 == ID),]
  if(nrow(DriverGene)>0){
    F1a.plot<-F1a.plot+ geom_vline(data=DriverGene,aes(xintercept = V3,color=V4),linetype=2)
    F1a.plot<-F1a.plot+ geom_text(data=DriverGene,aes(x = V3,y = Ceil*0.9,label=V2),size=3,color="black")
  }
  F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                            plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color='black'),text=element_text(size=12,face='bold.italic'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=10,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,face="plain",color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
  F1a.plot<-F1a.plot+ggtitle(paste0(ID," ",sample$V3[i]," ",MutCount," somatic mutations"))+xlab("Allele frequency")+ylab('Density')+scale_fill_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[2],"#1aa3d9",gg_color_hue(5)[3],gg_color_hue(5)[5]))
  AF.mean=median(gene.table$V7)
  F1a.plot<-F1a.plot+geom_vline(xintercept= AF.mean,linetype=5, size=0.5,color="grey30") 
  F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits = c(0,Ceil))+scale_x_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0,100,20))#
  F1a.plot<-F1a.plot + scale_color_manual(NULL,values = c(GoF="#d01c8b",LoF="#4dac26"))
  grid.arrange(F1a.plot)
}

dev.off() 

#figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
#ggsave(file="./test_C002.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)
