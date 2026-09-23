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
sample <- read.table('all_eroll_patients_2025.txt',header = T)

ID<-as.character(sample[,1])


case <- unique(sample$ID)


yyy <- rep("controlDNA",length(sample$ID))
data0<-data.frame(sample$ID,yyy,sample$Cohort,sample$Order)
colnames(data0)<-c("xxx","yyy","his","Order")
F0b.plot<-ggplot()#+theme_classic()
F0b.plot<-F0b.plot+geom_tile(data = data0,aes(x=reorder(xxx,Order),y="Cohort",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F0b.plot<-F0b.plot+scale_fill_manual(name=NULL,values=c(Xinhua='#dfc27d',CBTTC='#b3de69',StJude='#8dd3c7',Nature2014='#80b1d3',BASIC3='#80b1d3'),guide = guide_legend(override.aes=list(size=0.1),nrow=1),na.translate = F)
F0b.plot<-F0b.plot+theme(plot.margin=unit(c(0.15,4,0,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold'),
                         text=element_text(size=16,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='top',axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F0b.plot<-F0b.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
F0b.plot


#order <- c(1:nrow(sample))
yyy <- rep("Histology",length(sample$ID))
data1<-data.frame(sample$ID,yyy,sample$Gender,sample$Order)
colnames(data1)<-c("xxx","yyy","his","Order")
F1b.plot<-ggplot()#+theme_classic()
F1b.plot<-F1b.plot+geom_tile(data = data1,aes(x=reorder(xxx,Order),y="Gender",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F1b.plot<-F1b.plot+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5',Yes='black',No='grey90'),guide = NULL)
F1b.plot<-F1b.plot+theme(plot.margin=unit(c(0.15,4,0,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position=c(1.2,0.7),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.plot<-F1b.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
F1b.plot


data1c<-data.frame(sample$ID,sample$Site,sample$Order)
colnames(data1c)<-c("xxx","yyy","Order")
data1c$his <-    data1c$yyy
data1c$his[which(data1c$yyy=="No")] <- "*Not specified"
data1c$his[which(data1c$yyy=="NOS")] <- "*Not specified"
data1c$his[which(data1c$yyy=="Dys")] <- "*Not specified"
data1c$his[which(data1c$yyy=="ventricle")] <- "Intracranial"
data1c$his[which(data1c$yyy=="suprasellar")] <- "Intracranial"
#data1c$his[which(data1c$yyy=="Ovary")] <- "Ovarian/Pelvic"
#data1c$his[which(data1c$yyy=="Testis")] <- "Testicular"
#data1c$his[which(data1c$yyy=="Pelvis")] <- "Ovarian/Pelvic"
data1c$his[which(data1c$yyy=="Ovary")] <- "Extracranial"
data1c$his[which(data1c$yyy=="Testis")] <- "Extracranial"
data1c$his[which(data1c$yyy=="Pelvis")] <- "Extracranial"
data1c$his[which(data1c$yyy=="Parotid")] <- "Extracranial"
data1c$his[which(data1c$yyy=="Retroperitoneum")] <- "Extracranial"
F1c.plot<-ggplot()#+theme_classic()
F1c.plot<-F1c.plot+geom_tile(data = data1c,aes(x=reorder(xxx,Order),y="Primary site",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F1c.plot<-F1c.plot+scale_fill_manual(name=NULL,values=c('grey90',"#a6dba0","#9970ab"),guide = NULL)
F1c.plot<-F1c.plot+theme(plot.margin=unit(c(0.15,4,0,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position=c(1.2,0.7),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1c.plot<-F1c.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
F1c.plot

####
yyy <- rep("Histology",length(sample$ID))
data1x<-data.frame(sample$ID,yyy,sample$Histology,sample$Order)
colnames(data1x)<-c("xxx","yyy","his","Order")

F1x.plot<-ggplot()#+theme_classic()
F1x.plot<-F1x.plot+geom_tile(data = data1x,aes(x=reorder(xxx,Order),y="Histology",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F1x.plot<-F1x.plot+scale_fill_manual(name=NULL,values=c(Germinoma="#227C70",IMT="#227C70",Teratoma="#227C70",YST="#FF9FFF",Mix="#227C70",ChorioC="#227C70",EmbryonalC="#227C70",GCT="grey90",NGGCT="grey90"),guide = NULL)
#F1x.plot<-F1x.plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Teratoma="#1aa3d9",YST="#FF9FFF",Mix="#00E6AA",ChorioC="#227C70",EmbryonalC="#BDBAD7",GCT="grey90",NGGCT="grey90"),guide = NULL)
F1x.plot<-F1x.plot+theme(plot.margin=unit(c(0.15,4,0.15,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position=c(1.2,0.7),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1x.plot<-F1x.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
F1x.plot


yyy <- rep("Histology",length(sample$ID))
data2<-data.frame(sample$ID,yyy,sample$FatherBlood_WGS,sample$Order)
colnames(data2)<-c("xxx","yyy","his","Order")
F2b.plot<-ggplot()#+theme_classic()
F2b.plot<-F2b.plot+geom_tile(data = data2,aes(x=reorder(xxx,Order),y="Father normal",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F2b.plot<-F2b.plot+scale_fill_manual(name=NULL,values=c(Yes='#ff7f00',No='grey90'),guide = NULL)
F2b.plot<-F2b.plot+theme(plot.margin=unit(c(0.15,4,0,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position=c(1.2,0.7),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F2b.plot<-F2b.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)


yyy <- rep("Histology",length(sample$ID))
data3<-data.frame(sample$ID,yyy,sample$MotherBlood_WGS,sample$Order)
colnames(data3)<-c("xxx","yyy","his","Order")
F3b.plot<-ggplot()#+theme_classic()
F3b.plot<-F3b.plot+geom_tile(data = data3,aes(x=reorder(xxx,Order),y="Mother normal",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F3b.plot<-F3b.plot+scale_fill_manual(name=NULL,values=c(Yes='#ff7f00',No='grey90'),guide = NULL)
F3b.plot<-F3b.plot+theme(plot.margin=unit(c(0.15,4,0,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position=c(1.2,0.7),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F3b.plot<-F3b.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)



data4<-data.frame(sample$ID,sample$PatientBlood_WGS,sample$Order)
colnames(data4)<-c("xxx","WGS","Order")
data6<-data.frame(sample$ID,sample$PatientBlood_WES)
colnames(data6)<-c("xxx","WES")
data4 <- merge(data4,data6,1,1)

data4$his <- "No"
data4$his[which(data4$WES=="Yes")] <- "WES"
data4$his[which(data4$WGS=="Yes")] <- "WGS"

F4b.plot<-ggplot()#+theme_classic()
F4b.plot<-F4b.plot+geom_tile(data = data4,aes(x=reorder(xxx,Order),y="Patient normal",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F4b.plot<-F4b.plot+scale_fill_manual(name=NULL,values=c(WGS='#ff7f00',WES='#d53e4f',No='grey90'),guide = NULL)
F4b.plot<-F4b.plot+theme(plot.margin=unit(c(0.15,4,0,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position=c(1.2,0.7),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.plot<-F4b.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)




data5<-data.frame(sample$ID,sample$PatientTumor_WGS,sample$Order)
colnames(data5)<-c("xxx","WGS","Order")
data7<-data.frame(sample$ID,sample$PatientTumor_WES)
colnames(data7)<-c("xxx","WES")
data5 <- merge(data5,data7,1,1)

data5$his <- "No"
data5$his[which(data5$WES=="Yes")] <- "WES"
data5$his[which(data5$WGS=="Yes")] <- "WGS"

F5b.plot<-ggplot()#+theme_classic()
F5b.plot<-F5b.plot+geom_tile(data = data5,aes(x=reorder(xxx,Order),y="Patient tumor",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F5b.plot<-F5b.plot+scale_fill_manual(name=NULL,values=c(WGS='#ff7f00',WES='#d53e4f',No='grey90'),guide = NULL)
F5b.plot<-F5b.plot+theme(plot.margin=unit(c(0.15,4,0,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position=c(1.2,0.7),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.plot<-F5b.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)


yyy <- rep("Histology",length(sample$ID))
data8<-data.frame(sample$ID,yyy,sample$RNAseq,sample$Order)
colnames(data8)<-c("xxx","yyy","his","Order")
F8b.plot<-ggplot()#+theme_classic()
F8b.plot<-F8b.plot+geom_tile(data = data8,aes(x=reorder(xxx,Order),y="Patient tumor",fill=his),color=NA,width=1.1,height=1.2,stat='identity')
F8b.plot<-F8b.plot+scale_fill_manual(name=NULL,values=c(Yes='#2166ac',No='grey90'),guide = NULL)
F8b.plot<-F8b.plot+theme(plot.margin=unit(c(0.15,4,2,5),'lines'),panel.grid = element_blank(),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position=c(1.2,0.7),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F8b.plot<-F8b.plot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)




gene.scale1 <- 1/14
figure_1<-rbind(ggplotGrob(F0b.plot),ggplotGrob(F1b.plot),ggplotGrob(F1c.plot),ggplotGrob(F1x.plot),ggplotGrob(F2b.plot),ggplotGrob(F3b.plot),ggplotGrob(F4b.plot),ggplotGrob(F5b.plot),ggplotGrob(F8b.plot),size="first")
panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
figure_1$heights[panels][1] <- unit(gene.scale1,'null')
figure_1$heights[panels][2] <- unit(gene.scale1,'null')
figure_1$heights[panels][3] <- unit(gene.scale1,'null')
figure_1$heights[panels][4] <- unit(gene.scale1,'null')
figure_1$heights[panels][5] <- unit(gene.scale1,'null')
figure_1$heights[panels][6] <- unit(gene.scale1,'null')
figure_1$heights[panels][7] <- unit(gene.scale1,'null')
figure_1$heights[panels][8] <- unit(gene.scale1,'null')
figure_1$heights[panels][9] <- unit(gene.scale1,'null')


#grid.draw(figure_1)
ggsave(file="Fig1b_quartet.pdf", plot=figure_1,bg = 'white', width = 34, height = 8, units = 'cm', dpi = 600)





##########################################################################################################################
