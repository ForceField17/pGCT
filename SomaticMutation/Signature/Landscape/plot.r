# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library("gridExtra")
library(ggplot2)
library("reshape2")
library(ggbeeswarm)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l =80, c = 90)[1:n]
}


data_summary <- function(x) {
  m <- mean(x)
  ymin <- max(m-sd(x),0)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


WGS.table <- read.table("./bootstrap/Final_contribution.txt",header = T,row.names = "CaseID")
Info <- read.table("../AFdistribution/results/Summary_all_SomaticAlter.WGS.txt",header = F)
WGS.table <- WGS.table[which(rownames(WGS.table) != "T135"),]
Info <- Info[which(Info$V1 != "T135"),]

WGS.new <- WGS.table
WGS.new$xx <- rownames(WGS.new)

#################################################plot
myTable <- melt(WGS.new,id = c("xx"))

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=myTable,aes(x=xx,y=value,fill=variable),alpha=1,width=0.75, color = "black",size=0.2,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.001),breaks = seq(0,1,0.2))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0,2,2,2),'lines'),axis.ticks.x = element_blank(),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",axis.line = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("relative contribution")+scale_x_discrete(position = "bottom")#+coord_flip()
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(SBS1="#BBF0D2",SBS5="#81D3A2",SBS17a="#92c5de",SBS17b="#2166ac",SBS18="#D57CBE",SBS31="#BCBD45",SBS39="#8D69B8",SBSyst="#FF7F0F"))
PartB <- revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="WGS_proportion.pdf", plot=figure_2,bg = 'white', width = 28, height = 10, units = 'cm', dpi = 600)



#################################################plot
Tumor <- Info[which(Info$V4=="IMT"),]
Table <- myTable[which(myTable$xx %in% Tumor$V1),]
Order <- Tumor[,c(1,2)]; colnames(Order) <- c("xx","MutNum")
Table <- merge(Table,Order,1,1)

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=Table,aes(x=reorder(xx,MutNum),y=value,fill=variable),alpha=1,width=1,color="grey40",size=0.2,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.001),breaks = seq(0,1,0.5))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0,0.2,2,2),'lines'),axis.ticks.x = element_blank(),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",axis.line = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("relative contribution")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(SBS1="#92c5de",SBS5="grey80",SBS17a="#BBF0D2",SBS17b="#81D3A2",SBS18="#D57CBE",SBS31="#BCBD45",SBS39="#8D69B8",SBSyst="#FF7F0F"))
PartA <- revp_plot; TscaleA <- nrow(Tumor)/82

###
Tumor <- Info[which(Info$V4=="Teratoma"),]
Table <- myTable[which(myTable$xx %in% Tumor$V1),]
Order <- Tumor[,c(1,2)]; colnames(Order) <- c("xx","MutNum")
Table <- merge(Table,Order,1,1)
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=Table,aes(x=reorder(xx,MutNum),y=value,fill=variable),alpha=1,width=1,color="grey40",size=0.2,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.001),breaks = seq(0,1,0.5))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0,0.2,2,0),'lines'),axis.ticks.x = element_blank(),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",axis.line = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("relative contribution")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(SBS1="#92c5de",SBS5="grey80",SBS17a="#BBF0D2",SBS17b="#81D3A2",SBS18="#D57CBE",SBS31="#BCBD45",SBS39="#8D69B8",SBSyst="#FF7F0F"))
PartB <- revp_plot; TscaleB <- nrow(Tumor)/82

###
Tumor <- Info[which(Info$V4=="Germinoma"),]
Table <- myTable[which(myTable$xx %in% Tumor$V1),]
Order <- Tumor[,c(1,2)]; colnames(Order) <- c("xx","MutNum")
Table <- merge(Table,Order,1,1)
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=Table,aes(x=reorder(xx,MutNum),y=value,fill=variable),alpha=1,width=1,color="grey40",size=0.2,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.001),breaks = seq(0,1,0.5))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0,0.2,2,0),'lines'),axis.ticks.x = element_blank(),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",axis.line = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("relative contribution")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(SBS1="#92c5de",SBS5="grey80",SBS17a="#BBF0D2",SBS17b="#81D3A2",SBS18="#D57CBE",SBS31="#BCBD45",SBS39="#8D69B8",SBSyst="#FF7F0F"))
PartC <- revp_plot; TscaleC <- nrow(Tumor)/82

###
Tumor <- Info[which(Info$V4=="Mix"),]
Table <- myTable[which(myTable$xx %in% Tumor$V1),]
Order <- Tumor[,c(1,2)]; colnames(Order) <- c("xx","MutNum")
Table <- merge(Table,Order,1,1)
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=Table,aes(x=reorder(xx,MutNum),y=value,fill=variable),alpha=1,width=1,color="grey40",size=0.2,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.001),breaks = seq(0,1,0.5))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0,0.2,2,0),'lines'),axis.ticks.x = element_blank(),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",axis.line = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("relative contribution")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(SBS1="#92c5de",SBS5="grey80",SBS17a="#BBF0D2",SBS17b="#81D3A2",SBS18="#D57CBE",SBS31="#BCBD45",SBS39="#8D69B8",SBSyst="#FF7F0F"))
PartD <- revp_plot; TscaleD <- nrow(Tumor)/82

###
Tumor <- Info[which(Info$V4=="YST"),]
Table <- myTable[which(myTable$xx %in% Tumor$V1),]
Order <- Tumor[,c(1,2)]; colnames(Order) <- c("xx","MutNum")
Table <- merge(Table,Order,1,1)
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=Table,aes(x=reorder(xx,MutNum),y=value,fill=variable),alpha=1,width=1,color="grey40",size=0.2,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.001),breaks = seq(0,1,0.5))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(0,1,2,0),'lines'),axis.ticks.x = element_blank(),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",axis.line = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("relative contribution")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(SBS1="#92c5de",SBS5="grey80",SBS17a="#BBF0D2",SBS17b="#81D3A2",SBS18="#D57CBE",SBS31="#BCBD45",SBS39="#8D69B8",SBSyst="#FF7F0F"))
PartE <- revp_plot; TscaleE <- nrow(Tumor)/82


#figure_1<-cbind(cbind(ggplotGrob(PartA),ggplotGrob(PartB),ggplotGrob(PartC),ggplotGrob(PartD) ,ggplotGrob(PartE)),size="last")
#panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
#
#figure_1$widths[5]  <- unit(TscaleA,'null')
#figure_1$widths[14] <- unit(TscaleB,'null')
#figure_1$widths[23] <- unit(TscaleC,'null')
#figure_1$widths[32] <- unit(TscaleD,'null')
#figure_1$widths[41] <- unit(TscaleE,'null')
##figure_1$heights[panels][1] <-   unit(1/2,'null')
##figure_1$heights[panels][2] <-   unit(1,'null')
#ggsave(file="final_WGS.pdf", plot=figure_1,bg = 'white', width = 32, height = 10, units = 'cm', dpi = 600)

#################################################plot
Tumor <- Info[which(Info$V4=="IMT"),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V2,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank)+1,max(Tumor$Rank)-1),c(median(Tumor$V2),median(Tumor$V2)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V2,fill=V4,color=V4),alpha=0.75,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(3,3500),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,0.6,0.2,2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks.x = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of SNVs")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
Part1 <- revp_plot

Tumor <- Info[which(Info$V4=="Teratoma"),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V2,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank)+1,max(Tumor$Rank)-1),c(median(Tumor$V2),median(Tumor$V2)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V2,fill=V4,color=V4),alpha=0.75,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(3,3500),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,0.6,0.2,0),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of SNVs")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
Part2 <- revp_plot

Tumor <- Info[which(Info$V4=="Germinoma"),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V2,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank)+1,max(Tumor$Rank)-1),c(median(Tumor$V2),median(Tumor$V2)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V2,fill=V4,color=V4),alpha=0.75,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(3,3500),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,0.6,0.2,0),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of SNVs")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
Part3 <- revp_plot

Tumor <- Info[which(Info$V4=="Mix"),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V2,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank)+1,max(Tumor$Rank)-1),c(median(Tumor$V2),median(Tumor$V2)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V2,fill=V4,color=V4),alpha=0.75,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(3,3500),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,0.6,0.2,0),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of SNVs")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
Part4 <- revp_plot

Tumor <- Info[which(Info$V4=="YST"),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V2,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank)+1,max(Tumor$Rank)-1),c(median(Tumor$V2),median(Tumor$V2)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V2,fill=V4,color=V4),alpha=0.75,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(3,3500),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,0.2,0),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),
                            axis.text.x=element_text(size=9,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of SNVs")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
Part5 <- revp_plot


figure_1<-rbind(cbind(ggplotGrob(Part1),ggplotGrob(Part2),ggplotGrob(Part3),ggplotGrob(Part4) ,ggplotGrob(Part5)),
                cbind(ggplotGrob(PartA),ggplotGrob(PartB),ggplotGrob(PartC),ggplotGrob(PartD) ,ggplotGrob(PartE)),size="last")
panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]

figure_1$widths[7]  <- unit(TscaleA,'null')
figure_1$widths[20] <- unit(TscaleB,'null')
figure_1$widths[33] <- unit(TscaleC,'null')
figure_1$widths[46] <- unit(TscaleD,'null')
figure_1$widths[59] <- unit(TscaleE,'null')
scaleA <- 1 ; scaleB <- 1
figure_1$heights[panels][1] <-   unit(scaleA,'null')
figure_1$heights[panels][2] <-   unit(scaleB,'null')
figure_1$heights[panels][3] <-   unit(scaleA,'null')
figure_1$heights[panels][4] <-   unit(scaleB,'null')
figure_1$heights[panels][5] <-   unit(scaleA,'null')
figure_1$heights[panels][6] <-   unit(scaleB,'null')
figure_1$heights[panels][7] <-   unit(scaleA,'null')
figure_1$heights[panels][8] <-   unit(scaleB,'null')
figure_1$heights[panels][9] <-   unit(scaleA,'null')
figure_1$heights[panels][10] <-  unit(scaleB,'null')
ggsave(file="final_WGS.pdf", plot=figure_1,bg = 'white', width = 28, height = 13, units = 'cm', dpi = 600)


#################################################################################################################
Ploidy <- read.table("../../Landscape/results/Final_ploidy.txt",header = T)
Ploidy <- Ploidy[,c(1,2,3)]; colnames(Ploidy) <- c("xx","purity","ploidy")

xx <- WGS.new[,c("SBS1","SBS5","xx")];  YY <- Info[,c(1,2,4,6)]; colnames(YY) <- c("xx","SNV","His","Age");
qq <- merge(xx,YY,by.y="xx");  qq$Num <-  qq$SNV;  qq <- qq#[which(qq$His %in% c("IMT","Teratoma") ),]
qq <- merge(qq,Ploidy,by.y="xx")
cor.test(formula = ~ Age + Num, data = qq,method = "spearman")
qq$clock <- qq$Num * (qq$SBS1+qq$SBS5)
cor.test(formula = ~ Age + clock, data = qq,method = "spearman")
qq$clock.adj <- log10(qq$clock / qq$ploidy +1)
cor.test(formula = ~ Age + clock.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + clock.adj, data = qq,method = "spearman")

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = Age, y = clock.adj),method = "lm", se = F,formula= y ~ x,color="grey40")
new.plot <- new.plot + geom_point(data = qq,  aes(x = Age, y = clock.adj,color = His),alpha=0.9,size=2)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = Age, y = clock.adj,label = xx),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=9,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("Number of SNVs per monoploid")+xlab("Onset age of patients ")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(-0.5,20),breaks = seq(0,20,5))+scale_y_continuous(expand=c(0,0),limits = c(0,3),breaks=seq(0,2000,1))
new.plot <- new.plot + scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_Monoploid_SBS1_All.pdf", plot=plot7,bg = 'white', width =10.5, height = 10, units = 'cm', dpi = 600)



xx <- WGS.new[,c("SBS1","SBS5","xx")];  YY <- Info[,c(1,2,4,6)]; colnames(YY) <- c("xx","SNV","His","Age");
qq <- merge(xx,YY,by.y="xx");  qq$Num <-  qq$SNV;  qq <- qq[which(qq$His %in% c("IMT","Teratoma") ),]
qq <- merge(qq,Ploidy,by.y="xx")
cor.test(formula = ~ Age + Num, data = qq,method = "spearman")
qq$clock <- qq$Num * (qq$SBS1+qq$SBS5)
cor.test(formula = ~ Age + clock, data = qq,method = "spearman")
qq$clock.adj <- qq$clock / qq$ploidy
cor.test(formula = ~ Age + clock.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + clock.adj, data = qq,method = "spearman")

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = Age, y = clock.adj),method = "lm", se = F,formula= y ~ x,color="grey40")
new.plot <- new.plot + geom_point(data = qq,  aes(x = Age, y = clock.adj,fill = His),shape=21,alpha=0.75,size=2.5)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = Age, y = clock.adj,label = xx),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("Number of SNVs per monoploid")+xlab("Onset age of patients ")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(-0.5,20),breaks = seq(0,20,5))+scale_y_continuous(expand=c(0,0),limits = c(-5,275),breaks=seq(0,2000,50))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_Monoploid_Age_Teratomas.pdf", plot=plot7,bg = 'white', width =10.5, height = 10, units = 'cm', dpi = 600)

#################

xx <- WGS.new[,c("SBS1","SBS5","xx")];  YY <- Info[,c(1,2,4,6)]; colnames(YY) <- c("xx","SNV","His","Age");
qq <- merge(xx,YY,by.y="xx"); qq$Num <-  qq$SNV ; qq <- qq[which(qq$His %in% c("Germinoma") ),]
qq <- merge(qq,Ploidy,by.y="xx")
cor.test(formula = ~ Age + Num, data = qq,method = "spearman")
qq$clock <- qq$Num * (qq$SBS1+qq$SBS5)
cor.test(formula = ~ Age + clock, data = qq,method = "spearman")
#qq$ploidy[which(qq$ploidy<2)] <- 2 
qq$num.adj <- qq$Num / qq$ploidy 
cor.test(formula = ~ Age + num.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + num.adj, data = qq,method = "spearman")

qq$clock.adj <- qq$clock / qq$ploidy 
cor.test(formula = ~ Age + clock.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + clock.adj, data = qq,method = "spearman")


new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = Age, y = clock.adj),method = "lm", se = F,formula= y ~ x,color="grey40")
new.plot <- new.plot + geom_point(data = qq,  aes(x = Age, y = clock.adj,fill = His),shape=21,alpha=0.75,size=2.5)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = Age, y = clock.adj,label = xx),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("Number of clock SNVs per monoploid")+xlab("Onset age of patients ")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(-0.5,20),breaks = seq(0,20,5))+scale_y_continuous(expand=c(0,0),limits = c(-5,275),breaks=seq(0,2000,50))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_Monoploid_Age_Germinoma.pdf", plot=plot7,bg = 'white', width =10.5, height = 10, units = 'cm', dpi = 600)


#################
xx <- WGS.new[,c("SBS1","SBS5","xx")];  YY <- Info[,c(1,2,4,6)]; colnames(YY) <- c("xx","SNV","His","Age");
qq <- merge(xx,YY,by.y="xx"); qq$Num <-  qq$SNV ; qq <- qq[which(qq$His %in% c("YST") ),]
qq <- merge(qq,Ploidy,by.y="xx")
cor.test(formula = ~ Age + Num, data = qq,method = "spearman")
qq$clock <- qq$Num * (qq$SBS5+qq$SBS1)
cor.test(formula = ~ Age + clock, data = qq,method = "spearman")
#qq$ploidy[which(qq$ploidy<2)] <- 2 
qq$clock.adj <- qq$clock / qq$ploidy # /qq$purity
cor.test(formula = ~ Age + clock.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + clock.adj, data = qq,method = "spearman")

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = Age, y = clock.adj),method = "lm", se = F,formula= y ~ x,color="grey40")
new.plot <- new.plot + geom_point(data = qq,  aes(x = Age, y = clock.adj,color = His),alpha=0.9,size=2)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = Age, y = clock.adj,label = xx),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=9,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("Number of clock SNVs per monoploid")+xlab("Onset age of patients ")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(-0.5,20),breaks = seq(0,20,5))+scale_y_continuous(expand=c(0,0),limits = c(-10,275),breaks=seq(0,2000,50))
new.plot <- new.plot + scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_Monoploid_Age_YST.pdf", plot=plot7,bg = 'white', width =10.5, height = 10, units = 'cm', dpi = 600)


#################
xx <- WGS.new[,c("SBS1","SBS5","xx")];  YY <- Info[,c(1,2,4,6)]; colnames(YY) <- c("xx","SNV","His","Age");
qq <- merge(xx,YY,by.y="xx"); qq$Num <-  qq$SNV ; qq <- qq[which(qq$His %in% c("Mix") ),]
qq <- merge(qq,Ploidy,by.y="xx")
cor.test(formula = ~ Age + Num, data = qq,method = "spearman")
qq$clock <- qq$Num * (qq$SBS5+qq$SBS1)
cor.test(formula = ~ Age + clock, data = qq,method = "spearman")
#qq$ploidy[which(qq$ploidy<2)] <- 2 
qq$clock.adj <- qq$clock / qq$ploidy 
cor.test(formula = ~ Age + clock.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + clock.adj, data = qq,method = "spearman")

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = Age, y = clock.adj),method = "lm", se = F,formula= y ~ x,color="grey40")
new.plot <- new.plot + geom_point(data = qq,  aes(x = Age, y = clock.adj,color = His),alpha=0.9,size=2)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = Age, y = clock.adj,label = xx),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=9,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("Number of clock SNVs per monoploid")+xlab("Onset age of patients ")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(3,20),breaks = seq(0,20,5))+scale_y_continuous(expand=c(0,0),limits = c(-10,485),breaks=seq(0,2000,100))
new.plot <- new.plot + scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_Monoploid_Age_Mix.pdf", plot=plot7,bg = 'white', width =10.5, height = 10, units = 'cm', dpi = 600)


#################
xx <- WGS.new[,c("SBS1","SBS5","xx")];  YY <- Info[,c(1,2,4,6)]; colnames(YY) <- c("xx","SNV","His","Age");
qq <- merge(xx,YY,by.y="xx"); qq$Num <-  qq$SNV #; qq <- qq[which(qq$His %in% c("Mix") ),]
qq <- merge(qq,Ploidy,by.y="xx")
cor.test(formula = ~ Age + Num, data = qq,method = "spearman")
qq$clock <- qq$Num * (qq$SBS5+qq$SBS1)
cor.test(formula = ~ Age + clock, data = qq,method = "spearman")
#qq$ploidy[which(qq$ploidy<2)] <- 2 
qq$clock.adj <- qq$clock / qq$ploidy 
cor.test(formula = ~ Age + clock.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + clock.adj, data = qq,method = "spearman")

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = Age, y = clock.adj,color=His),method = "lm", se = F,formula= y ~ x)
new.plot <- new.plot + geom_point(data = qq,  aes(x = Age, y = clock.adj,color = His),alpha=0.9,size=2)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = Age, y = clock.adj,label = xx),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=9,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("Number of clock SNVs per monoploid")+xlab("Onset age of patients ")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(-0.5,20),breaks = seq(0,20,5))+scale_y_continuous(expand=c(0,0),limits = c(-10,325),breaks=seq(0,2000,100))
new.plot <- new.plot + scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_Monoploid_Age_All.pdf", plot=plot7,bg = 'white', width =10.5, height = 10, units = 'cm', dpi = 600)




#################
xx <- WGS.new[,c("SBS1","SBS5","xx")];  YY <- Info[,c(1,2,4,6)]; colnames(YY) <- c("xx","SNV","His","Age");
qq <- merge(xx,YY,by.y="xx"); qq$Num <-  qq$SNV ; qq <- qq[which(qq$His %in% c("YST","Germinoma","Mix") & qq$xx!="SJ01" ),]
cor.test(formula = ~ Age + Num, data = qq,method = "spearman")

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = Age, y = Num),method = "lm", se = F,formula= y ~ x,color="grey40")
new.plot <- new.plot + geom_point(data = qq,  aes(x = Age, y = Num,color = His),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=9,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("Number of SNVs contribute from\nclock-like signatures (SBS1/SBS5)")+xlab("Onset age of patients ")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(-0.5,20),breaks = seq(0,20,5))+scale_y_continuous(expand=c(0,0),limits = c(-20,2200),breaks=seq(0,2000,500))
new.plot <- new.plot + scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_SNV_Age_YST.pdf", plot=plot7,bg = 'white', width =11.6, height = 10, units = 'cm', dpi = 600)

##################
xx <- WGS.new[,c("SBS18","SBSyst","xx")];  YY <- Info[,c(1,2,4,6)]; colnames(YY) <- c("xx","SNV","His","Age");
qq <- merge(xx,YY,by.y="xx"); qq$Num <-  qq$SNV 
cor.test(formula = ~ Age + Num, data = qq,method = "spearman")
qq$theGCT <- (qq$SBS18+qq$SBSyst) * qq$SNV

qq <- qq[which(qq$His=="YST"),]
cor.test(formula = ~ theGCT + SNV, data = qq,method = "spearman")
new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = SNV, y = theGCT),method = "lm", se = F,formula= y ~ x,color="grey40")
new.plot <- new.plot + geom_point(data = qq,  aes(x = SNV, y = theGCT,fill = His),shape=21,alpha=0.75,size=2.5)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("SBS-GCT/SBS18 mutations in YST")+xlab("Number of SNVs")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(-30,2800),breaks = seq(0,3000,500))+scale_y_continuous(expand=c(0,0),limits = c(-15,850),breaks=seq(0,2000,200))
new.plot <- new.plot + scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_SNV_SBSyst.pdf", plot=plot7,bg = 'white', width =10.5, height = 10, units = 'cm', dpi = 600)



# qq <- qq[which(qq$His %in% c("YST") ),]
qq$chr1p36 <- "aWT"
qq$chr1p36[which(qq$xx %in% c("C005","F017","F026","F036","F055","SJ01","SJ10","SJ25","C009","F031","F042","SJ32"))] <- "Loss"
qq$SBSxx <- qq$SBS18 +qq$SBSyst
wilcox.test(formula = SBSxx ~ chr1p36,data = qq,conf.level = 0.95 ,alternative = "two.sided" )

NAN_plot <- ggplot(data=qq,aes(c)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(data=qq,aes(x=chr1p36,y=SBSxx*100),width=0.3,size=0.5,alpha=0,outlier.shape = NA)+
                       geom_quasirandom(data=qq,aes(x=chr1p36,y=SBSxx*100,fill=His),shape=21,width = 0.3,size=2,alpha=0.9,stroke=0.3, varwidth = T,stat = "identity",dodge.width = 0)
NAN_plot <- NAN_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_text(size=9,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot+ggtitle(NULL)+ylab("Percentage of SBS18 and SBS-YST")+xlab(NULL)
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2,80),breaks = seq(0,200,20) )
NAN_plot <- NAN_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
NAN_plot <- NAN_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))

plotxxx<-cbind(ggplotGrob(NAN_plot ),size="first")
ggsave(file=paste0("./SBS18_and_SBSyst.pdf"), plot=plotxxx,bg = 'white', width = 6, height = 9, units = 'cm', dpi = 600)

the <- compare_means(SBSxx ~ chr1p36,data = qq,p.adjust.method ="BH" )
the 


qq <- qq[which(qq$His %in% c("YST") ),]
qq$chr1p36 <- "aWT"
qq$chr1p36[which(qq$xx %in% c("C005","F017","F026","F036","F055","SJ01","SJ10","SJ25","C009","F031","F042","SJ32"))] <- "Loss"
qq$SBSxx <- qq$SBS18 +qq$SBSyst
wilcox.test(formula = SBSxx ~ chr1p36,data = qq,conf.level = 0.95 ,alternative = "two.sided" )

NAN_plot <- ggplot(data=qq,aes(c)) + theme_classic() 
NAN_plot <- NAN_plot + geom_boxplot(data=qq,aes(x=chr1p36,y=SBSxx*100),width=0.4,size=0.5,alpha=0,outlier.shape = NA)+
  geom_quasirandom(data=qq,aes(x=chr1p36,y=SBSxx*100,fill=His),shape=21,width = 0.3,size=2,alpha=0.9,stroke=0.3, varwidth = T,stat = "identity",dodge.width = 0)
NAN_plot <- NAN_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                            axis.text.x=element_text(size=9,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot+ggtitle(NULL)+ylab("Percentage of SBS18 and SBS-YST")+xlab(NULL)
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2,80),breaks = seq(0,200,20) )
NAN_plot <- NAN_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))
NAN_plot <- NAN_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#BEBAD8",Teratoma="#1aa3d9",YST="#FF9FFF"))

plotxxx<-cbind(ggplotGrob(NAN_plot ),size="first")
ggsave(file=paste0("./YST_SBS18_and_SBSyst.pdf"), plot=plotxxx,bg = 'white', width = 6, height = 8, units = 'cm', dpi = 600)

the <- compare_means(SBSxx ~ chr1p36,data = qq,p.adjust.method ="BH" )
the 


