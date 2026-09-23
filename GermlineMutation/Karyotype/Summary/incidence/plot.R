# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
table2 <- read.table('trisomy21.txt',header = T)

table2$Fraction <- table2$Num/table2$Total

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=His,y=Fraction,fill=Mut,color=Mut),width=0.5,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.05),breaks = seq(0,0.06,0.01))+
  scale_x_discrete(labels = c(Chro1="population",Chro2="pGCT patients",Chro22="",Chro23="",Chro3="male population",Chro4="male pGCT patients",Chro45="",Chro44="",Chro5="female population",Chro6="female pGCT patients"))
  
#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=14,angle=0,vjust=0.5,hjust=0.5,face='bold',color='black'),
                            axis.text.x=element_text(size=14,angle=25,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=16,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Incidence")
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=gg_color_hue(4)[1],xFragile=gg_color_hue(4)[3],yFragile=gg_color_hue(4)[2],zFragile=gg_color_hue(4)[4]),labels=c(wFragile='Trisomy21 ',xFragile='XXY        ',yFragile='XO    ',zFragile='Swyer  '),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=gg_color_hue(4)[1],xFragile=gg_color_hue(4)[3],yFragile=gg_color_hue(4)[2],zFragile=gg_color_hue(4)[4]),labels=c(wFragile='Trisomy21 ',xFragile='XXY        ',yFragile='XO    ',zFragile='Swyer  '),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)

revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="All4.pdf", plot=figure_2,bg = 'white', width = 12, height = 14, units = 'cm', dpi = 600)



####################

1-pbinom(3, size=229, prob=(1/1000))
1-binom.test(3,229,1/1000,alternative = "less")$p.value

table <-  table2[c(1,2),]
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table,aes(x=reorder(His,-Fraction),y=Fraction),fill=c("grey40","#9B8EE6"),color=c("grey40","#9B8EE6"),width=0.6,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.0136),breaks = seq(0,0.06,0.004))+
  scale_x_discrete(labels = c(Chro1="population",Chro2="pGCT",Chro22="",Chro23="",Chro3="male population",Chro4="male pGCT",Chro45="",Chro44="",Chro5="female population",Chro6="female pGCT"))

revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=14,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=14,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=16,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Incidence")+coord_flip()

revp_plot1 <- revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="Trisomy21.pdf", plot=figure_2,bg = 'white',width = 10, height = 4.5, units = 'cm', dpi = 600)

binom.test(table[2,3],table[2,4],table[1,3]/table[1,4],alternative = "two.sided")

#########################333
table <-  table2[c(4,5),]
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table,aes(x=reorder(His,-Fraction),y=Fraction),fill=c("grey40","#9B8EE6"),color=c("grey40","#9B8EE6"),width=0.6,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.044),breaks = seq(0,0.06,0.01))+
  scale_x_discrete(labels = c(Chro1="population",Chro2="pGCT",Chro22="",Chro23="",Chro3="population",Chro4="pGCT",Chro45="",Chro44="",Chro5="population",Chro6="pGCT"))

revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=14,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=14,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=16,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Incidence")+coord_flip()

revp_plot2 <- revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="XXY.pdf", plot=figure_2,bg = 'white',width = 10, height = 4.5, units = 'cm', dpi = 600)

binom.test(table[2,3],table[2,4],table[1,3]/table[1,4],alternative = "two.sided")



###############
table <-  table2[c(7,8),]

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table,aes(x=reorder(His,-Fraction),y=Fraction),fill=c("grey40","#9B8EE6"),color=c("grey40","#9B8EE6"),width=0.6,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.00925),breaks = seq(0,0.00925,0.003))+
  scale_x_discrete(labels = c(Chro1="population",Chro2="pGCT",Chro22="",Chro23="",Chro3="population",Chro4="pGCT",Chro45="",Chro44="",Chro5="population",Chro6="pGCT"))

revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=14,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=14,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=16,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Incidence")+coord_flip()

revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=gg_color_hue(4)[1],xFragile=gg_color_hue(4)[3],yFragile=gg_color_hue(4)[2]),labels=c(wFragile='Trisomy21 ',xFragile='XXY        ',yFragile='XO    '),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=gg_color_hue(4)[1],xFragile=gg_color_hue(4)[3],yFragile=gg_color_hue(4)[2]),labels=c(wFragile='Trisomy21 ',xFragile='XXY        ',yFragile='XO    '),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)

revp_plot3 <- revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="XO.pdf", plot=figure_2,bg = 'white', width = 10, height = 4.5, units = 'cm', dpi = 600)

binom.test(table[2,3],table[2,4],p=table[1,3]/table[1,4],alternative = "two.sided")




###############
table <-  table2[c(10,11),]

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table,aes(x=reorder(His,-Fraction),y=Fraction),fill=c("grey40","#9B8EE6"),color=c("grey40","#9B8EE6"),width=0.6,alpha=0.7,stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.0185),breaks = seq(0,0.0182,0.006))+
  scale_x_discrete(labels = c(Chro1="population",Chro2="pGCT",Chro22="",Chro23="",Chro3="population",Chro4="pGCT",Chro45="",Chro44="",Chro7="population",Chro8="pGCT"))

revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=14,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=14,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=16,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Incidence")+coord_flip()

revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(wFragile=gg_color_hue(4)[1],xFragile=gg_color_hue(4)[3],yFragile=gg_color_hue(4)[2]),labels=c(wFragile='Trisomy21 ',xFragile='XXY        ',yFragile='XO    '),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(wFragile=gg_color_hue(4)[1],xFragile=gg_color_hue(4)[3],yFragile=gg_color_hue(4)[2]),labels=c(wFragile='Trisomy21 ',xFragile='XXY        ',yFragile='XO    '),guide=guide_legend(override.aes=list(size=0.5),nrow=1),na.translate = F)

revp_plot4 <- revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="XY.pdf", plot=figure_2,bg = 'white', width = 10, height = 4.5, units = 'cm', dpi = 600)

binom.test(table[2,3],table[2,4],table[1,3]/table[1,4],alternative = "two.sided")





  

figure_2<-rbind(ggplotGrob(revp_plot1),ggplotGrob(revp_plot2),ggplotGrob(revp_plot3),ggplotGrob(revp_plot4),size="first")
ggsave(file="All4_bar.pdf", plot=figure_2,bg = 'white', width = 16, height = 17, units = 'cm', dpi = 600)



