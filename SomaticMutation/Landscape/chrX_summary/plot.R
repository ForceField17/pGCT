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
sample <- read.table('../results/All_XY.txt',header = T)

ID<-as.character(sample[,1])
case <- unique(sample$AID)



############################################################3
gene.table2 <- sample
comp.table <- data.frame(gene.table2$chrX,gene.table2$Gender)
colnames(comp.table) <- c("Site","Gender")
table2 <- as.data.frame(table(comp.table))
table2 <- table2 %>%
  mutate(Pos = cumsum(Freq) - 0.5*Freq)
table2

table2$xxx[1]  <- table2$Freq[1] + table2$Freq[2] + table2$Freq[3]
table2$xxx[2]  <- table2$Freq[1] + table2$Freq[2] + table2$Freq[3]
table2$xxx[3]  <- table2$Freq[1] + table2$Freq[2] + table2$Freq[3]
table2$xxx[4]  <- table2$Freq[4] + table2$Freq[5] + table2$Freq[6]
table2$xxx[5]  <- table2$Freq[4] + table2$Freq[5] + table2$Freq[6]
table2$xxx[6]  <- table2$Freq[4] + table2$Freq[5] + table2$Freq[6]

table2$per <- table2$Freq/table2$xxx
table2$CNV <- c("a.Gain","c.Loss","b.WT","a.Gain","c.Loss","b.WT")

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=reorder(Gender,as.integer(as.character.factor(Gender))),y=per,fill=CNV),alpha=0.75,width=0.65, color = "black",stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.005),breaks = seq(0,1,0.2))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,2,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="right",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=25,vjust=1,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Propotion of cases")+scale_x_discrete(position = "bottom")#+coord_flip()

revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(c.Loss='#2166ac',b.WT="grey",a.Gain='#b2182b'))  #c(Male='#41b6c4',Female='#e78ac3'))
revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="chrX_CNV.pdf", plot=figure_2,bg = 'white', width = 8, height = 9, units = 'cm', dpi = 600)



############################################################3
gene.table2 <- sample
comp.table <- data.frame(paste0(gene.table2$Karyotype),gene.table2$Gender)
colnames(comp.table) <- c("Site","Gender")
table2 <- as.data.frame(table(comp.table))
table2 <- table2 %>%
  mutate(Pos = cumsum(Freq) - 0.5*Freq)
table2 <- table2[which(table2$Freq>0),]
table2

a  <- sum(table2$Freq[which(table2$Gender=="Female")])
b  <- sum(table2$Freq[which(table2$Gender=="Male")])

table2$xxx <- c(rep(a,length(which(table2$Gender=="Female"))),rep(b,length(which(table2$Gender=="Male"))))


table2$per <- table2$Freq/table2$xxx
table2$CNV <- paste0(c("j","i","h","g","f","c","b","a","e","d"),table2$Site) #c("a.Gain","c.Loss","b.WT","a.Gain","c.Loss","b.WT")

revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=reorder(Gender,as.integer(as.character.factor(Gender))),y=per,fill=CNV),alpha=0.7,width=0.65, color = "black",stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.005),breaks = seq(0,1,0.2))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,2,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="right",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=25,vjust=1,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Propotion of cases")+scale_x_discrete(position = "bottom")#+coord_flip()

revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c("#b2182b", "#b2182b","#b2182b", "#d6604d", "#f4a582", "#fddbc7", "grey99", "#4393c3","#c51b7d","#de77ae","grey99","#4393c3"))  #c(Male='#41b6c4',Female='#e78ac3'))
revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="chrX_genotype.pdf", plot=figure_2,bg = 'white', width = 8, height = 9, units = 'cm', dpi = 600)

