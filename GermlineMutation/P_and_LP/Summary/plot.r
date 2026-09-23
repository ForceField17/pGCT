# pGCT
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
library(oncoprint)
library(pheatmap)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

table1 <- read.table("./results/Final_potential_UKBB.txt",sep="\t",  head=T)

table2 <- read.table("./results/SummaryMut_UKBB.txt",sep="\t",  head=T)

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot +
  geom_bar(data=table2,aes(x=reorder(Gene,rank),y=case,fill=Effect),alpha=1,color="grey20",size=0.4,width=0.6,stat='identity',position=position_stack()) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='bottom',legend.text=element_text(size=10,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='italic',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='italic',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_discrete() + scale_y_continuous(expand=c(0,0),limits=c(0,7.1),breaks = seq(0,7,1) )
NAN_plot<- NAN_plot +ylab("Number of cases") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values=c(A.Missense='#66c2a5',B.Inframe='#a6d854',E.Frameshift='#fc8d62',C.Splice=gg_color_hue(6)[5],D.Stop='#df65b0'),
                                                    labels=c(A.Missense='Missense',B.Inframe='In Frame Deletion',E.Frameshift='Frame Shift',C.Splice="Splice Site",D.Stop='Nonsense'),guide = guide_legend(override.aes=list(size=0.4),nrow=3),na.translate = F)
#NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot3 <- NAN_plot

figure_2<-rbind(ggplotGrob(NAN_plot3),size="last")

ggsave(file="Summary_of_Putative_Pathogenic_variants.pdf", plot=figure_2,bg = 'white', width =11.5, height = 11, units = 'cm', dpi = 600)


##############################################################################33
##########################################
tableA <- table1[which(table1$Alt_Tum!="-"),]
nrow(tableA)

pValue_Chiq <- rep(1,nrow(tableA))
pValue_Bino <- rep(1,nrow(tableA))
pValue_BinoB <- rep(1,nrow(tableA))
for(i in 1:nrow(tableA)){
  table <- data.frame(c(0,0),c(0,0))
  table[1,1] <- as.numeric(tableA$Ref_Nor[i]) 
  table[1,2] <- as.numeric(tableA$Alt_Nor[i])
  table[2,1] <- as.numeric(tableA$Ref_Tum[i])
  table[2,2] <- as.numeric(tableA$Alt_Tum[i])
  
  LinearR <- fisher.test(table)
  pValue_Chiq[i] <- LinearR$p.value
  
  pValue_Bino[i] <- binom.test(as.numeric(tableA$Alt_Tum[i]), (as.numeric(tableA$Alt_Tum[i]) + as.numeric(tableA$Ref_Tum[i])), 0.5,alternative = "two.sided")$p.value 
  pValue_BinoB[i] <- binom.test(as.numeric(tableA$Alt_Nor[i]), (as.numeric(tableA$Alt_Nor[i]) + as.numeric(tableA$Ref_Nor[i])), 0.5,alternative = "two.sided")$p.value 
  
}

tableA$pValue_Chiq <- pValue_Chiq
tableA$pValue_Bino <- pValue_Bino
tableA$pValue_BinoB <- pValue_BinoB
tableA$ad.p_Chiq <- p.adjust(tableA$pValue_Chiq,method ="bonferroni")
tableA$ad.p_Bino <- p.adjust(tableA$pValue_Bino,method = "bonferroni")
tableA$MAF_B <- as.numeric(tableA$Alt_Nor)/(as.numeric(tableA$Alt_Nor) + as.numeric(tableA$Ref_Nor))
tableA$MAF_T <- as.numeric(tableA$Alt_Tum)/(as.numeric(tableA$Alt_Tum) + as.numeric(tableA$Ref_Tum))

tableB <- as.data.frame(tableA)#[which(tableA$ad.p_Bino <= 0.05),])
tableB$rank <- rank(tableB$MAF_T)

tableC <- data.frame(c(rep("Blood",nrow(tableB)),rep("Tumor",nrow(tableB))) ,c(tableB$Gene_Name, tableB$Gene_Name), c(tableB$CaseID,tableB$CaseID), 
                     c(tableB$Amino_Acid_Change, tableB$Amino_Acid_Change), c(tableB$MAF_B, tableB$MAF_T) , c(tableB$rank,tableB$rank),c(tableB$ad.p_Bino,tableB$ad.p_Bino) )

colnames(tableC) <- c("Tissue","Gene","Patient","Mutation","MAF","Rank","Bino_Adjust_P")

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot +
  geom_bar(data=tableC,aes(x=reorder(paste0(Gene," (",Patient,") ",Mutation),-Rank),y=MAF*100,fill=Tissue),alpha=0.9,color="grey20",size=0.4,width=0.6,stat='identity',position=position_dodge()) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',legend.text=element_text(size=10,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='italic',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='italic',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_discrete() + scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks = seq(0,100,25))
NAN_plot<- NAN_plot +ylab("MAF (%)") +xlab(NULL)+ geom_hline(yintercept = 65,color="black",linetype=2)
NAN_plot<- NAN_plot +ylab("MAF (%)") +xlab(NULL)+ geom_hline(yintercept = 50,color="black",linetype=2)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values=c(Blood='#7fbc41',Tumor='#de77ae'),
                                         guide = guide_legend(override.aes=list(size=0.4),nrow=1),na.translate = F)
NAN_plot3 <- NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot3),size="last")
ggsave(file="Germline_Blood_Tumor_LOH_and_ROH.pdf", plot=figure_2,bg = 'white', width =19, height = 15, units = 'cm', dpi = 600)





table2 <- read.table("./results/Cellular_process.txt",sep="\t",  head=T)

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot +
  geom_point(data=table2[1,],aes(x=case,y=variant,color=Term,fill=Term),alpha=0.95,shape=21,stroke=1,size=22) +
  geom_point(data=table2[2,],aes(x=case,y=variant,color=Term,fill=Term),alpha=0.95,shape=21,stroke=1,size=8) +
  geom_point(data=table2[4,],aes(x=case,y=variant,color=Term,fill=Term),alpha=0.95,shape=21,stroke=1,size=8) +
  geom_point(data=table2[3,],aes(x=case,y=variant,color=Term,fill=Term),alpha=0.95,shape=21,stroke=1,size=6) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='bottom',legend.text=element_text(size=10,hjust=0,face='plain'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='italic',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='italic',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(0,25),breaks = seq(0,25,5) ) + scale_y_continuous(expand=c(0,0),limits=c(0,25),breaks = seq(0,25,5) )
NAN_plot<- NAN_plot +ylab("Number of variants") +xlab("Number of cases")
NAN_plot <- NAN_plot + scale_size(range = c(5,25),breaks = c(2,3,3,11),guide =F,transform = "identity")
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values=c("#B02418","#68349A","#2F6EBA","#595959"),guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =c("#E1C5C3","#D2C8DC","#DEE4EB","#EFEFEF"))
NAN_plot3 <- NAN_plot

figure_2<-rbind(ggplotGrob(NAN_plot3),size="last")

ggsave(file="CellularProcess.pdf", plot=figure_2,bg = 'white', width =11.5, height = 11, units = 'cm', dpi = 600)

#######################
count.data <- table2#[rev(order(table2$variant)),]
count.data$xxx <- c("d","a","b","c")
count.data <- count.data %>%
  mutate(percent = variant/sum(count.data$variant )) %>%
  arrange(xxx)

count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = round(percent,3), fill = xxx)) +theme_classic() + 
  geom_bar( width=1,size=1.2,stat = "identity", color = "grey10",alpha=1) +
  coord_polar(theta ="y", start = 0)+
  #geom_text(aes(x=1.8, y = percent, label = paste0(variant," (",round(percent,digits = 3)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("#dfc27d","#A69CD0","grey90","#fcc5c0")) +
  theme_void()

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="pie_CellularProcess.pdf", plot=figure_1,bg = 'white', width = 12, height = 10, units = 'cm', dpi = 600)
