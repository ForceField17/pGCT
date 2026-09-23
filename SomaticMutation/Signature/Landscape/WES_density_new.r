# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
library(ggplot2)
library("gridExtra")
library(MutationalPatterns)
library(BSgenome)
library(ggpubr)
head(available.genomes())
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)


### 1. SNP vcf data loading
vcf_files <- list.files("../WES/GCT",pattern = ".vcf", full.names = TRUE)

sample_names <- c("G1","G10","G11","G12","G2","G4","G5","G7","G8","G9","M1","M2","M3","M4","NG10",
                  "NG11","NG12","NG2","NG3","NG5","NG6","NG7","NG8","NG9","SJ05","SJ06","SJ07","SJ15","SJ16","SJ22",
                  "SJ23","SJ24","SJ26","SJ27","SJ29","SJ35","SJ36","SJ37","SJ38","SJ39","SJ40")

# summary(vcfs)
Info <- read.table("../AFdistribution/results/Summary_all_SomaticAlter.WES.txt",header = F)
Info <- Info[which(Info$V1 %in% sample_names),]
tissue <- Info$V4

##################
Tumor <- Info[which(Info$V4 %in% c("IMT","Teratoma")),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V3,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank),max(Tumor$Rank)),c(median(Tumor$V3),median(Tumor$V3)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V3,fill=V4,color=V4),alpha=0.9,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(0.1,100),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,0,1,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks.x = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of point mutations")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
Part1 <- revp_plot; Tscale1 <- nrow(Tumor)/41


Tumor <- Info[which(Info$V4=="Germinoma"),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V3,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank),max(Tumor$Rank)),c(median(Tumor$V3),median(Tumor$V3)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V3,fill=V4,color=V4),alpha=0.9,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(0.1,100),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,0,2,0),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of SNVs")+scale_x_discrete(position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
Part3 <- revp_plot; Tscale3 <- nrow(Tumor)/41

Tumor <- Info[which(Info$V4=="YST" ),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V3,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank),max(Tumor$Rank)),c(median(Tumor$V3),median(Tumor$V3)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V3,fill=V4,color=V4),alpha=0.9,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(0.1,100),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,1,0),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of SNVs")+scale_x_discrete(expand=c(0,1),position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF",EmbryonalC="#BDBAD7"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF",EmbryonalC="#BDBAD7"))
Part4 <- revp_plot; Tscale4 <- nrow(Tumor)/41


Tumor <- Info[which(Info$V4=="Mix"),]; rownames(Tumor) <- c(1:nrow(Tumor))
Tumor <- Tumor[order(Tumor$V3,decreasing = F),]; Tumor$Rank <- c(1:nrow(Tumor))
xxx <- data.frame(c(min(Tumor$Rank),max(Tumor$Rank)),c(median(Tumor$V3),median(Tumor$V3)))
colnames(xxx) <- c("X","Y")
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_point(data=Tumor,aes(x=Rank,y=V3,fill=V4,color=V4),alpha=0.9,size=2,shape=21)+
  geom_path(data=xxx,aes(x=X,y=Y),alpha=0.85,size=1,color="#4d4d4d")+
  scale_y_continuous(expand=c(0,0),limits=c(0.1,100),breaks = c(1,10,100,1000),trans = "log10")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,0,1,0),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",axis.line = element_blank(),axis.ticks = element_blank(),
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_blank(),
                            axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_blank())
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of SNVs")+scale_x_discrete(expand=c(0,1),position = "bottom")#+coord_flip()
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF",EmbryonalC="#BDBAD7"))
revp_plot<-revp_plot+scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF",EmbryonalC="#BDBAD7"))
Part2 <- revp_plot; Tscale2 <- nrow(Tumor)/41

figure_1<-rbind(cbind(ggplotGrob(Part1),ggplotGrob(Part2),ggplotGrob(Part3),ggplotGrob(Part4) ),size="last")
panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
figure_1$widths[7]  <- unit(1,'null') #unit(Tscale1,'null') #
figure_1$widths[20] <- unit(1,'null') #unit(Tscale2,'null') #
figure_1$widths[33] <- unit(1,'null') #unit(Tscale3,'null') #
figure_1$widths[46] <- unit(1,'null') #unit(Tscale4,'null') #

ggsave(file="MutDensity_WES_new.pdf", plot=figure_1,bg = 'white', width = 12, height = 11, units = 'cm', dpi = 600)

Info$XX <- Info$V4
Info$XX[which(Info$V4=="IMT" | Info$V4=="Teratoma")] <- "Cluster1"
Info$XX[which(Info$V4=="Mix")] <- "Cluster2"
Info$XX[which(Info$V4=="Germinoma")] <- "Cluster3"
Info$XX[which(Info$V4=="YST" )] <- "Cluster4"

the <- compare_means(V3 ~ XX,data = Info,p.adjust.method ="BH" )
the 



Ploidy <- read.table("../../Landscape/results/Final_ploidy.txt",header = T)
Ploidy <- Ploidy[,c(1,2,3)]; colnames(Ploidy) <- c("xx","purity","ploidy")

YY <- Info[,c(1,2,3,4,5,6)]; colnames(YY) <- c("xx","SNV","Mut","His","Dom.His","Age");
qq <- YY[which(YY$His %in% c("IMT","Teratoma") ),]
qq <- merge(qq,Ploidy,by.y="xx"); qq <- qq[which(qq$xx %in% sample_names),]
cor.test(formula = ~ Age + Mut, data = qq,method = "spearman")
cor.test(formula = ~ Age + SNV, data = qq,method = "spearman")
qq$SNV.adj <- qq$SNV / qq$ploidy 
cor.test(formula = ~ Age + SNV.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + SNV.adj, data = qq,method = "spearman")
qq$clock.adj <- qq$Mut / qq$ploidy 
cor.test(formula = ~ Age + clock.adj, data = qq,method = "pearson")
cor.test(formula = ~ Age + clock.adj, data = qq,method = "spearman")

new.plot<-ggplot() + theme_classic()
new.plot <- new.plot + geom_smooth(data = qq, aes(x = Age, y = SNV.adj),method = "lm", se = F,formula= y ~ x,color="grey40")
new.plot <- new.plot + geom_point(data = qq,  aes(x = Age, y = SNV.adj,color = His),alpha=0.9,size=2)  
#new.plot <- new.plot + geom_text(data = qq,  aes(x = Age, y = clock.adj,label = xx),alpha=0.9,size=2)  
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,1,1,1),'lines'),title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='none',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+ylab("Number of SNV per monoploid")+xlab("Onset age of patients ")
new.plot <- new.plot + scale_x_continuous(expand=c(0,0),limits = c(-0.5,20),breaks = seq(0,20,5))+scale_y_continuous(expand=c(0,0),limits = c(-0.2,6.7),breaks=seq(0,2000,2))
new.plot <- new.plot + scale_color_manual(name=NULL,values=c(Germinoma="#FFA49C",IMT="#CBCE00",Mix="#00E6AA",Teratoma="#1aa3d9",YST="#FF9FFF"))
plot7<-cbind(ggplotGrob(new.plot),size="first")
ggsave(file="linear_Monoploid_Teratoma_WES.pdf", plot=plot7,bg = 'white', width =10.5, height = 10, units = 'cm', dpi = 600)
cor.test(formula = ~ Age + SNV.adj, data = qq,method = "spearman")
