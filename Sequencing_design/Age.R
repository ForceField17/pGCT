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


##########################################################################################################################

sample$His    <- sample$Histology
sample$His[which(sample$Histology=="NGGCT")] <- "*Notspecified"
sample$His[which(sample$Histology=="GCT")] <- "*Notspecified"





#######
comp.table <- sample
count.data <- as.data.frame(table(comp.table$Cohort))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) %>%
  arrange(percent)
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=0.3,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  #geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c(Xinhua='#dfc27d',CBTTC='#b3de69',StJude='#8dd3c7',Nature2014='#80b1d3',BASIC3='#80b1d3')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="Fig_number/pie_cohort.pdf", plot=figure_1,bg = 'white', width = 7, height = 5, units = 'cm', dpi = 600)

#########################
count.data <- data.frame(c("Yes","No"),c(50,229-50))
colnames(count.data) <- c("Var1","Freq")
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) %>%
  arrange(percent)
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1,linewidth = 1.2, size=1,stat = "identity", color = "white",alpha=0.75) +
  coord_polar("y", start = 0,direction = -1)+
  #geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c(No="grey80",Yes='#461DAC',StJude='#8dd3c7',Nature2014='#80b1d3',BASIC3='#80b1d3')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="Fig_number/pie_germline_mut.pdf", plot=figure_1,bg = 'white', width = 7, height = 5, units = 'cm', dpi = 600)


#########################
count.data <- data.frame(c("Intronic","Intergenic"),c(17,48))
colnames(count.data) <- c("Var1","Freq")
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) %>%
  arrange(percent)
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1,linewidth = 1.2, size=1,stat = "identity", color = "white",alpha=0.75) +
  coord_polar("y", start = 0,direction = -1)+
  #geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c(No="grey80",Yes='#461DAC',Intronic='#8dd3c7',Intergenic='#80b1d3',BASIC3='#80b1d3')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="Fig_number/pie_RT_region.pdf", plot=figure_1,bg = 'white', width = 7, height = 5, units = 'cm', dpi = 600)


#######
comp.table <- sample
count.data <- as.data.frame(table(comp.table$His))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=0.5,stat = "identity", color = "white",alpha=0.75) +
  coord_polar("y", start = 0)+
 # geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey","#fdb462","#bebada",gg_color_hue(5)[1],gg_color_hue(5)[2],"#1aa3d9",gg_color_hue(5)[3],gg_color_hue(5)[5] )) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="Fig_number/pie_Histology.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)

#######
comp.table <- sample[which(sample$PrimarySite == "Intracranial"),]
count.data <- as.data.frame(table(comp.table$His))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=0.5,stat = "identity", color = "white",alpha=0.75) +
  coord_polar("y", start = 0)+
  geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey",gg_color_hue(5)[1],gg_color_hue(5)[2],gg_color_hue(5)[4],gg_color_hue(5)[3] ,gg_color_hue(5)[5]  )) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="Fig_number/pie_Histology_Intracranial.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)

#######
comp.table <- sample[which(sample$PrimarySite == "Extracranial"),]
count.data <- as.data.frame(table(comp.table$His))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=0.5,stat = "identity", color = "white",alpha=0.75) +
  coord_polar("y", start = 0)+
  geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c(gg_color_hue(5)[1] ,gg_color_hue(5)[3] ,gg_color_hue(5)[4],gg_color_hue(5)[5]  )) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="Fig_number/pie_Histology_Extracranial.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)





### pie plot for primary site (all pGCT)
comp.table <- data.frame(sample$Site,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
comp.table$His <- comp.table$Histology
comp.table$His[which(comp.table$Histology=="No")] <- "*Notspecified"
comp.table$His[which(comp.table$Histology=="NOS")] <- "*Notspecified"
comp.table$His[which(comp.table$Histology=="Dys")] <- "*Notspecified"
comp.table$His[which(comp.table$Histology=="ventricle")] <- "Intracranial"
comp.table$His[which(comp.table$Histology=="suprasellar")] <- "Intracranial"
comp.table$His[which(comp.table$Histology=="Ovary")] <- "Ovarian"
comp.table$His[which(comp.table$Histology=="Testis")] <- "Testicular"
comp.table$His[which(comp.table$Histology=="Pelvis")] <- "Pelvic"
comp.table$His[which(comp.table$Histology=="Parotid")] <- "Parotid"
comp.table$His[which(comp.table$Histology=="Retroperitoneum")] <- "Retroperitoneal"

#data1c$his[which(data1c$yyy=="Parotid")] <- "Extracranial"
#data1c$his[which(data1c$yyy=="Retroperitoneum")] <- "Extracranial"

count.data <- as.data.frame(table(comp.table$His))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=0.1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  geom_text(aes(x=2, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey","#fb9a99","#a6dba0","#9970ab","#8da0cb","#ffd92f","#e78ac3")) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="Fig_number/pie_PrimarySite.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)


############################################################3

#####intracranial
gene.table2 <- sample[which(sample$PrimarySite!="No"),]
comp.table <- data.frame(gene.table2$PrimarySite,gene.table2$Gender)
colnames(comp.table) <- c("Site","Gender")
table2 <- as.data.frame(table(comp.table))
table2 <- table2 %>%
  mutate(Pos = cumsum(Freq) - 0.5*Freq)
table2

table2$xxx[1]  <- table2$Freq[1]
table2$xxx[2] <- table2$Freq[4]
table2$xxx[3]  <- table2$Freq[3]
table2$xxx[4] <- table2$Freq[2]

table2$per[1]  <- table2$Freq[1]/(table2$Freq[1]+table2$Freq[3])
table2$per[2] <-  table2$Freq[2]/(table2$Freq[2]+table2$Freq[4])
table2$per[3]  <- table2$Freq[3]/(table2$Freq[1]+table2$Freq[3])
table2$per[4] <-  table2$Freq[4]/(table2$Freq[2]+table2$Freq[4])

table2$Pos[1] <-       table2$per[3] + 0.5 * table2$per[1] 
table2$Pos[2] <-       0.5 * table2$per[4]
table2$Pos[3] <-       0.5 * table2$per[3] 
table2$Pos[4] <-       table2$per[4] + 0.5 * table2$per[2]


revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=reorder(Site,Freq),y=per,fill=Gender),alpha=0.7,width=0.65, color = "black",stat='identity',position=position_stack())+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.005),breaks = seq(0,1,0.2))
#revp_plot<-revp_plot+geom_text(data = table2,aes(x=reorder(Site,-Freq),y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,2,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="top",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='plain',color='black'),
                            axis.text.x=element_text(size=12,angle=25,vjust=1,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Propotion of cases")+scale_x_discrete(position = "bottom")#+coord_flip()
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'))  #c(Male='#41b6c4',Female='#e78ac3'))
revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="Fig_number/Site_gender.pdf", plot=figure_2,bg = 'white', width = 4.5, height = 12, units = 'cm', dpi = 600)


#fisher exact test
Mu.FEtest <- cbind(c(0,0),c(0,0))
Mu.FEtest[1,1] <- table2$Freq[1]
Mu.FEtest[1,2] <- table2$Freq[2]
Mu.FEtest[2,1] <- table2$Freq[3]
Mu.FEtest[2,2] <- table2$Freq[4]
fisher.test(Mu.FEtest,alternative ="two.sided")
chisq.test(Mu.FEtest)


########################
#gene.table <- sample[which(sample$Age != "No"),]

gene.table <- sample[which(sample$Age != "No"),]
gene.table$age <-  as.numeric(gene.table$Age)
gene.table <- gene.table[which(gene.table$age<21),]
#My <- data.frame(gene.table$ID,gene.table$Gender,gene.table$Histology,gene.table$Position,gene.table$age,gene.table$Cohort)
#colnames(My) <- c("ID","Gender","His","Site","age","Cohort")
#Japan <- read.table("Japanese_51_pGCT_cohort.txt",header=T,sep="\t")
#Jap <- data.frame(Japan$sample_name,Japan$gender,Japan$Diagnosis,rep("Extracranial",nrow(Japan)),Japan$Age,rep("Japan",nrow(Japan)))
#colnames(Jap) <- c("ID","Gender","His","Site","age","Cohort")
#gene.table <-  rbind(My,Jap)

gene.table <- gene.table[which(gene.table$Histology!="EmbryonalC" & gene.table$Histology!="GCT"  & gene.table$Histology!="NGGCT" & gene.table$Histology!="ChorioC" & gene.table$Cohort!="Japan"),]

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_histogram(data=gene.table,aes(x=age,y=(..density..)/3,fill=His),binwidth=1,color="grey20",size=0.2,alpha=0.5,position='stack')
F1a.plot<-F1a.plot+ geom_density(data=gene.table ,aes(x=age,color=His),size=1,alpha=1)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="bottom",legend.text=element_text(size=10,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab("Age of disease diagnosis")+ylab('Density')
cdat <- ddply(gene.table, "His", summarise, AF.mean=median(age))
cdat
F1a.plot<-F1a.plot+geom_vline(data=cdat, aes(xintercept= AF.mean,colour=His),linetype=2, size=0.5) 
F1a.plot<-F1a.plot+scale_fill_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[2],gg_color_hue(5)[3],"#1aa3d9",gg_color_hue(5)[5]))
F1a.plot<-F1a.plot+scale_color_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[2],gg_color_hue(5)[3],"#1aa3d9",gg_color_hue(5)[5]))

F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(0,0.2))+scale_x_continuous(expand=c(0,0),limits=c(-0.5,20.5),breaks=seq(0,20,5))#

F1a.plot
figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="./Fig_Age/Age_vs_Histology_all.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)

the2 <- compare_means(age ~ Gender, data = gene.table)
the2
the2 <- compare_means(age ~ Site, data = gene.table)
the2
the2 <- compare_means(age ~ His, data = gene.table)
the2





##################
gene.table <- sample[which(sample$Age != "No"),]
gene.table$age <-  as.numeric(gene.table$Age)
gene.table <- gene.table[which(gene.table$age<21),]
gene.table <- gene.table[which(gene.table$Histology!="EmbryonalC" & gene.table$Histology!="GCT"  & gene.table$Histology!="NGGCT" & gene.table$Histology!="ChorioC" & gene.table$Cohort!="Japan"),]
gene.table2 <- gene.table[which(gene.table$PrimarySite=="Intracranial"),]

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_histogram(data=gene.table2,aes(x=age,y=(..density..)/3,fill=His),binwidth=1,color="grey20",size=0.3,alpha=0.5,position='stack')
F1a.plot<-F1a.plot+ geom_density(data=gene.table2 ,aes(x=age,color=His),size=0.75,alpha=0.9)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="right",legend.text=element_text(size=10,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab("Age of disease diagnosis")+ylab('Density')+scale_fill_manual(name=NULL,values=c(gg_color_hue(5) ))
cdat <- ddply(gene.table2, "His", summarise, AF.mean=mean(age))
cdat
F1a.plot<-F1a.plot+geom_vline(data=cdat, aes(xintercept= AF.mean,colour=His),linetype=2, size=0.5) 
F1a.plot<-F1a.plot+scale_fill_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[2],gg_color_hue(5)[3],"#1aa3d9",gg_color_hue(5)[5]))
F1a.plot<-F1a.plot+scale_color_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[2],gg_color_hue(5)[3],"#1aa3d9",gg_color_hue(5)[5]))

F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(0,0.25))+scale_x_continuous(expand=c(0,0),limits=c(-0.5,20.5),breaks=seq(0,20,5))#

F1a.plot
figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="./Fig_Age/Age_vs_Histology_Intracranial.pdf", plot=figure_1,bg = 'white', width = 22.5, height = 9, units = 'cm', dpi = 600)

the2 <- compare_means(age ~ Gender, data = gene.table2)
the2

the2 <- compare_means(age ~ His, data = gene.table2, method = "wilcox")
the2


#testin matrix




##################
gene.table <- sample[which(sample$Age != "No"),]
gene.table$age <-  as.numeric(gene.table$Age)
gene.table <- gene.table[which(gene.table$age<21),]
gene.table <- gene.table[which(gene.table$Histology!="EmbryonalC" & gene.table$Histology!="GCT"  & gene.table$Histology!="NGGCT" & gene.table$Histology!="ChorioC" & gene.table$Cohort!="Japan"),]
gene.table2 <- gene.table[which(gene.table$PrimarySite=="Extracranial"),]

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_histogram(data=gene.table2,aes(x=age,y=(..density..)/3,fill=His),binwidth=1,color="grey20",size=0.2,alpha=0.5,position='stack')
F1a.plot<-F1a.plot+ geom_density(data=gene.table2 ,aes(x=age,color=His),size=1,alpha=1)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="right",legend.text=element_text(size=10,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab("Age of disease diagnosis")+ylab('Density')+scale_fill_manual(name=NULL,values=c(gg_color_hue(5) ))
cdat <- ddply(gene.table2, "His", summarise, AF.mean=median(age))
cdat
F1a.plot<-F1a.plot+geom_vline(data=cdat, aes(xintercept= AF.mean,colour=His),linetype=2, size=0.5) 
F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(0,0.25))+scale_x_continuous(expand=c(0,0),limits=c(-0.5,20.5),breaks=seq(0,20,5))#
F1a.plot<-F1a.plot+scale_fill_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[3],gg_color_hue(5)[4],gg_color_hue(5)[5]))
F1a.plot<-F1a.plot+scale_color_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[3],gg_color_hue(5)[4],gg_color_hue(5)[5]))

F1a.plot
figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="./Fig_Age/Age_vs_Histology_Extracranial.pdf", plot=figure_1,bg = 'white', width = 14, height = 8, units = 'cm', dpi = 600)

the2 <- compare_means(age ~ Gender, data = gene.table2)
the2

the2 <- compare_means(age ~ His, data = gene.table2)
the2


library(ggpubr)


##################
gene.table <- sample[which(sample$Age != "No" & sample$PrimarySite!="No" & sample$His!="IMT" & sample$His!="Teratoma"),]
gene.table$age <-  as.numeric(gene.table$Age)
gene.table <- gene.table[which(gene.table$age<21),]
gene.table2 <- gene.table[which( gene.table$Histology!="GCT"  & gene.table$Histology!="NGGCT"  & gene.table$Cohort!="Japan"),]
gene.table2$class <- paste0(gene.table2$His,gene.table2$PrimarySite)

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_boxplot(data=gene.table2,aes(x=class,y=age),width=0.3,size=0.5,alpha=0,outlier.shape = NA)
F1a.plot<-F1a.plot+ geom_quasirandom(data=gene.table2,aes(x=class,y=age,color=His,fill=His),width = 0.25,size=1.7,alpha=0.5,stroke=0.8, varwidth = T)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="none",legend.text=element_text(size=14,hjust=0,face='italic'),
                          axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                          axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("Age"))+xlab(NULL)
F1a.plot<-F1a.plot + scale_y_continuous(expand=c(0,0),limits=c(0,21),breaks = seq(0,21,5)) 
the <- compare_means(age ~ class,  data = gene.table2)
the 

F1a.plot <- F1a.plot + scale_fill_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[3],gg_color_hue(5)[5]))
F1a.plot <- F1a.plot + scale_color_manual(name=NULL,values=c(gg_color_hue(5)[1],gg_color_hue(5)[3],gg_color_hue(5)[5]))
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./Fig_Age/All_His_Site.pdf"), plot=plotxxx,bg = 'white', width = 11, height = 10, units = 'cm', dpi = 600)

