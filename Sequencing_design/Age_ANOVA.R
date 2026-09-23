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
library(ggpubr)
library(ggbeeswarm)

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
sample$His[which(sample$Histology=="NGGCT")] <- "Notspecified"
sample$His[which(sample$Histology=="GCT")] <- "Notspecified"
sample$Site    <- sample$PrimarySite
sample$Site[which(sample$PrimarySite=="No")] <- "aNo"
###############3
comp.table <- data.frame(sample$His ,sample$Site)
colnames(comp.table) <- c("Histology","Site")
table1 <- as.data.frame(table(comp.table))
table1$count <- c(c(table(comp.table$Histology)),c(table(comp.table$Histology)),c(table(comp.table$Histology)))
table1$count[which(table1$Histology=="Notspecified")] <- 0
therank <- rank(table1$count)

##################
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table1,aes(x=reorder(Histology,-count),y=Freq,fill=Site),width=0.7, alpha=0.75,color = "black",stat='identity',position=position_stack())+scale_y_continuous(expand=c(0,0),limits=c(0,74),breaks = seq(0,80,10))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,2,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",
                            legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='bold',color='black'),
                            axis.text.x=element_text(size=12,angle=25,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=16,vjust=-4,hjust=0.5,face='bold',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Number of cases")#+scale_x_discrete(position = "bottom")
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(Extracranial="#a6dba0",Intracranial="#9970ab",aNo='grey70'))
revp_plot<-revp_plot+scale_x_discrete(labels=c(Germinoma='Germinoma',YST='Yolk Sac Tumor',IMT='Immature Teratoma' ,Mix='Mixed Germ Cell Tumor', Teratoma="Mature Teratoma",EmbryonalC="Embryonal Carcinoma",ChorioC= "Choriocarcinoma",nsGCT="NOS Germ Cell Tumor"))
revp_plot
figure_3<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="./Fig_number/Histology_Site_number.pdf", plot=figure_3,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)

##########################################################################3
comp.table <- data.frame(sample$His,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
count.data <- as.data.frame(table(comp.table$Gender[which(comp.table$Histology=="Germinoma")]))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),legend.position="none",
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="./Fig_number/pie_gender_Germinoma.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)

##########################################################################3
comp.table <- data.frame(sample$His,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
count.data <- as.data.frame(table(comp.table$Gender[which(comp.table$Histology=="YST")]))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),legend.position="none",
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="./Fig_number/pie_gender_YST.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)


##########################################################################3
comp.table <- data.frame(sample$His,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
count.data <- as.data.frame(table(comp.table$Gender[which(comp.table$Histology=="Mix")]))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),legend.position="none",
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="./Fig_number/pie_gender_Mix.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)


##########################################################################3
comp.table <- data.frame(sample$His,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
count.data <- as.data.frame(table(comp.table$Gender[which(comp.table$Histology=="IMT")]))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),legend.position="none",
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="./Fig_number/pie_gender_IMT.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)


##########################################################################3
comp.table <- data.frame(sample$His,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
count.data <- as.data.frame(table(comp.table$Gender[which(comp.table$Histology=="Teratoma")]))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),legend.position="none",
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="./Fig_number/pie_gender_Teratoma.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)



##########################################################################3
comp.table <- data.frame(sample$His,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
count.data <- as.data.frame(table(comp.table$Gender[which(comp.table$Histology=="EmbryonalC")]))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),legend.position="none",
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="./Fig_number/pie_gender_EmbryonalC.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)



##########################################################################3
comp.table <- data.frame(sample$His,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
count.data <- as.data.frame(table(comp.table$Gender[which(comp.table$Histology=="ChorioC")]))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),legend.position="none",
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="./Fig_number/pie_gender_ChorioC.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)


##########################################################################3
comp.table <- data.frame(sample$His,sample$Gender)
colnames(comp.table) <- c("Histology","Gender")
count.data <- as.data.frame(table(comp.table$Gender[which(comp.table$Histology=="Notspecified")]))
count.data <- count.data[rev(order(count.data$Var1)),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=1) +
  coord_polar("y", start = 0)+
  scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female')) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),legend.position="none",
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="./Fig_number/pie_gender_Notspecified.pdf", plot=figure_1,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)




gene.table <- sample[which(sample$Age != "No" & sample$PrimarySite!="No" ),]
gene.table$age <-  as.numeric(gene.table$Age)
gene.table <- gene.table[which(gene.table$age<21),]

gene.table2 <- gene.table
res.aov3 <- aov(age ~ Gender * PrimarySite * His, data = gene.table2)
summary(res.aov3)


