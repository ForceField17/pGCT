 # WatchDHL
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


Clinical <- read.table("../../../lib/Clinical.txt",sep = "\t",header = T)
Focal <- read.table("../../Landscape/results/Final_Focal.txt",header=T)

Interaction <- merge(Clinical,Focal)
Interaction <- Interaction[,c("AID","Histology","Domin.His","chr1p36")]
Interaction$chr1p36[which(Interaction$chr1p36=="Gain" | Interaction$chr1p36=="WT")] <- "aWT"

############################## pie chart
theTable <- Interaction[which(Interaction$Histology=="YST"),]
count.data <- as.data.frame(table(theTable$chr1p36))
count.data <- count.data[rev(order(count.data$Var1)),]
library("dplyr") 
library("plyr")
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "grey20",alpha=0.9) +
  coord_polar("y", start = 0)+
 # geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey95","#2166ac")) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="pie_YST.pdf", plot=figure_1,bg = 'white', width = 10, height = 6, units = 'cm', dpi = 600)

############################## pie chart
theTable <- Interaction[which(Interaction$Histology!="YST"),]
count.data <- as.data.frame(table(theTable$chr1p36))
count.data <- count.data[rev(order(count.data$Var1)),]
library("dplyr") 
library("plyr")
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "grey20",alpha=0.9) +
  coord_polar("y", start = 0)+
  #geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey95","#2166ac")) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="pie_nonYST.pdf", plot=figure_1,bg = 'white', width = 10, height =6, units = 'cm', dpi = 600)


##########################################
japan <- read.table("../../../Transcriptome_subtyping/RNA/Genetic_profiles_All.txt",header = T)
japan <- japan[which(japan$Cohort=="Japan"),]
Interaction <- japan[,c("case","Histology","Gender","Dominant.His","PRDM2_loss")]
Interaction$chr1p36 <- Interaction$PRDM2_loss

Interaction$chr1p36[which(Interaction$PRDM2_loss=="Gain" | Interaction$PRDM2_loss=="WT")] <- "aWT"

############################## pie chart
theTable <- Interaction[which(Interaction$Histology=="YST"),]
count.data <- as.data.frame(table(theTable$chr1p36))
count.data <- count.data[rev(order(count.data$Var1)),]
library("dplyr") 
library("plyr")
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "grey20",alpha=0.9) +
  coord_polar("y", start = 0)+
#  geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey95","#2166ac")) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="pie_YST_japan.pdf", plot=figure_1,bg = 'white', width = 10, height = 6, units = 'cm', dpi = 600)

############################## pie chart
theTable <- Interaction[which(Interaction$Histology!="YST"),]

count.data <- as.data.frame(table(theTable$chr1p36))
count.data$Freq[which(count.data$Var1=="aWT")] <- 20 #SJ24,SJ12 no RNAseq
count.data <- count.data[rev(order(count.data$Var1)),]
library("dplyr") 
library("plyr")
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "grey20",alpha=0.9) +
  coord_polar("y", start = 0)+
 # geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey95","#2166ac")) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="pie_nonYST_japan.pdf", plot=figure_1,bg = 'white', width = 10, height = 6, units = 'cm', dpi = 600)



############################## pie chart
Xiehe <- read.table("../../../CNV/Independent_cohort/Xiehe_results.txt",header = T)
Interaction <- Xiehe[which(!is.na(Xiehe$Log2R_PRDM2)),]
Interaction$chr1p36 <- "aWT"
Interaction$chr1p36[which(Interaction$Log2R_PRDM2 <= -0.3)] <- "Loss"

theTable <- Interaction
count.data <- as.data.frame(table(theTable$chr1p36))
count.data <- count.data[rev(order(count.data$Var1)),]

count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "grey20",alpha=0.9) +
  coord_polar("y", start = 0)+
  # geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey95","#2166ac")) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="pie_YST_xiehe.pdf", plot=figure_1,bg = 'white', width = 10, height = 6, units = 'cm', dpi = 600)

#######################################################
theTable <- Interaction[which(Interaction$Age > 20),]
count.data <- as.data.frame(table(theTable$chr1p36))
count.data <- count.data[rev(order(count.data$Var1)),]

count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
rvp_plot <- ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "grey20",alpha=0.9) +
  coord_polar("y", start = 0)+
  # geom_text(aes(x=1.8, y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2)*100,"%)") ), color = "black")+
  scale_fill_manual(name=NULL,values=c("grey95","#2166ac")) +
  theme_classic() +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.title.x=element_blank(),
        axis.ticks = element_blank(),plot.title = element_blank(),axis.title.y=element_blank())

figure_1<-rbind(ggplotGrob(rvp_plot),size="first")
ggsave(file="pie_YST_xiehe_adult.pdf", plot=figure_1,bg = 'white', width = 10, height = 6, units = 'cm', dpi = 600)
#######################################################


