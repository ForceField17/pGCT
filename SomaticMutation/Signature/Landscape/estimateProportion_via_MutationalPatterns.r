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
head(available.genomes())
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)


### 1. SNP vcf data loading
vcf_files <- list.files("../WGS/GCT",pattern = ".vcf", full.names = TRUE)

sample_names <- c("C001","C002","C003","C004","C005","C006","C007","C008","C009","C010","C010r","C011","C013","C014",
                  "F001","F002","F002r","F004","F005","F006","F007","F008","F010","F011","F012","F013","F014","F015",
                  "F016","F017","F018","F019","F020","F021","F022","F023","F024","F026","F027","F029","F031","F034",
                  "F036","F037","F038","F039","F041","F042","F043","F044","F045","F046","F048","F050","F052","F053",
                  "F054","F055","F056","F057","F058","F059","F059r","SJ01","SJ02","SJ03","SJ04","SJ08","SJ09","SJ10",
                  "SJ12","SJ13","SJ14","SJ17","SJ18","SJ19","SJ20","SJ25","SJ28","SJ31","SJ32","SJ33","T135")

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,type="snv")

# summary(vcfs)
Info <- read.table("../AFdistribution/results/Summary_all_SomaticAlter.WGS.txt",header = F)
Info$V1 == sample_names
tissue <- Info$V3


### 2. Base substitution types ( 6 types of single base substitution and 96 substitution profiles)
muts = mutations_from_vcf(grl[[1]])
head(muts, 12)

 # convert 12 to the 6 types
types = mut_type(grl[[1]])
head(types, 12)

 #retrieve one base upstream and one base downstream
context = mut_context(grl[[1]], ref_genome)
head(context, 12)

 #retrieve all
type_context = type_context(grl[[1]], ref_genome)
lapply(type_context, head, 12)


 #annotate CpG sites
type_occurrences <- mut_type_occurrences(grl, ref_genome)
type_occurrences


 ### 96 trinucleodide mutation profile 4 * 6 * 4
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)


### 4. Find optimal contribution of known signatures
mutName <- read.table("../StandardName.txt",header = T)
the <- read.table("../Purify_denovoSig/mySpectrum_finalized.txt",header = T,row.names = "MutationType")
Order <- mutName[order(mutName$Type3),] # put in standard order
mySpectrum <- as.matrix(the[Order$Type1,])
#rownames(mySpectrum) <- Order$Type1
#signatures <- get_known_signatures()

fit_res <- fit_to_signatures(mut_mat, mySpectrum)
pdf(file.path( "./bootstrap/Stricter_contribution1.pdf"), width = 24, height = 7)
plot_contribution(fit_res$contribution,coord_flip = FALSE,mode = "absolute")
dev.off()
plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, y_intercept = 0.95)


### Stricter refitting
strict_refit <- fit_to_signatures_strict(mut_mat, mySpectrum, max_delta = 0.01, method = "best_subset")
fig_list <- strict_refit$sim_decay_fig
fig_list[[1]]

fit_res_strict <- strict_refit$fit_res
pdf(file.path( "./bootstrap/Stricter_contribution2.pdf"), width = 24, height = 7)
plot_contribution(fit_res_strict$contribution,coord_flip = FALSE, mode = "absolute")
dev.off()
plot_original_vs_reconstructed(mut_mat, fit_res_strict$reconstructed, y_intercept = 0.95)



###############################################################################################################3333
WGS <- read.table("../WGS/bootstrap_contribution.txt",header = T,row.names = "CaseID")
Info <- read.table("../AFdistribution/results/Summary_all_SomaticAlter.WGS.txt",header = F)

WGS.table <- WGS[!(rownames(WGS) %in% c("SJ01","T135")),]
SJ01 <- WGS[rownames(WGS) %in% c("SJ01","T135"),]
WGS.tableA <- WGS.table[which( rownames(WGS.table) %in% Info$V1[Info$V2>100]),]
WGS.tableB <- WGS.table[which( rownames(WGS.table) %in% Info$V1[Info$V2<=100]),]

################### fitting for SNV <= 100, maxDelta = 0.1
strict_refit <- fit_to_signatures_strict(mut_mat[,rownames(WGS.tableB)], mySpectrum[,c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS39","SBSyst")], max_delta = 0.1, method = "best_subset")
fig_list <- strict_refit$sim_decay_fig
fit_res_strict <- strict_refit$fit_res
pdf(file.path( "./bootstrap/Stricter_contribution_LOWdensity_round1.pdf"), width = 24, height = 7)
plot_contribution(fit_res_strict$contribution,coord_flip = FALSE, mode = "absolute")
dev.off()
plot_original_vs_reconstructed(mut_mat[,rownames(WGS.tableB)], fit_res_strict$reconstructed, y_intercept = 0.95)

#Considering only SBS1 & SBS5 show cosine similarity > 0.6 in reconstruction, so SBS1 & SBS5 were remained
fit_res <- fit_to_signatures(mut_mat[,rownames(WGS.tableB)], mySpectrum[,c("SBS1","SBS5")])
pdf(file.path( "./bootstrap/Stricter_contribution_LOWdensity_round2.pdf"), width = 24, height = 7)
plot_contribution(fit_res$contribution,coord_flip = FALSE, mode = "absolute")
dev.off()
plot_original_vs_reconstructed(mut_mat[,rownames(WGS.tableB)], fit_res$reconstructed, y_intercept = 0.95)

yy <- t(fit_res$contribution)
AA <- rep(0,2)
for(i in c(1:2)){
  AA[i] <- yy[1,i]/rowSums(yy)[1]
}
for(j in c(2:ncol(mut_mat[,rownames(WGS.tableB)]))){
  BB <- rep(0,2)
  for(i in c(1:2)){
    BB[i] <- yy[j,i]/rowSums(yy)[j]
  }
  AA <- rbind(AA,BB)
}
AA <- as.data.frame(AA)
AA$VV3 <- 0 ;AA$VV4 <- 0 ;AA$VV5 <- 0 ;AA$VV6 <- 0 ;AA$VV7 <- 0 ;AA$VV8 <- 0 
rownames(AA) <- rownames(yy)
colnames(AA) <- c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS31","SBS39","SBSyst")

PartB <- rbind(AA,SJ01)


################### fitting for SNV > 100, maxDelta = 0.01
strict_refit <- fit_to_signatures_strict(mut_mat[,rownames(WGS.tableA)], mySpectrum[,c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS39","SBSyst")], max_delta = 0.01, method = "best_subset")
fig_list <- strict_refit$sim_decay_fig
fit_res_strict <- strict_refit$fit_res
pdf(file.path( "./bootstrap/Stricter_contribution_Highdensity_round1.pdf"), width = 24, height = 7)
plot_contribution(fit_res_strict$contribution,coord_flip = FALSE, mode = "absolute")
dev.off()
plot_original_vs_reconstructed(mut_mat[,rownames(WGS.tableA)], fit_res_strict$reconstructed, y_intercept = 0.95)

### Bootstrapped refitting for 1000 times
contri_boots <- fit_to_signatures_bootstrapped(mut_mat[,rownames(WGS.tableA)],
  mySpectrum[,c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS39","SBSyst")],n_boots = 1000,max_delta = 0.01, method = "strict")
pdf(file.path( "./bootstrap/Bootstrapped_refitting_for_tableA.pdf"), width = 4, height = 120)
plot_bootstrapped_contribution(contri_boots)
dev.off()

ID <- colnames(mut_mat[,rownames(WGS.tableA)])
a <- contri_boots[which(grepl(paste0(ID[1],"_"),rownames(contri_boots) )),]
yy <- rep(0,7)
for(i in c(1:7)){
  if(length(which(a[,i]>0)) >= 950){
    yy[i] <- mean(a[,i]/rowSums(a))
  }
  else{
    yy[i] <- 0
  }
}
for(j in c(2:ncol(mut_mat[,rownames(WGS.tableA)]))){
  a <- contri_boots[which(grepl(paste0(ID[j],"_"),rownames(contri_boots) )),]
  xx <- rep(0,7)
  for(i in c(1:7)){
    if(length(which(a[,i]>0)) >= 950){
       xx[i] <- mean(a[,i]/rowSums(a))
    }
    else{
       xx[i] <- 0
    }
  }
  yy <- rbind(yy,xx)
}

rownames(yy) <- ID
colnames(yy) <- c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS39","SBSyst")

##############################3333
AA <- rep(0,7)
for(i in c(1:7)){
    AA[i] <- yy[1,i]/rowSums(yy)[1]
}
for(j in c(2:ncol(mut_mat[,rownames(WGS.tableA)]))){
  BB <- rep(0,7)
  for(i in c(1:7)){
    BB[i] <- yy[j,i]/rowSums(yy)[j]
  }
  AA <- rbind(AA,BB)
}

rownames(AA) <- ID
colnames(AA) <- c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS39","SBSyst")

AA<- as.data.frame(AA)
PartA <- data.frame(AA$SBS1,AA$SBS5,AA$SBS17a,AA$SBS17b,AA$SBS18,rep(0,nrow(AA)),AA$SBS39,AA$SBSyst)
colnames(PartA) <- c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS31","SBS39","SBSyst")
rownames(PartA) <- rownames(AA)
  
Final <- rbind(PartA,PartB)
write.table(Final,file = "./bootstrap/Final_contribution.txt",quote = F,sep = "\t",row.names = T)


## Similarity between mutational profiles and signatures
#mySpectrum2 <- mySpectrum[,c("SBS39","SBS18","SBSyst")]
#cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, mySpectrum2)
#cos_sim_samples_signatures[1:3, 1:3]
#
#plot_cosine_heatmap(cos_sim_samples_signatures, cluster_rows = TRUE, cluster_cols = TRUE)
#
#
#cos_sim_samples <- cos_sim_matrix(mut_mat, mut_mat)
#plot_cosine_heatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE)




