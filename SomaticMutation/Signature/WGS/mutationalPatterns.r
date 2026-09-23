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


# overall Mutation spectrum
p1 <- plot_spectrum(type_occurrences);p1
p2 <- plot_spectrum(type_occurrences, CT = TRUE);p2
p3 <- plot_spectrum(type_occurrences, CT = TRUE, indv_points = TRUE, legend = FALSE);p3

grid.arrange(p1, p2, p3, ncol=3, widths=c(3,3,1.75))

  # Mutation spectrum by subtypes
p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE);p4
palette <- c("pink", "orange", "blue", "lightblue", "green", "red", "purple")
p5 <- plot_spectrum(type_occurrences, CT=TRUE, legend=TRUE, colors=palette, error_bars = "stdev");p5
grid.arrange(p4, p5, ncol=2, widths=c(4,2.3))


 ### 96 trinucleodide mutation profile 4 * 6 * 4
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)

## plot first 9 sampkes
plot_96_profile(mut_mat[,c(1:7)])
plot_96_profile(mut_mat[,c(20:30)], condensed = T,ymax = 0.06)

plot_compare_profiles(mut_mat[,1],mut_mat[,2],profile_names = colnames(mut_mat[,c(1,2)]),
                      profile_ymax = 0.08, diff_ylim = c(-0.08, 0.08),condensed = TRUE)

mut_mat_ext_context <- mut_matrix(grl, ref_genome, extension = 2)
head(mut_mat_ext_context)
plot_profile_heatmap(mut_mat_ext_context, by = tissue)
plot_river(mut_mat_ext_context[,c(1,4)])


### 4. Find optimal contribution of known signatures
mutName <- read.table("../StandardName.txt",header = T)
the <- read.table("../Purify_denovoSig/mySpectrum_finalized.txt",header = T,row.names = "MutationType")
Order <- mutName[order(mutName$Type3),] # put in standard order
mySpectrum <- as.matrix(the[Order$Type1,])
#rownames(mySpectrum) <- Order$Type1

signatures <- get_known_signatures()
fit_res <- fit_to_signatures(mut_mat, mySpectrum)

plot_contribution(fit_res$contribution,coord_flip = FALSE,mode = "absolute")

plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, y_intercept = 0.95)


### Stricter refitting
strict_refit <- fit_to_signatures_strict(mut_mat, mySpectrum, max_delta = 0.01, method = "best_subset")
fig_list <- strict_refit$sim_decay_fig
fig_list[[1]]

fit_res_strict <- strict_refit$fit_res
pdf(file.path( "../WGS/MutationPatterns/Stricter_contribution.pdf"), width = 24, height = 7)
plot_contribution(fit_res_strict$contribution,coord_flip = FALSE, mode = "absolute")
dev.off()
plot_original_vs_reconstructed(mut_mat, fit_res_strict$reconstructed, y_intercept = 0.95)


#merged_signatures <- merge_signatures(mySpectrum, cos_sim_cutoff = 0.8)


### Bootstrapped refitting.
contri_boots <- fit_to_signatures_bootstrapped(mut_mat[,],mySpectrum,n_boots = 1000,max_delta = 0.01, method = "strict")
pdf(file.path( "../WGS/MutationPatterns/Bootstrapped_refitting.pdf"), width = 4, height = 120)
plot_bootstrapped_contribution(contri_boots)
dev.off()

ID <- colnames(mut_mat)
a <- contri_boots[which(grepl(paste0(ID[1],"_"),rownames(contri_boots) )),]
yy <- rep(0,8)
for(i in c(1:8)){
  if(length(which(a[,i]>5)) >= 950){
    yy[i] <- mean(a[,i]/rowSums(a))
  }
  else{
    yy[i] <- 0
  }
}
for(j in c(2:ncol(mut_mat))){
  a <- contri_boots[which(grepl(paste0(ID[j],"_"),rownames(contri_boots) )),]
  xx <- rep(0,8)
  for(i in c(1:8)){
    if(length(which(a[,i]>5)) >= 950){
       xx[i] <- mean(a[,i]/rowSums(a))
    }
    else{
       xx[i] <- 0
    }
  }
  yy <- rbind(yy,xx)
}

rownames(yy) <- ID
colnames(yy) <- colnames(mySpectrum)


##############################3333
AA <- rep(0,8)
for(i in c(1:8)){
    AA[i] <- yy[1,i]/rowSums(yy)[1]
}
for(j in c(2:ncol(mut_mat))){
  BB <- rep(0,8)
  for(i in c(1:8)){
    BB[i] <- yy[j,i]/rowSums(yy)[j]
  }
  AA <- rbind(AA,BB)
}

rownames(AA) <- ID
colnames(AA) <- colnames(mySpectrum)



write.table(AA,file = "./bootstrap_contribution.txt",quote = F,sep = "\t",row.names = T)


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




