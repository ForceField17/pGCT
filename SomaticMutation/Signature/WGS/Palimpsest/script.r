# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library("gridExtra")
library(ggplot2)
# Load Palimpsest & reference genome packages
library(Palimpsest)
library(BSgenome.Hsapiens.UCSC.hg38) # Use this package for hg38 data


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

Cosmic3.3 <- read.table("../../standard_cosmic_v3.3.txt",header = T,row.names = "Type")
Cosmic3.3 <- t(Cosmic3.3)

GCT <- read.table("./allGCT.vcf",header = T)
GCT <- annotate_VCF(vcf = GCT, ref_genome = BSgenome.Hsapiens.UCSC.hg38,  ref_fasta = "/Users/songdong/Dropbox/Dropbox/germ_cell_tumor/pGCT/Signature/hg38.fa", add_ID_cats = TRUE)

resdir <- file.path("denovo_allGCT/"); if(!file.exists(resdir))	dir.create(resdir) ## set results directory
# Extract de novo SBS signatures
SBS_input <- palimpsest_input(vcf = GCT, Type = "SBS")
SBS_denovo_sigs <- NMF_Extraction(input_matrices = SBS_input, range_of_sigs = 1:10,num_of_sigs = 7, nrun = 20,  method = "brunet", resdir = resdir)


# Compare the de novo signatures with published COSMIC signatures
compare_tab <- compare_results(reference_sigs = Cosmic3.3, extraction_1 = SBS_denovo_sigs); compare_tab
readr::write_delim(compare_tab, path = file.path(resdir,"Comparison_table.txt"))

pdf(file.path(resdir, "Cosine_Similarity_Heatmap.pdf"), width = 11, height = 10)
SBS_cosine_similarities <- deconvolution_compare(SBS_denovo_sigs, Cosmic3.3)
dev.off()
save(SBS_cosine_similarities, file = file.path(resdir, "Cosine_Similarity_matrix.RData"))


# Define signature colours for plotting
SBS_col <- signature_colour_generator(rownames(SBS_denovo_sigs))

# Calculate and plot the exposure of the signatures across the series (Plotting the contribution of SBS signatures in each sample)
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input, input_signatures = SBS_denovo_sigs,
                                        resdir = resdir, signature_colours = SBS_col, input_vcf = GCT)

pdf(file.path(resdir, "signature_content_plot.pdf"), width=12, height=5)
deconvolution_exposure(signature_colours = SBS_col , signature_contribution = SBS_signatures_exp)
dev.off()



#-------------------------------------------------------------------------------------------------
# 3] Extract with published SBS COSMIC signatures
#-------------------------------------------------------------------------------------------------
resdir <- file.path("Cosmic3.3_Extraction/"); if(!file.exists(resdir))	dir.create(resdir) 

# select desired COSMIC SBS reference signatures 
SBS_liver_names <- c("SBS1","SBS5","SBS17a","SBS17b","SBS18","SBS31") 

for(new_name in c(SBS_liver_names[SBS_liver_names %!in% names(sig_cols)]))
  sig_cols[new_name] <- signature_colour_generator(new_name)  ## generate colours for new signatures 




SBS_liver_sigs <- Cosmic3.3[rownames(Cosmic3.3) %in% SBS_liver_names,]
denovoYST <- as.data.frame(SBS_denovo_sigs[c(3,4),])
colnames(denovoYST) <- colnames(SBS_liver_sigs)
SBS_liver_sigs <- rbind(SBS_liver_sigs,denovoYST)
SBS_liver_sigs <- SBS_liver_sigs[which(rownames(SBS_liver_sigs)!="SBS_denovo_4"),]

write.table(file = "SBS_liver_sigs.txt",SBS_liver_sigs ,quote = F,sep="\t")
# calculate and plot the exposure of the signatures across the series
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input, input_signatures = SBS_liver_sigs,
                                        threshold = 6, signature_colours = gg_color_hue(8), resdir = resdir, 
                                        input_vcf = vcf)

pdf(file.path(resdir, "SBS_signature_content_plot.pdf"), width=12, height=5)
deconvolution_exposure(signature_contribution = SBS_signatures_exp, signature_colours = gg_color_hue(5))
dev.off()



