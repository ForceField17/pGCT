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

mutName <- read.table("../StandardName.txt",header = T)

library(sigfit)
data("cosmic_signatures_v3.3")


################################### purify SBSyst via simulation data
YSTsig = read.table("../WGS/April6/SBS96/All_Solutions/SBS96_7_Signatures/Signatures/SBS96_S7_Signatures.txt", header=T )
rownames(YSTsig)<- YSTsig$MutationType
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigA = YSTsig[Order$Type1,]
xx <- as.data.frame(YSTsigA$SBS96G)
rownames(xx) <- Order$Type2
colnames(xx) <- "sig"

probs <- c(1) %*% as.matrix(t(xx))
Simu_mutations <- matrix(rmultinom(1, 60000, probs), nrow = 1)
colnames(Simu_mutations) <- rownames(xx)

pdf(file.path( "./SigProfiler/SimulatedMut_for_denovoSignature.pdf"), width = 12, height = 6)
plot_spectrum(Simu_mutations, name = "Simulated counts")
dev.off()



############################### Fit-Ext models to check if the SBSyst can be further decomposed by COSMIC
known_signatures <- Cosmic3.3
mcmc_samples_extr_1 <- fit_extract_signatures(counts = Simu_mutations,
                                              signatures = known_signatures,
                                              num_extra_sigs=1)
names(mcmc_samples_extr_1)
signatures <- retrieve_pars(mcmc_samples_extr_1,par = "signatures")
rownames(signatures$mean)

pdf(file.path( "./SigProfiler/SBSA_FitExt_signature.pdf"), width = 12, height = 10)
par(mfrow = c(2, 1))
plot_spectrum(signatures)
dev.off()

exposures <- retrieve_pars(mcmc_samples_extr_1 , "exposures")

pdf(file.path( "./SigProfiler/exposure_round1.pdf"), width = 28, height = 10)
plot_exposures(mcmc_samples = mcmc_samples_extr_1)
dev.off()


############################### Fit-Ext models for SBS1
known_signatures <- Cosmic3.3[c("SBS1","SBS17a","SBS17b"),]
mcmc_samples_extr_3 <- fit_extract_signatures(counts = Simu_mutations,
                                              signatures = known_signatures,
                                              num_extra_sigs=1)
names(mcmc_samples_extr_3)
signatures <- retrieve_pars(mcmc_samples_extr_3,par = "signatures")
rownames(signatures$mean)

pdf(file.path( "./SigProfiler/SBSA_FitExt_signature_Final.pdf"), width = 12, height = 12)
par(mfrow = c(4, 1))
plot_spectrum(signatures)
dev.off()

myDenovo <- t(signatures$mean)
write.table(myDenovo,file = "SigProfiler/FitExt_denove.txt",sep = "\t",quote = F,row.names = T)

exposures <- retrieve_pars(mcmc_samples_extr_3 , "exposures")

pdf(file.path( "./SigProfiler/exposure_round2.pdf"), width = 12, height = 10)
plot_exposures(mcmc_samples = mcmc_samples_extr_3)
dev.off()






################################### purify SBSyst via simulation data
YSTsig = read.table("../WGS/DenoveSignaturePalimpsest.txt", header=T )
rownames(YSTsig)<- YSTsig$MutationType
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigC = YSTsig[Order$Type1,]
xx <- as.data.frame(YSTsigC$SBS96F)
rownames(xx) <- Order$Type2
colnames(xx) <- "sig"

probs <- c(1) %*% as.matrix(t(xx))
Simu_mutations <- matrix(rmultinom(1, 60000, probs), nrow = 1)
colnames(Simu_mutations) <- rownames(xx)

pdf(file.path( "./Palimpsest/SimulatedMut_for_denovoSignature.pdf"), width = 12, height = 6)
plot_spectrum(Simu_mutations, name = "Simulated counts")
dev.off()

############################### Fit-Ext models to check if the SBSyst can be further decomposed by COSMIC
known_signatures <- Cosmic3.3
mcmc_samples_extr_1 <- fit_extract_signatures(counts = Simu_mutations,
                                              signatures = known_signatures,
                                              num_extra_sigs=1)
names(mcmc_samples_extr_1)
signatures <- retrieve_pars(mcmc_samples_extr_1,par = "signatures")
rownames(signatures$mean)

pdf(file.path( "./Palimpsest/SBSA_FitExt_signature.pdf"), width = 12, height = 10)
par(mfrow = c(2, 1))
plot_spectrum(signatures)
dev.off()

exposures <- retrieve_pars(mcmc_samples_extr_1 , "exposures")

pdf(file.path( "./Palimpsest/exposure_round1.pdf"), width = 28, height = 10)
plot_exposures(mcmc_samples = mcmc_samples_extr_1)
dev.off()


############################### Fit-Ext models for SBS1
known_signatures <- Cosmic3.3[c("SBS1"),]
mcmc_samples_extr_3 <- fit_extract_signatures(counts = Simu_mutations,
                                              signatures = known_signatures,
                                              num_extra_sigs=1)
names(mcmc_samples_extr_3)
signatures <- retrieve_pars(mcmc_samples_extr_3,par = "signatures")
rownames(signatures$mean)

pdf(file.path( "./Palimpsest/SBSA_FitExt_signature_Final.pdf"), width = 12, height = 12)
par(mfrow = c(4, 1))
plot_spectrum(signatures)
dev.off()

myDenovo <- t(signatures$mean)
write.table(myDenovo,file = "Palimpsest/FitExt_denove.txt",sep = "\t",quote = F,row.names = T)

exposures <- retrieve_pars(mcmc_samples_extr_3 , "exposures")

pdf(file.path( "./Palimpsest/exposure_round2.pdf"), width = 12, height = 10)
plot_exposures(mcmc_samples = mcmc_samples_extr_3)
dev.off()






################################### purify SBSyst via simulation data
YSTsig = read.table("../WGS/SigFit/original_denovo.txt", header=T )
rownames(YSTsig)<- YSTsig$MutationType
Order <- mutName[order(mutName$Type3),] # put in standard order
YSTsigB = YSTsig[Order$Type2,]
xx <- as.data.frame(YSTsigB$SignatureA)
rownames(xx) <- Order$Type2
colnames(xx) <- "sig"

probs <- c(1) %*% as.matrix(t(xx))
Simu_mutations <- matrix(rmultinom(1, 60000, probs), nrow = 1)
colnames(Simu_mutations) <- rownames(xx)

pdf(file.path( "./SigFit/SimulatedMut_for_denovoSignature.pdf"), width = 12, height = 6)
plot_spectrum(Simu_mutations, name = "Simulated counts")
dev.off()

############################### Fit-Ext models to check if the SBSyst can be further decomposed by COSMIC
known_signatures <- Cosmic3.3
mcmc_samples_extr_1 <- fit_extract_signatures(counts = Simu_mutations,
                                              signatures = known_signatures,
                                              num_extra_sigs=1)
names(mcmc_samples_extr_1)
signatures <- retrieve_pars(mcmc_samples_extr_1,par = "signatures")
rownames(signatures$mean)

pdf(file.path( "./SigFit/SBSA_FitExt_signature.pdf"), width = 12, height = 10)
par(mfrow = c(2, 1))
plot_spectrum(signatures)
dev.off()

exposures <- retrieve_pars(mcmc_samples_extr_1 , "exposures")

pdf(file.path( "./SigFit/exposure_round1.pdf"), width = 28, height = 10)
plot_exposures(mcmc_samples = mcmc_samples_extr_1)
dev.off()


############################### Fit-Ext models for SBS1
known_signatures <- Cosmic3.3[c("SBS1"),]
mcmc_samples_extr_3 <- fit_extract_signatures(counts = Simu_mutations,
                                              signatures = known_signatures,
                                              num_extra_sigs=1)
names(mcmc_samples_extr_3)
signatures <- retrieve_pars(mcmc_samples_extr_3,par = "signatures")
rownames(signatures$mean)

pdf(file.path( "./SigFit/SBSA_FitExt_signature_Final.pdf"), width = 12, height = 12)
par(mfrow = c(4, 1))
plot_spectrum(signatures)
dev.off()

myDenovo <- t(signatures$mean)
write.table(myDenovo,file = "SigFit/FitExt_denove.txt",sep = "\t",quote = F,row.names = T)

exposures <- retrieve_pars(mcmc_samples_extr_3 , "exposures")

pdf(file.path( "./SigFit/exposure_round2.pdf"), width = 12, height = 10)
plot_exposures(mcmc_samples = mcmc_samples_extr_3)
dev.off()






