# pGCT
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
library(ggplot2)
library("gridExtra")
library(reticulate)
#use_python("/Users/songdong/miniconda3/envs/spmg_r/bin/")
use_python("/Users/dongsong/miniconda3/envs/spmg_r/bin/")
### conda activate spmg_r  [run this commend in terminal]
library("SigProfilerMatrixGeneratorR")
library("SigProfilerPlottingR")
library("SigProfilerExtractorR")

#matrices_Control <- SigProfilerMatrixGeneratorR("test_SBS_Control", "GRCh38", "Control", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
matrices_GCT <- SigProfilerMatrixGeneratorR("test_SBS_Control", "GRCh38", "GCT", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)

#plotSBS(matrix_path = "GCT/output/SBS/test_SBS_Control.SBS96.all",project = "GCT" ,output_path="theGCT",plot_type = "96", percentage = F)
#plotSBS("GCT/output/SBS/test_SBS_Control.SBS96.all", output_path="theGCT", "GCT", "96", percentage=TRUE)


#matrices <- SigProfilerMatrixGeneratorR("test", "GRCh38", "Merge", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F )
#matrices_all <- SigProfilerMatrixGeneratorR("test", "GRCh38", "GCT", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F )

#plotSBS("Merge/output/SBS/test.SBS96.all", output_path="Comparison", "GCT", "96", percentage=TRUE)

sigprofilerextractor("matrix", "SBS2026","./GCT/output/SBS/test_SBS_Control.SBS96.all", reference_genome="GRCh38",
                     opportunity_genome = "GRCh38", context_type = "default",
                     exome = F, minimum_signatures=1, maximum_signatures=15,
                     nmf_replicates=200, resample = T, batch_size=1, cpu=-1,
                     gpu=F, nmf_init="random", precision= "single",
                     matrix_normalization= "gmm", seeds= "random",
                     min_nmf_iterations= 10000, max_nmf_iterations=2000000,
                     nmf_test_conv= 20000, nmf_tolerance= 1e-15,
                     nnls_add_penalty=0.05, nnls_remove_penalty=0.01,
                     initial_remove_penalty=0.05, get_all_signature_matrices= F)


sigprofilerextractor("matrix", "ID2026","./GCT/output/ID/test_SBS_Control.ID83.all", reference_genome="GRCh38",
                     opportunity_genome = "GRCh38", context_type = "default",
                     exome = F, minimum_signatures=1, maximum_signatures=10,
                     nmf_replicates=200, resample = T, batch_size=1, cpu=-1,
                     gpu=F, nmf_init="random", precision= "single",
                     matrix_normalization= "gmm", seeds= "random",
                     min_nmf_iterations= 1000, max_nmf_iterations=100000,
                     nmf_test_conv= 20000, nmf_tolerance= 1e-15,
                     nnls_add_penalty=0.05, nnls_remove_penalty=0.01,
                     initial_remove_penalty=0.05, get_all_signature_matrices= F)

sigprofilerextractor("matrix", "DBS2026","./GCT/output/DBS/test_SBS_Control.DBS78.all", reference_genome="GRCh38",
                     opportunity_genome = "GRCh38", context_type = "default",
                     exome = F, minimum_signatures=1, maximum_signatures=10,
                     nmf_replicates=200, resample = T, batch_size=1, cpu=-1,
                     gpu=F, nmf_init="random", precision= "single",
                     matrix_normalization= "gmm", seeds= "random",
                     min_nmf_iterations= 1000, max_nmf_iterations=100000,
                     nmf_test_conv= 20000, nmf_tolerance= 1e-15,
                     nnls_add_penalty=0.05, nnls_remove_penalty=0.01,
                     initial_remove_penalty=0.05, get_all_signature_matrices= F)



sigprofilerextractor("matrix", "April7","./GCT/output/SBS/test.SBS96.all", reference_genome="GRCh38",
                     min_nmf_iterations= 10, max_nmf_iterations=20,
                     opportunity_genome = "GRCh38", context_type = "default",cosmic_version=3.3,
                     minimum_signatures=7,maximum_signatures=7)



estimate_solution(base_csvfile="./April6/SBS96/All_solutions_stat.csv", 
                  All_solution="./April6/SBS96/All_Solutions/", 
                  genomes="./April6/SBS96/Samples.txt", 
                  output="./April6/results", 
                  title="Selection_Plot",
                  stability=0.8, 
                  min_stability=0.6, 
                  combined_stability=1.25)

library(SigProfilerAssignmentR)




cosmic_fit(samples="./GCT/output/SBS/test.SBS96.all",  output="./April6/MySelection",
           input_type="matrix",signature_database = "April7/SBS96/All_Solutions/SBS96_6_Signatures/Signatures/SBS96_S6_Signatures.txt",
           context_type="96",exome=FALSE,
           exclude_signature_subgroups=NULL, export_probabilities=T, 
           genome_build="GRCh38",export_probabilities_per_mutation=T, make_plots=TRUE,
           sample_reconstruction_plots="pdf")

#cosmic_fit <- function(samples,
#                       output,
#                       signatures=NULL,
#                       signature_database=NULL,
#                       nnls_add_penalty=0.05,
#                       nnls_remove_penalty=0.01,
#                       initial_remove_penalty=0.05,
#                       genome_build="GRCh37",
#                       cosmic_version=3.3,
#                       make_plots=T,
#                       collapse_to_SBS96=T,
#                       connected_sigs=T,
#                       verbose=F,
#                       devopts=NULL,
#                       exclude_signature_subgroups=NULL,
#                       exome=F,
#                       input_type='matrix',
#                       context_type="96",
#                       export_probabilities=T,
#                       export_probabilities_per_mutation=F,
#                       sample_reconstruction_plots=F
