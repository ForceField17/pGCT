# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
# ==============================================================
# Disentangling 1p36 loss from L1HS derepression
# Per-cancer-type negative binomial regression of L1 insertion burden
#
# Question: after holding L1HS methylation and genomic instability
# constant, does 1p36 loss retain an independent association with
# somatic L1 insertion burden?
# ==============================================================

library(MASS)      # glm.nb
library(car)       # vif
library(dplyr)
library(tibble)

# --------------------------------------------------------------
# 1. INPUT
# --------------------------------------------------------------
# One row per TUMOR SAMPLE. Required columns:
#
#   sample_id    chr
#   cancer_type  chr      "HNSC", "LUSC", "UCEC", "TGCT"
#   insertions   int      somatic L1 insertion count
#   loss1p36     0/1      1 = 1p36 loss
#   meth_L1HS    num      PER-SAMPLE L1HS promoter methylation (beta)
#                         NOT the cancer-type mean
#   purity       num      tumor purity (optional, recommended)
#
# dat <- readr::read_tsv("l1_pancancer_sample_table.tsv")

TCGA <- read.table("../../Pancancer_L1/Final_CNV_RT.txt",header = T)
TCGA$insertions <- TCGA$RT_burden
Methy <- read.table("../TCGA_samples_with_methy_and_CNV.txt",header = T)
dat <- merge(TCGA,Methy)
dat$sample_id <- dat$patient_id
dat$cancer_type <- dat$subtype
dat$loss1p36 <- 0; dat$loss1p36[which(dat$chr1p36=="Loss")] <- 1
dat$meth_L1HS <- dat$MeanL1HSmethy

TYPES <- c("HNSC", "LUSC", "UCEC", "TGCT")


USE_PURITY <- FALSE   # set TRUE if a purity column is present

# --------------------------------------------------------------
# 2. SANITY CHECKS  -- read this table before trusting any model
# --------------------------------------------------------------
check_type <- function(d, type) {
  # Spearman between the two predictors: this IS the collinearity
  # question. With only two predictors, VIF = 1 / (1 - r^2).
  r <- cor(d$meth_L1HS, d$loss1p36, method = "spearman")
  tibble(
    cancer_type          = type,
    n                    = nrow(d),
    n_loss               = sum(d$loss1p36 == 1),
    n_neutral            = sum(d$loss1p36 == 0),
    pct_zero_insertions  = round(mean(d$insertions == 0) * 100, 1),
    median_ins_loss      = median(d$insertions[d$loss1p36 == 1]),
    median_ins_neutral   = median(d$insertions[d$loss1p36 == 0]),
    var_mean_ratio       = round(var(d$insertions) / mean(d$insertions), 1),
    spearman_meth_1p36   = round(r, 3),
    approx_vif           = round(1 / (1 - r^2), 2)
  )
}

sanity <- bind_rows(lapply(TYPES, function(t)
  check_type(filter(dat, cancer_type == t), t)))
print(as.data.frame(sanity))

# Interpretation:
#   n_loss or n_neutral < ~10   -> too unstable to model; report
#                                  descriptively instead
#   pct_zero_insertions > ~50   -> consider zero-inflated NB (sec. 7)
#   var_mean_ratio >> 1         -> negative binomial is the right family
#   approx_vif > 5              -> the two predictors are largely
#                                  redundant in this cancer type;
#                                  expect wide CIs (see sec. 4)

# --------------------------------------------------------------
# 3. MODEL FITTING
# --------------------------------------------------------------
build_formula <- function() {
  rhs <- c("loss1p36", "meth_z")
  if (USE_PURITY) rhs <- c(rhs, "purity")
  as.formula(paste("insertions ~", paste(rhs, collapse = " + ")))
}

fit_one <- function(d, type) {
  
  d <- d %>%
    mutate(
      loss1p36 = factor(loss1p36, levels = c(0, 1),
                        labels = c("neutral", "loss")),
      # z-score WITHIN cancer type: beta values span only ~0.74-0.88,
      # so an unscaled coefficient would be uninterpretably large,
      # and scaling makes coefficients comparable across types
      meth_z = as.numeric(scale(meth_L1HS))
    )
  
  m <- tryCatch(
    suppressWarnings(glm.nb(build_formula(), data = d)),
    error = function(e) {
      message(type, ": glm.nb failed -- ", e$message); NULL
    })
  if (is.null(m)) return(NULL)
  
  ci <- suppressMessages(suppressWarnings(confint(m)))  # profile likelihood
  
  tab <- tibble(
    cancer_type = type,
    n           = nrow(d),
    term        = names(coef(m)),
    rate_ratio  = exp(coef(m)),
    ci_low      = exp(ci[, 1]),
    ci_high     = exp(ci[, 2]),
    p           = coef(summary(m))[, "Pr(>|z|)"],
    theta       = m$theta,
    theta_se    = m$SE.theta
  ) %>% filter(term != "(Intercept)")
  
  list(model = m, table = tab)
}

fits <- lapply(TYPES, function(t)
  fit_one(filter(dat, cancer_type == t), t))
names(fits) <- TYPES

main_table <- bind_rows(lapply(fits, function(x) x$table))
print(as.data.frame(main_table), digits = 3)

# --------------------------------------------------------------
# 4. HOW TO READ THE loss1p36 ROW
# --------------------------------------------------------------
# rate_ratio is a RATE RATIO: 1.8 means that, at equal L1HS
# methylation, 1p36-loss tumors carry 1.8x the insertion rate.
#
#   CI excludes 1
#       -> 1p36 retains an association independent of methylation
#
#   CI spans 1 and is NARROW (e.g. 0.90-1.15)
#       -> genuinely no meaningful residual effect. Informative null.
#
#   CI spans 1 and is WIDE (e.g. 0.4-3.5), approx_vif high
#       -> the data cannot separate the two effects.
#          Do NOT write "no effect". Write that the contributions
#          could not be estimated independently because copy-number
#          status and L1HS methylation are collinear in this cohort,
#          and quote the VIF and CI as evidence.

# --------------------------------------------------------------
# 5. NESTED MODELS + LIKELIHOOD RATIO TESTS   -- the crux
# --------------------------------------------------------------
compare_models <- function(d, type) {
  d <- d %>% mutate(
    loss1p36 = factor(loss1p36, levels = c(0, 1)),
    meth_z   = as.numeric(scale(meth_L1HS)))
  
  m_null <- glm.nb(insertions ~ 1,                  data = d)
  m_meth <- glm.nb(insertions ~ meth_z,             data = d)
  m_1p36 <- glm.nb(insertions ~ loss1p36,           data = d)
  m_full <- glm.nb(insertions ~ loss1p36 + meth_z,  data = d)
  
  tibble(
    cancer_type = type,
    p_meth_alone      = anova(m_null, m_meth)$`Pr(Chi)`[2],
    p_1p36_alone      = anova(m_null, m_1p36)$`Pr(Chi)`[2],
    # does 1p36 add anything ON TOP OF methylation?
    p_1p36_given_meth = anova(m_meth, m_full)$`Pr(Chi)`[2],
    # does methylation add anything ON TOP OF 1p36?
    p_meth_given_1p36 = anova(m_1p36, m_full)$`Pr(Chi)`[2],
    aic_meth_only = AIC(m_meth),
    aic_1p36_only = AIC(m_1p36),
    aic_full      = AIC(m_full)
  )
}

lrt <- bind_rows(lapply(TYPES, function(t)
  tryCatch(suppressWarnings(
    compare_models(filter(dat, cancer_type == t), t)),
    error = function(e) {message(t, ": LRT failed"); NULL})))
print(as.data.frame(lrt), digits = 3)

# Decision table:
#   p_1p36_given_meth NS  +  p_meth_given_1p36 sig
#       -> methylation is doing the work in this cancer type
#   both significant
#       -> independent contributions
#   both NS, but p_1p36_alone AND p_meth_alone both significant
#       -> textbook collinearity: each explains the same variance.
#          Report as inseparable, not as absent.

# --------------------------------------------------------------
# 6. DESCRIPTIVE BACKUP  (always report alongside the models)
# --------------------------------------------------------------
# Cliff's delta of insertion burden by 1p36 status, per type.
# Ties are common with zero-inflated counts and shrink |delta|
# toward 0, so treat this as supporting, not primary.
cliffs_delta <- function(x, y) {
  n <- outer(x, y, ">") ; m <- outer(x, y, "<")
  (sum(n) - sum(m)) / (length(x) * length(y))
}
desc <- bind_rows(lapply(TYPES, function(t) {
  d <- filter(dat, cancer_type == t)
  tibble(cancer_type = t,
         cliffs_delta = round(cliffs_delta(
           d$insertions[d$loss1p36 == 1],
           d$insertions[d$loss1p36 == 0]), 3),
         wilcox_p = wilcox.test(insertions ~ loss1p36, data = d)$p.value)
}))
print(as.data.frame(desc), digits = 3)

# --------------------------------------------------------------
# 7. IF glm.nb MISBEHAVES
# --------------------------------------------------------------
# Warnings about theta usually mean many zeros or small n.
# library(glmmTMB)
# glmmTMB(insertions ~ loss1p36 + meth_z, ziformula = ~1,
#         family = nbinom2, data = d)

# --------------------------------------------------------------
# 8. TO ADD BEFORE REVISION: genomic instability burden
# --------------------------------------------------------------
# Sources to check for per-sample values:
#   - cBioPortal: "Fraction Genome Altered" is provided as a clinical
#     attribute for TCGA studies and can be downloaded per sample
#   - Taylor et al. 2018 (TCGA aneuploidy paper) published arm-level
#     copy-number calls and aneuploidy scores; because calls are
#     arm-level you can RECOMPUTE the score with chr1p EXCLUDED
#
# *** Whatever measure you use, exclude 1p (or at minimum 1p36).
#     Otherwise you are partially adjusting 1p36 for itself, which
#     absorbs part of the effect you are trying to estimate. ***
#
# Then set:
#   dat$log_scna_z <- as.numeric(scale(log(dat$scna_burden + 1e-3)))
# and add "log_scna_z" to build_formula() and to every model in
# compare_models().

# --------------------------------------------------------------
# 9. EXPORT
# --------------------------------------------------------------
write.csv(sanity,     "S16_sanity_checks.csv",     row.names = FALSE)
write.csv(main_table, "S16_nb_coefficients.csv",   row.names = FALSE)
write.csv(lrt,        "S16_model_comparison.csv",  row.names = FALSE)
write.csv(desc,       "S16_descriptive.csv",       row.names = FALSE)
