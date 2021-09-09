#' lmerTestContrast
#'
#' @description evaluate model-based comparison given a linear mixed-model as
#' fit by \insertRef{Kuznetsova2017}{lmerTest}
#'
#' @param fm - linear mixed model fit by `lmerTest::lmer`.
#'
#' @param contrast - a vector indicating a contrast between model coefficients.
#' See the `getContrast` function to create a contrast.
#'
#' @param df_prior - prior degrees of freedom for moderated comparison; see also
#' `limma::squeezeVar` to compute prior degrees of freedom.
#'
#' @param s2_prior - prior sigma squared for moderated comparison; see also 
#' `limma::squeezeVar` to compute prior sigma squared.
#'
#' @return a `data.frame` containing the result for the comparison and the
#' following columns:
#'
#' @return `Contrast` - indicates the comparison between model coefficients.
#'
#' @return `log2FC` - the fold change estimated from the fit model.
#'
#' @return `percentControl` - the fold-change converted to a percent relative to
#'   the control (positive) coefficient.
#'
#' @return `SE` - the standard error of the comparison computed as the square
#'   root of the variance.
#'
#' @return `Tstatistic` - the t-statistic for the comparison, calculated using
#'   equation (11) from Kuznetsova et al., 2017.
#'
#' @return `Pvalue` - the p-value for the comparison, calculated from the t-value
#'   and degrees of freedom using the student's t-distribution, `pt`.
#'
#' @return `DF` - the degrees of freedom for the comparison.
#'
#' @return `S2` - sigma squared -- the estimated standard deviation of the errors
#'   or residual standard deviation (sigma) squared -- see also `?sigma`.
#'
#' @return `isSingular` - see also `lme4::isSingular`. If `TRUE` this indicates
#' that one or more of the model's parameters explains little to no variance.
#'
#' @import lmerTest
#'
#' @importFrom dplyr %>%
#'
#' @export lmerTestContrast

lmerTestContrast <- function(fm, contrast,
                             df_prior = 0, s2_prior = 0) {

#      # example:
#      library(dplyr)
#      library(tidyProt)
#      library(SwipProteomics)
#  
#      data(swip) # "Q3UMB9"
#      data(swip_tmt) # Courtland et al., 2020
#  
#      fx <- formula("Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)")
#      fm <- lmerTest::lmer(fx, swip_tmt %>% filter(Protein == swip))
#  
#      # create a contrast!
#      geno <- swip_tmt$Genotype
#      biof <- swip_tmt$BioFraction
#      conditions <- levels(interaction(geno, biof))
#      contrast <- vector("numeric", length(conditions)) %>%
#   	    setNames(nm = names(lme4::fixef(fm)))
#     contrast[grepl(".*Mutant.*F4",names(contrast))] <- -1
#     contrast[grepl(".*Control.*F4",names(contrast))] <- +1
#  
#     # test the comparison!
#     lmerTestContrast(fm, contrast)

  stopifnot(all(names(contrast) == names(lme4::fixef(fm))))

  pos_group <- names(contrast[contrast > 0])
  neg_group <- names(contrast[contrast < 0])
  comparison <- paste(pos_group, neg_group, sep = "-")

  # compute Satterthwaite degrees of freedom
  model_summary <- summary(fm, ddf = "Satterthwaite")

  # variance associated with mixed effects
  #mixef_var <- as.data.frame(lme4::VarCorr(fm, comp = "Variance"))

  # collect model's degrees of freedom and coefficients (beta)
  s2_df <- as.numeric(model_summary$coefficients[,"df"][pos_group])[1]
  coeff <- model_summary$coeff[, "Estimate"] # == lme4::fixef(fm)

  # compute the unscaled covar matrix
  # all(unscaled_vcov * sigma^2 == vcov)
  unscaled_vcov <- fm@pp$unsc()

  # compute scaled variance-covariance matrix
  vcov <- as.matrix(model_summary$vcov) # == vcov(fm) == fm@vcov_beta

  # compute standard error^2
  se2 <- as.numeric(contrast %*% vcov %*% contrast) # == variance

  # extract asymtoptic var-covar matrix from fit model
  A <- fm@vcov_varpar

  # calculate gradient from gradient matrices
  # consider using lmerTest::calcSatterth(tv, L)
  g <- sapply(fm@Jac_list, function(gm) contrast %*% gm %*% contrast)

  # given gradient and asymptoptic var-covar, compute posterior df
  denom <- as.numeric(g %*% A %*% g)
  df_post <- 2 * se2^2 / denom + df_prior

  # calculate posterior s2
  s2 <- sigma(fm)^2
  s2_post <- (s2_prior * df_prior + s2 * s2_df) / (df_prior + s2_df)

  # compute variance given s2_post
  vcov_post <- unscaled_vcov * s2_post # same as vcov
  variance <- as.numeric(contrast %*% vcov_post %*% contrast)

  # compute fold change and the t-statistic [lmerTest eq 11]
  FC <- (contrast %*% coeff)[, 1]
  t <- FC / sqrt(variance)

  # compute the p-value given t-statistic and posterior degrees of freedom
  p <- 2 * pt(-abs(t), df = df_post)

  # collect stats
  prot_stats <- data.frame(
    Contrast = comparison,
    log2FC = FC,
    percentControl = 100*(2^FC),
    SE = sqrt(variance),
    Tstatistic = t,
    Pvalue = p,
    DF = df_post,
    S2 = s2,
    isSingular = lme4::isSingular(fm)
  )

  return(prot_stats)
} # EOF
