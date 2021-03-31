#' lmTestContrast
#'
#' @description assess a contrast for a given linear model as fit by lm
#'
#' @param fm - fit linear model 
#'
#' @param contrast - a vector indicating a contrast between model coefficients.
#' See the `getContrast` function to create a contrast.
#'
#' @param df_prior - prior degrees of freedom for moderated comparison; see also
#' `limma::squeezeVar` to compute prior degrees of freedom.
#'
#' @param s2_prior - prior sigma squared for moderated comparison; see also 
#'
#' @return a data.table giving the result of the model-based comparison,
#' contains the following columns:
#' \itemize{
#'   \item{Contrast - }{A vector defining a comparison between positive
#'   and negative coefficients in a model. See also, `getContrast`.}
#'   \item{log2FC - }{The fold change estimated from the fit model.}
#'   \item{percentControl}{The fold-change converted to a percentage of the
#'   control.}
#'   \item{SE - }{The standard error of the comparison computed as the square
#'   root of the variance.}
#'   \item{Tstatistic - }{The t-statistic for the comparison, calculated using
#'   equation 11 from Kuznetsova et al., 2017.}
#'   \item{Pvalue - }{The p-value for the comparison, calculated from the t-value
#'   and degrees of freedom using the student's t-distribution, `pt`.}
#'   \item{DF - }{The degrees of freedom for the comparison.}
#'   \item{S2 - }{sigma squared -- the estimated standard deviation of the errors
#'   or residual standard deviation (sigma) squared -- see also `?sigma`.}
#'
#' @export lmTestContrast

lmTestContrast <- function(fm, LT, 
			   s2_prior = 0, 
			   df_prior = 0) {

  ## !/usr/bin/env Rscript
  #
  ## library(SwipProteomics)
  #
  # data(swip)
  # data(msstats_prot)
  #
  # library(dplyr)
  #
  # fx <- "Abundance ~ Condition"
  # fm <- lm(fx, msstats_prot %>% 
  #            subset(Protein == swip))
  #
  # LT <- getContrast(fm,"ConditionMutant.F6",
  #                      "ConditionControl.F6")
  #
  # lmTestContrast(fm,LT) %>% knitr::kable()

  # check input
  stopifnot(inherits(fm, "lm"))
  stopifnot(all(names(LT) %in% names(coef(fm))))

  # create comparison
  pos_coef <- names(LT)[LT > 0]
  neg_coef <- names(LT)[LT < 0]
  comparison <- paste(pos_coef, neg_coef, sep = "-")

  # extract model coefficients
  coeff <- fm$coefficients

  # compute FC
  FC <- as.numeric(LT %*% coeff)
  coeff <- fm$coefficients

  # compute S2 and DF
  av <- anova(fm)
  s2 <- av["Residuals", "Mean Sq"]
  s2_df <- av["Residuals", "Df"]

  # compute posterior S2
  s2_post <- (s2_prior * df_prior + s2 * s2_df) / (df_prior + s2_df)

  # compute variance
  variance <- diag(t(LT) %*% summary(fm)$cov.unscaled %*% LT) * s2_post

  # compute posterior DF
  df_post <- s2_df + df_prior

  # compute t-statistic
  t <- FC / sqrt(variance)

  # compute p-value
  p <- 2 * pt(-abs(t), df = df_post)

  # collect stats
  contrast_stats <- data.frame(
    Contrast = comparison,
    log2FC = FC,
    percentControl = 100*(2^FC),
    SE = sqrt(variance),
    Tstatistic = t,
    Pvalue = p,
    DF = df_post
  )

  # return the statistics
  return(contrast_stats)
}
