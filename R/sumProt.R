#' sumProt
#'
#' summarize to protein level
#'
#' @export sumProt
#'
#' @importFrom dplyr %>% ungroup group_by summarize mutate

sumProt <- function(tp) {
  # Sum to protein level
  # any NA become 0
  tp$Intensity[is.na(tp$Intensity)] <- 0
  tp <- dplyr::ungroup(tp)
  proteins <- tp %>%
    dplyr::group_by(
      Mixture, Condition, BioFraction, Accession
    ) %>%
    dplyr::summarize(
      Peptides = length(Intensity),
      Intensity = sum(Intensity, na.rm = TRUE), .groups="drop"
    )
  proteins$Intensity[proteins$Intensity == 0] <- NA
  return(proteins)
}
