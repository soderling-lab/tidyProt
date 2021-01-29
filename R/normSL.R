#' normSL
#'
#' perform sample loading normalization
#'
#' @export normSL
#'
#' @importFrom dplyr mutate group_by %>% summarize

normSL <- function(tp) {
  sample_sums <- tp %>% dplyr::group_by(Sample) %>% 
	  # sum(Intensity) for all Samples
	  dplyr::summarize(Total = sum(Intensity, na.rm=TRUE), 
			   .groups="drop") %>%
          # calculate normalization factors
          dplyr::mutate(mean = mean(Total)) %>%
	  dplyr::mutate(normFactor = mean / Total) %>%
	  dplyr::mutate(norm = Total * normFactor)
  # collect normalization factors as named vector
  r <- sample_sums$normFactor
  names(r) <- as.character(sample_sums[["Sample"]])
  # normalize measurements by sample factors
  SL_tp <- tp %>% 
	  dplyr::mutate(Intensity = Intensity * r[as.character(Sample)])
  return(SL_tp)
}
