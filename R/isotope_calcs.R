# Here we will include functions for making isotope corrections for del 13C
# Manipulate DisplayDelta1 column

# Might want to use the df with nmol/g calculated in order to take weighted average
# of biomarkers

#' Process del13C values
#'
#' This function corrects the raw del13C output for methylation and normalizes
#' to the del13C standard USGS40 (L-glutamic acid). This replaces
#' \code{correct_iso()} as of 8/10/2021.
#'
#' @param \code{df} peak list dataframe or tibble. This should contain the
#' following columns at a minimum: Batch, DataFileName, RetTimeSecs,
#' MajorHeightnA, TotalPeakArea1, DisplayDelta1, Name. This can be the output
#' from \code{import_batch()} or \code{import_batch_multi()}.
#'
#' @param \code{c13_correction} difference between measured USGS40 value and
#' expected USGS40 value (-26.39 per mil).
#'
#' @param \code{methanol_13c} 13C signature of the methanol used for
#' methylation.
#'
#' @return Returns a tibble with columns from input dataframe and a corrected
#' d13C column.
#'
#' @examples
#'
correct_iso_tidy <- function(df, d13c_correction){
  corrected_df <- df %>%
    left_join(lipid_reference, by = c('Name' = 'fame')) %>%
    mutate(d13C_corrected = ((DisplayDelta1 - !!d13c_correction)*(1) -
                               d13CMethanol))

}

