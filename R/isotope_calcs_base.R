# Here we will include functions for making isotope corrections for del 13C
# Manipulate DisplayDelta1 column

# Might want to use the df with nmol/g calculated in order to take weighted average
# of biomarkers

#' Process del13C values
#'
#' This function corrects the raw del13C output for methylation and normalizes
#' to the del13C standard USGS40 (L-glutamic acid).
#'
#' @param \code{df} peak list dataframe or tibble. This should contain the
#' following columns at a minimum: Batch, DataFileName, RetTimeSecs,
#' MajorHeightnA, TotalPeakArea1, DisplayDelta1, Name. This can be the output
#' from import_batch() or import_batch_multi().
#'
#' @param \code{d13d_correction} difference between measured USGS40 value and
#' expected USGS40 value (-26.39 per mil). If IonVantage output was already
#' corrected, the enter 0.
#'
#' @param \code{methanol_13c} 13C signature of the methanol used for
#' methylation.
#'
#' @return Returns a tibble with columns from input dataframe and a corrected
#' d13C column.
#'
#' @examples
#'
correct_iso_base <- function(df, d13c_correction, methanol_13c){
  df[['DisplayDelta1']] <- as.numeric(df[['DisplayDelta1']])

  corrected_df <- merge(df, lipid_reference, by.x = 'Name', by.y = 'fame',
                        all.x = TRUE)
  corrected_df <- transform(corrected_df, d13C_corrected =
                              ((C_PLFA + 1) * (DisplayDelta1- d13c_correction) -
                                 methanol_13c) / C_PLFA)

}
# Does the methanol correction have to be corrected for too?
# Or is this number already include instrument correction?
