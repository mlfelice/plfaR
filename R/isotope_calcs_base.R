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
#' @param \code{c13_correction} difference between measured USGS40 value and
#' expected USGS40 value (-26.39 per mil). If IonVantage output has already
#' been corrected, then enter 0.
#'
#' @param \code{methanol_13c} 13C signature of the methanol used for
#' methylation.
#'
#' @return Returns a tibble with columns from input dataframe and a corrected
#' d13C column.
#'
#' @examples
#'
correct_iso_base <- function(df, d13c_correction, methanol_13c,
                             min_height = NULL){
  df[['DisplayDelta1']] <- as.numeric(df[['DisplayDelta1']])

  corrected_df <- merge(df, lipid_reference, by.x = 'Name', by.y = 'fame',
                        all.x = TRUE)
  corrected_df <- transform(corrected_df, d13C_corrected =
                              ((C_PLFA + 1) * (DisplayDelta1 - d13c_correction) -
                                 methanol_13c) / C_PLFA)

  corrected_df <- filter_by_height(corrected_df, min_height = min_height)

  return(corrected_df)

}

correct_iso_base2 <- function(df, d13c_correction, meth_13c_df,
                              min_height = NULL){
  # This version uses a dataframe matching methanol lot to batch so
  # that we can use different methanol lots for an experiment
  #
  # meth_13c_df should have 2 columns: Batch and Meth13C

  df <- merge(df, meth_13c_df, by = 'Batch', all.x = TRUE)

  df[['DisplayDelta1']] <- as.numeric(df[['DisplayDelta1']])

  corrected_df <- merge(df, lipid_reference, by.x = 'Name', by.y = 'fame',
                        all.x = TRUE)

  corrected_df <- transform(corrected_df, d13C_corrected =
                              ((C_PLFA + 1) * (DisplayDelta1 - d13c_correction) -
                                 Meth13C) / C_PLFA)

  corrected_df <- filter_by_height(corrected_df, min_height = min_height)

  return(corrected_df)

}

# Does the methanol correction have to be corrected for too?
# Or is this number already include instrument correction?


indic_iso_base <- function(df, sample_metadata){

  indicator_list <- list(f_lipids <- c('16:1 w5c', '18:1 w9c', '18:2 w6,9c'),
                         b_lipids <- c('13:0 iso', '13:0 anteiso', '14:0 3OH', '15:0 iso', '15:0 anteiso',
                                       '16:0 iso', '16:1 w7c', '16:0 10me', '17:0 iso', '17:0 anteiso',
                                       '18:1 w9t', '18:1 w7c', '18:0 10me'),
                         gram_pos_lipids <- c('13:0 iso', '13:0 anteiso', '15:0 iso', '15:0 anteiso',
                                              '17:iso', '17:0 anteiso'),
                         gram_neg_lipids <- c('14:0 3OH', '16:1 w7c', '16:1 w9c', '18:1 w7c',
                                              '18:1 w9t'),
                         actino_lipids <- c('16:0 10me', '18:0 10me'),
                         amf_lipids <- c('16:1 w5c'),
                         s_fungi_lipids <- c('18:1 w9c', '18:1 w6,9c'),
                         anaerobe_lipids <- c('19:0 cyclo')
  )
  # Naming list facilitates naming the output of calcs
  names(indicator_list) <- c('f_lipids', 'b_lipids', 'gram_pos_lipids', 'gram_neg_lipids',
                             'actino_lipids', 'amf_lipids', 's_fungi_lipids',
                             'anaerobe_lipids')

  wide_df <- reshape(data = df[c('BatchDataFileName', 'DataFileName', 'Batch',
                                 'Name', 'd13C_corrected')],
                     timevar = 'Name',
                     idvar = c('BatchDataFileName', 'DataFileName', 'Batch'),  # w/o 2+ idvars, some vals are equal and get dropped
                     direction = 'wide')
  colnames(wide_df)<- gsub('d13C_corrected.', '', colnames(wide_df))
  #wide_df[is.na(wide_df)] <- 0


  indic_df <- cbind(wide_df['BatchDataFileName'],
                    wide_df['DataFileName'],
                    lapply(indicator_list,
                           function(x){
                             tmp_df <- wide_df[, colnames(wide_df)[colnames(wide_df) %in% x]]
                             rowMeans(as.data.frame(tmp_df), na.rm = TRUE)  # Rows with all NA will return NaN
                           }
                    )
  )

  indic_df_long <- reshape(indic_df, varying = names(indic_df)[c(3:ncol(indic_df))],
                           v.names = 'avg_d13C_corrected',
                           timevar = 'Indicator', idvar = 'BatchDataFileName',
                           times = names(indic_df)[c(3:ncol(indic_df))], direction = 'long')

  indic_df_long <- merge(indic_df_long, sample_metadata, by = 'DataFileName', all.x = TRUE)

  indic_df_long[is.na(indic_df_long[['avg_d13C_corrected']]), 'avg_d13C_corrected'] <- NA  #replace NaN's resulting from rows containing all NA

  return(indic_df_long)
}

filter_by_height <- function(df, min_height){

  if (!is.null(min_height)){
    df_filtered <- df[df[['MajorHeightnA']] >= min_height, ]

    n_filtered <- nrow(df) - nrow(df_filtered)
  } else{

    n_filtered <- 0

    df_filtered <- df

  }

  cat(n_filtered, 'of', nrow(df), 'peaks removed from batch',
      unique(df[['Batch']])[1], '\n')

  return(df_filtered)

}
