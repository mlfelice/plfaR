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
#' from \code{import_batch()} or \code{import_batch_multi()}.
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

# TODO: allow for use of single vmethanol d13C alue or dataframe of methanol
# d13C values

#' Process del13C values
#'
#' This function corrects the raw del13C output for methylation and normalizes
#' to the del13C standard USGS40 (L-glutamic acid). This function is an update
#' of the correct_iso_base() function that uses a dataframe matching methanol
#' lot to batch so that we can use different methanol lots for an experiment.
#' Replaces \code{correct_iso_base2()} as of 8/10/2021.
#'
#' @param \code{df} peak list dataframe or tibble. This should contain the
#' following columns at a minimum: Batch, DataFileName, RetTimeSecs,
#' MajorHeightnA, TotalPeakArea1, DisplayDelta1, Name. This can be the output
#' from import_batch() or import_batch_multi().
#'
#' @param \code{d13c_correction} difference between measured USGS40 value and
#' expected USGS40 value (-26.39 per mil). If IonVantage output has already
#' been corrected, then enter 0.
#'
#' @param \code{meth_13c_df} Dataframe identifying the 13C signature of the
#' methanol used for each batch. This dataframe should have 2 columns: Batch
#' and Meth13C.
#'
#' @return Returns a tibble with columns from input dataframe and a corrected
#' d13C column.
#'
#' @examples
#'
#'
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


#' Calculate d13C for indicator groups
#'
#' This function aggregates individual lipid d13C signatures by indicator
#' group. Current implementation uses a straight average. Future versions
#' may offer the option to select between average and weighted average.
#'
#' @param \code{df} dataframe containing corrected d13C values for individual
#' lipids. Normally, this would be the output of correct_iso_base().
#' This should contain the following columns at a minimum: Batch, DataFileName,
#' BatchDataFileName, d13C_corrected, Name.
#'
#' @param \code{sample metadata}
#'
#' @return Returns a tibble with corrected d13C values for each indicator
#' group.
#'
#' @examples
#'
#'

# Does the methanol correction have to be corrected for too?
# Or is this number already include instrument correction?

indic_iso_base <- function(df, sample_metadata){

  indicator_list <- list(f_lipids = c('18:1 w9c', '18:2 w6,9c'),
                         # Lipids representing total bacteria for f to b ratio not shown in Cameron's data, so for now, I will use lipids for Gram +/-, actino, and anaerobe
                         b_lipids = c('15:0 iso', '15:0 anteiso',
                                      '16:1 w7c', '16:0 10me',
                                      '18:1 w9c', '18:1 w9t',
                                      '18:2 w6,9c', '18:0 10me', '19:0 cyclo'),
                         gram_pos_lipids = c('15:0 iso', '15:0 anteiso'),
                         gram_neg_lipids = c('16:1 w7c', '18:1 w9t'),
                         actino_lipids = c('16:0 10me', '18:0 10me'),
                         anaerobe_lipids = c('19:0 cyclo'),
                         protozoa = c('20:4 w6,9,12,15'),
                         # May need to redefine the total_biomass. Right now, it's just all of the
                         # indicator lipids. Not sure if there should be others as well
                         total_biomass = c('8:0', '10:0', '10:0 2OH', '11:0',
                                           '12:0', '12:0 2OH', '12:0 3OH',
                                           '13:0', '14:0', '14:0 2OH',
                                           '14:0 3OH', '14:1', '15:0',
                                           '15:0 anteiso', '15:0 iso',
                                           '16:0 10me', '16:0 2OH', '16:0 iso',
                                           '16:1 w5c', '16:1 w7c', '16:1 w9c',
                                           '17:0', '17:0 iso', '17:0 anteiso',
                                           '17:0 cyclo', '17:1', '17:1 iso',
                                           '18:0', '18:0 10me', '18:1 w9c',
                                           '18:1 w9t', '18:2 w6,9c', '19:0',
                                           '19:0 cyclo', '19:1', '20:0')
  )
  # Naming list facilitates naming the output of calcs
  names(indicator_list) <- c('f_lipids', 'b_lipids', 'gram_pos_lipids',
                             'gram_neg_lipids', 'actino_lipids',
                             'anaerobe_lipids', 'protozoa', 'total_biomass')

  wide_df <- reshape(data = df[c('BatchDataFileName', 'DataFileName', 'Batch',
                                 'Name', 'd13C_corrected')],
                     timevar = 'Name',
                     idvar = c('BatchDataFileName', 'DataFileName', 'Batch'),  # w/o 2+ idvars, some vals are equal and get dropped
                     direction = 'wide')
  names(wide_df)<- gsub('d13C_corrected.', '', names(wide_df))
  #wide_df[is.na(wide_df)] <- 0


  indic_df <- cbind(wide_df['BatchDataFileName'],
                    wide_df['DataFileName'],
                    lapply(indicator_list,
                           function(x){
                             tmp_df <- wide_df[,
                                               names(wide_df) %in% x,
                                               drop = FALSE]
                             rowMeans(tmp_df, na.rm = TRUE)  # Rows with all NA will return NaN
                           }
                    )
  )

  indic_df_long <- reshape(indic_df, varying = names(indic_df)[c(3:ncol(indic_df))],
                           v.names = 'avg_d13C_corrected',
                           timevar = 'Indicator', idvar = 'BatchDataFileName',
                           times = names(indic_df)[c(3:ncol(indic_df))],
                           direction = 'long')

  indic_df_long <- merge(indic_df_long, sample_metadata, by = 'DataFileName',
                         all.x = TRUE)

  indic_df_long[is.na(indic_df_long[['avg_d13C_corrected']]), 'avg_d13C_corrected'] <- NA  #replace NaN's resulting from rows containing all NA

  return(indic_df_long)
}


#' Process del13C values
#'
#' This function corrects the raw del13C output for methylation and normalizes
#' to the del13C standard USGS40 (L-glutamic acid). This function is an update
#' of the correct_iso_base() function that uses a dataframe matching methanol
#' lot to batch so that we can use different methanol lots for an experiment.
#'
#' @param \code{df} peak list dataframe or tibble. This should contain the
#' following columns at a minimum: Batch, DataFileName, RetTimeSecs,
#' MajorHeightnA, TotalPeakArea1, DisplayDelta1, Name. This can be the output
#' from import_batch() or import_batch_multi().
#'
#' @param \code{d13c_correction} difference between measured USGS40 value and
#' expected USGS40 value (-26.39 per mil). If IonVantage output has already
#' been corrected, then enter 0.
#'
#' @param \code{meth_13c_df} Dataframe identifying the 13C signature of the
#' methanol used for each batch. This dataframe should have 2 columns: Batch
#' and Meth13C.
#'
#' @return Returns a tibble with columns from input dataframe and a corrected
#' d13C column.
#'
#' @examples
#'
#' @export
#'
#'
correct_d13c <- function(df, d13c_correction, meth_13c_df,
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

# TODO: use external dataframe to define indicator groups rather than defining
# within function.
# TODO: Add option for weighted or straight average

#' Calculate d13C for indicator groups
#'
#' This function aggregates individual lipid d13C signatures by indicator
#' group. Current implementation uses a straight average. Future versions
#' may offer the option to select between average and weighted average. This
#' replaces \code{indic_iso_base()} as of 8/10/2021
#'
#' @param \code{df} dataframe containing corrected d13C values for individual
#' lipids. Normally, this would be the output of correct_iso_base().
#' This should contain the following columns at a minimum: Batch, DataFileName,
#' BatchDataFileName, d13C_corrected, Name.
#'
#' @param \code{sample metadata}
#'
#' @return Returns a tibble with corrected d13C values for each indicator
#' group.
#'
#' @examples
#'
#' @export
#'

# Does the methanol correction have to be corrected for too?
# Or is this number already include instrument correction?

calculate_indicator_d13c <- function(df, sample_metadata){

  indicator_list <- list(f_lipids = c('18:1 w9c', '18:2 w6,9c'),
                         # Lipids representing total bacteria for f to b ratio not shown in Cameron's data, so for now, I will use lipids for Gram +/-, actino, and anaerobe
                         b_lipids = c('15:0 iso', '15:0 anteiso',
                                      '16:1 w7c', '16:0 10me',
                                      '18:1 w9c', '18:1 w9t',
                                      '18:2 w6,9c', '18:0 10me', '19:0 cyclo'),
                         gram_pos_lipids = c('15:0 iso', '15:0 anteiso'),
                         gram_neg_lipids = c('16:1 w7c', '18:1 w9t'),
                         actino_lipids = c('16:0 10me', '18:0 10me'),
                         anaerobe_lipids = c('19:0 cyclo'),
                         protozoa = c('20:4 w6,9,12,15'),
                         # May need to redefine the total_biomass. Right now, it's just all of the
                         # indicator lipids. Not sure if there should be others as well
                         total_biomass = c('8:0', '10:0', '10:0 2OH', '11:0',
                                           '12:0', '12:0 2OH', '12:0 3OH',
                                           '13:0', '14:0', '14:0 2OH',
                                           '14:0 3OH', '14:1', '15:0',
                                           '15:0 anteiso', '15:0 iso',
                                           '16:0 10me', '16:0 2OH', '16:0 iso',
                                           '16:1 w5c', '16:1 w7c', '16:1 w9c',
                                           '17:0', '17:0 iso', '17:0 anteiso',
                                           '17:0 cyclo', '17:1', '17:1 iso',
                                           '18:0', '18:0 10me', '18:1 w9c',
                                           '18:1 w9t', '18:2 w6,9c', '19:0',
                                           '19:0 cyclo', '19:1', '20:0')
  )
  # Naming list facilitates naming the output of calcs
  names(indicator_list) <- c('f_lipids', 'b_lipids', 'gram_pos_lipids',
                             'gram_neg_lipids', 'actino_lipids',
                             'anaerobe_lipids', 'protozoa', 'total_biomass')

  wide_df <- reshape(data = df[c('BatchDataFileName', 'DataFileName', 'Batch',
                                 'Name', 'd13C_corrected')],
                     timevar = 'Name',
                     idvar = c('BatchDataFileName', 'DataFileName', 'Batch'),  # w/o 2+ idvars, some vals are equal and get dropped
                     direction = 'wide')
  names(wide_df)<- gsub('d13C_corrected.', '', names(wide_df))
  #wide_df[is.na(wide_df)] <- 0


  indic_df <- cbind(wide_df['BatchDataFileName'],
                    wide_df['DataFileName'],
                    lapply(indicator_list,
                           function(x){
                             tmp_df <- wide_df[,
                                               names(wide_df) %in% x,
                                               drop = FALSE]
                             rowMeans(tmp_df, na.rm = TRUE)  # Rows with all NA will return NaN
                           }
                    )
  )

  indic_df_long <- reshape(indic_df, varying = names(indic_df)[c(3:ncol(indic_df))],
                           v.names = 'avg_d13C_corrected',
                           timevar = 'Indicator', idvar = 'BatchDataFileName',
                           times = names(indic_df)[c(3:ncol(indic_df))],
                           direction = 'long')

  indic_df_long <- merge(indic_df_long, sample_metadata, by = 'DataFileName',
                         all.x = TRUE)

  indic_df_long[is.na(indic_df_long[['avg_d13C_corrected']]), 'avg_d13C_corrected'] <- NA  #replace NaN's resulting from rows containing all NA

  return(indic_df_long)
}


#' Calculate d13C for indicator groups using weighted averages
#'
#' This function aggregates individual lipid d13C signatures by indicator
#' group. This function uses a weighted average, unlike
#' \code{calculate_indicator_d13c()}, which uses a straight average. The
#' weighting accounts for the lipid concentration and the carbon content of the
#' lipid. Future functions may allow a choice between straight and weighted
#' average.
#'
#' @param \code{d13c_df} dataframe containing corrected d13C values for
#' individual lipids. Normally, this would be the output of
#' \code{correct_iso_base()}. This should contain the following columns at a
#' minimum: Batch, DataFileName, BatchDataFileName, d13C_corrected, Name.
#'
#' @param \code{wt_df} dataframe containing individual lipid concentrations.
#' Normally, this would be the output of \code{process_peak_area()}.
#'
#'
#' @param \code{sample metadata}
#'
#' @return Returns a tibble with corrected d13C values for each indicator
#' group.
#'
#' @examples
#'
#' @export
#'
calc_wt_13C_indicators <- function(d13c_df, wt_df, sample_metadata){
  # TODO: would be nice to have some more checks built in to ensure that our
  # matrices match
  indicator_list <- list(f_lipids = c('18:1 w9c', '18:2 w6,9c'),
                         # Lipids representing total bacteria for f to b ratio not shown in Cameron's data, so for now, I will use lipids for Gram +/-, actino, and anaerobe
                         b_lipids = c('15:0 iso', '15:0 anteiso',
                                      '16:1 w7c', '16:0 10me',
                                      '18:1 w9c', '18:1 w9t',
                                      '18:2 w6,9c', '18:0 10me', '19:0 cyclo'),
                         gram_pos_lipids = c('15:0 iso', '15:0 anteiso'),
                         gram_neg_lipids = c('16:1 w7c', '18:1 w9t'),
                         actino_lipids = c('16:0 10me', '18:0 10me'),
                         anaerobe_lipids = c('19:0 cyclo'),
                         protozoa = c('20:4 w6,9,12,15'),
                         # May need to redefine the total_biomass. Right now, it's just all of the
                         # indicator lipids. Not sure if there should be others as well
                         total_biomass = c('8:0', '10:0', '10:0 2OH', '11:0',
                                           '12:0', '12:0 2OH', '12:0 3OH',
                                           '13:0', '14:0', '14:0 2OH',
                                           '14:0 3OH', '14:1', '15:0',
                                           '15:0 anteiso', '15:0 iso',
                                           '16:0 10me', '16:0 2OH', '16:0 iso',
                                           '16:1 w5c', '16:1 w7c', '16:1 w9c',
                                           '17:0', '17:0 iso', '17:0 anteiso',
                                           '17:0 cyclo', '17:1', '17:1 iso',
                                           '18:0', '18:0 10me', '18:1 w9c',
                                           '18:1 w9t', '18:2 w6,9c', '19:0',
                                           '19:0 cyclo', '19:1', '20:0')
  )

  # Because we use vector math to calc weighted means, df need to have same dimensions
  # Use intersect() to ensure that both df have the same lipids included
  # Lipids may differ between the two because we often filter IRMS peak height
  overlap <- intersect(wt_df[['Name']], d13c_df[['Name']])

  wt_df <- wt_df[which(wt_df[['Name']] %in% overlap), ]
  d13c_df <- d13c_df[which(d13c_df[['Name']] %in% overlap), ]



  # create wide version of lipid conc df (used as weights for average)
  df_wide_wt <- reshape(data = wt_df[c('BatchDataFileName', 'DataFileName', 'Batch', 'Name', 'nmol_g')],
                        timevar = 'Name', idvar = c('BatchDataFileName', 'DataFileName', 'Batch'),
                        direction = 'wide')

  names(df_wide_wt) <- gsub('nmol_g.', '', names(df_wide_wt))

  # Sort on BatchDataFileName to be consistent with wide_13c_df so that cells
  # in each df correspond to each other
  df_wide_wt <- df_wide_wt[order(df_wide_wt[['BatchDataFileName']]), ]


  cplfa_wide_df <- reshape(data = wt_df[c('BatchDataFileName', 'DataFileName', 'Batch', 'Name', 'C_PLFA')],
                           timevar = 'Name', idvar = c('BatchDataFileName', 'DataFileName', 'Batch'),
                           direction = 'wide')

  names(cplfa_wide_df)<- gsub('C_PLFA.', '', names(cplfa_wide_df))

  cplfa_wide_df <- cplfa_wide_df[order(cplfa_wide_df[['BatchDataFileName']]), ]


  # create wide version of d13C df
  wide_13c_df <- reshape(data = d13c_df[c('BatchDataFileName', 'DataFileName', 'Batch',
                                          'Name', 'd13C_corrected')],
                         timevar = 'Name',
                         idvar = c('BatchDataFileName', 'DataFileName', 'Batch'),  # w/o 2+ idvars, some vals are equal and get dropped
                         direction = 'wide')


  names(wide_13c_df)<- gsub('d13C_corrected.', '', names(wide_13c_df))



  # Merge the wide d13C df with the DataFileName col of the wide conc df
  # This is to ensure the number of rows and their order matches
  # Samples may have been dropped from the IRMS data if no peaks above cutoff
  wide_13c_df <- merge(df_wide_wt[c('BatchDataFileName')], wide_13c_df,
                       by = c('BatchDataFileName') ,all.x = T)

  # Order columns based on order of cols in other df to ensure cells correspond
  wide_13c_df <- wide_13c_df[names(df_wide_wt)]

  # Sort on BatchDataFileName to be consistent with df_wide_wt so that cells
  # in each df correspond to each other
  wide_13c_df <- wide_13c_df[order(wide_13c_df[['BatchDataFileName']]), ]

  # For each indicator group, run the weighted average calc
  wtd_13c_df <- cbind(df_wide_wt['BatchDataFileName'],
                      df_wide_wt['DataFileName'],
                      lapply(indicator_list,
                             function(x){
                               # Convert both df to matrix for easier math and
                               # select only lipids in the desired group
                               tmp_13c_df <- as.matrix(wide_13c_df[,
                                                                   names(wide_13c_df) %in% x,
                                                                   drop = FALSE])

                               # Data don't completely match between the d13C and abundance df because we
                               # filter peaks <1nA from isotope data, but not abundance data. This can
                               # cause issues, as the weighting factor will try to account for lipids that
                               # don't actully have isotope data to weight. The result is numbers smaller in
                               # magnitude than they should be. The following code chunk are to account for this
                               # and unsure that weights are only included for actually d13C peaks.
                               # We create a 'TRUE/FALSE matrix from the 13c data. We can then use this as a
                               #'mask' to eliminate cells from the abundance df that don't have corresponding
                               #'isotope data. We achieve this by simply multiplying the T/F matrix by the
                               #'abundance matrix. B/c T = 1 and F = 0, this preserves data points
                               #'corresponding to the isotope data and eliminating the rest.

                               tmp_mask <- !is.na(tmp_13c_df)

                               tmp_conc_df <- as.matrix(df_wide_wt[,
                                                                   names(df_wide_wt) %in% x,
                                                                   drop = FALSE]) * tmp_mask

                               tmp_cplfa_df <- as.matrix(cplfa_wide_df[,
                                                                       names(cplfa_wide_df) %in% x,
                                                                       drop = FALSE])

                               tmp_wt_df <- tmp_conc_df * tmp_cplfa_df
                               # For each sample (row), multiply each lipid d13C
                               # value by the corresponding lipid conc, then
                               # divide by the total conc for that inidcator group
                               rowSums(tmp_wt_df * tmp_13c_df, na.rm = TRUE) / rowSums(tmp_wt_df, na.rm = TRUE)
                             }
                      )
  )


  indic_df_long <- reshape(wtd_13c_df, varying = names(wtd_13c_df)[c(3:ncol(wtd_13c_df))],
                           v.names = 'avg_d13C_corrected',
                           timevar = 'Indicator', idvar = 'BatchDataFileName',
                           times = names(wtd_13c_df)[c(3:ncol(wtd_13c_df))],
                           direction = 'long')

  indic_df_long <- merge(indic_df_long, sample_metadata, by = 'DataFileName',
                         all.x = TRUE)

  indic_df_long[is.na(indic_df_long[['avg_d13C_corrected']]), 'avg_d13C_corrected'] <- NA  #replace NaN's resulting from rows containing all NA

  return(indic_df_long)
  #wtd_13c_df


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
