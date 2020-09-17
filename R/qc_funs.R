# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
# Add package dependency:
# devtools::use_package('package name')
#
# Include data in the package
# devtools::use_data('data file name') # set internal = TRUE to keep private
#

hello <- function() {
  print("Hello, world!")
}

###############################################################################
# This file contains tidyverse-based functions for assessing the quality of
# PLFA peak naming. This is not the corrent implementation, and may need
# debugging to work with more current versions of dplyr and other tidyverse
# functions.
###############################################################################

find_duplicates_tidy <- function(df){

  duplicates_df <- df %>%
    group_by(DataFileName, Name) %>%
    summarise(Count = n()) %>%
    filter(Count > 1)

  if(any(duplicates_df[['Count']] > 1)){
    cat('Warning: Duplicate peaks detected. Check naming before proceeding\n\n')

  }

  return(duplicates_df)
}

find_missing_tidy <- function(df, lipids = c('13:0', '16:0', '19:0')){
  # input df must be grouped by DataFileName for this to work properly
  # This should now work without going through the intermediate check_quality() step

  biomarkers <- unique(df[['Name']])

  missing_std_df <- df %>%
    group_by(DataFileName) %>%
    tidyr::complete(Name = !!biomarkers) %>% # Makes missing vals explicit
    filter(is.na(TotalPeakArea1) & Name %in% !!lipids) %>%
    select(Batch, DataFileName, Name)

  if(nrow(missing_std_df) > 0){
    cat('Warning: Standards missing from at least one sample.\n')
    cat('Check data file and/or chromatograms before proceeding\n\n')
    #cat(missing_std_df[which(missing_std_df[['Count']] > 1),][['DataFileName']])
  }

  return(missing_std_df)

}

#' Calculate distribution of lipids
#'
#' \code{count_lipids_tidy} tabulates the number of occurrences and frequency
#' of detection of each lipid and returns as a tibble.
#'
#' @param df Dataframe or tibble containing the following columns: Batch,
#' DataFileName, RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1,
#' Name.
#'
#' @return A tibble in long format containing the original columns as well as
#' Count and LipidFrequency.
#'
#' @examples
#'
count_lipids_tidy <- function(df){

  n_samples <- length(unique(df[['DataFileName']]))

  lipid_freq_df <- df %>%
    mutate(Count = 1) %>%
    group_by(Name) %>%
    summarise(LipidFrequency = sum(Count)/!!n_samples)

  return(lipid_freq_df)
}

#' Inspect PLFA peak list quality.
#'
#' \code{quality_check_tidy} returns a list of dataframes containing the following
#' information related to the quality of the data in the peak list:
#' Samples with duplicate peak names
#' Samples that are missing one or more of the standards peaks 13:0, 16:0, 19:0
#' The distribution of lipids among samples.
#' Function will also notifiy you of batches with these issues in the console
#' whether or not you store the return values
#'
#' @param df Dataframe or tibble containing the following columns: Batch,
#' DataFileName, RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1,
#' Name.
#'
#' @return List of dataframes indicating samples with duplicate peak names,
#' samples missing standard peaks, and the distribution of lipids among
#' samples.
#'
#' @examples
#' If you receive a warning, you can easily pull out relevant details using
#' lapply().
#'
#' \code{lapply(qc_stats_2016, function(x){ # pull out dup lipids to ID samples
#' x$duplicate_lipids
#' }
#' )
#' }
#'
#'
quality_check_tidy <- function(df){


  qa_funs <- list(duplicate_lipids = find_duplicates_tidy, missing_lipids = find_missing_tidy,
                  lipid_frequency = count_lipids_tidy)

  qa_stats_list <- lapply(split(df, f = df[['Batch']]), function(x){
                    cat(paste0('Batch: ', unique(x[['Batch']]), '\n------\n'))
                    lapply(qa_funs, function(f){f(x)})
                   })

  #qa_stats_list <- lapply(qa_funs, function(f){f(df)})
  #names(qa_stats_list) <- c('duplicate_peak_names', 'missing_standard_peaks',
  #                          'lipid_distribution')

  names(qa_stats_list) <- unique(as.character(df[['Batch']]))

  invisible(qa_stats_list) # Create list of for assignemnt w/o printing

}


## Include a function here that computes fractional difference in peak area
## between batches as something to determine whether or not normalization is
## appropriate.
