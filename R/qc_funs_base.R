###############################################################################
# This file includes the current implementation of functions used for assessing
# the quality of PLFA peak naming. Functions identify duplicate named peaks,
# identify peaks that are expected but not found, and provide information on
# frequency with which peaks appear in samples.
###############################################################################


#' Identify duplicate peak names within samples
#'
#' \code{find_duplicates} Identify samples with naming errors resulting in a
#' lipid name being given to multiple peaks in a sample.
#'
#' @param df Dataframe or tibble containing the following columns: Batch,
#' DataFileName, RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1,
#' Name.
#'
#' @return Dataframe indicating the sample containing duplicate lipids, and the
#' lipids within the sample that were duplicated.
#'
#' @examples
#'
#' @export
#'
find_duplicates <- function(df){

  tmp_df <- aggregate(TotalPeakArea1~BatchDataFileName+Name, data = df,
                      FUN = length, drop = TRUE)
  names(tmp_df)[names(tmp_df) == 'TotalPeakArea1'] <- 'Count'
  duplicates_df <- tmp_df[tmp_df[['Count']] > 1, ]

  if(any(duplicates_df[['Count']] > 1)){
    cat('Warning: Duplicate peaks detected. Check naming before proceeding\n\n')

  }

  return(duplicates_df)
}

#` Identify missing lipids that are expected in all samples
#'
#' \code{find_missing} Identifies samples that are missing lipids that are
#' expected to be in all samples. For example, 13:0 and 19:0 lipids are added
#' to all samples as standard and surrogate standard, so these lipids missing
#' may indicate issues with the run or with the naming. For peat samples, 16:0
#' is expected to appear in all samples, as this lipid is present in Sphagnum.
#'
#' @param df Dataframe or tibble containing the following columns: Batch,
#' DataFileName, RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1,
#' Name.
#'
#' @param lipids Character vector listing the lipids that should be found in
#' all samples.
#'
#' @return Dataframe indicating the sample containing missing the specified
#' lipids, and the lipids within the sample that were missing.
#'
#' @examples
#'
#' @export
#'
find_missing <- function(df, lipids = c('13:0', '16:0', '19:0')){

  data_files <- unique(df$DataFileName)  # list of all sample fnames in batch
  # '_' to avoid conflict btwn df
  lipids_df <- data.frame(FAME = lipids, Batch_ = df[['Batch']][1])
  # loop through each sample (filename) to create df of missing lipids
  tmp_list <- lapply(data_files,
                     function(x){
                       # for each samp left join named peaks to col of all microbe lipids
                       # if lipid absent from sample, all fields for FAME are NA
                       tmp <- merge(lipids_df, df[df[['DataFileName']] == x, ], by.x = 'FAME',
                                    by.y = 'Name', all.x = TRUE)

                       # fill NA DataFileName cols - uses first NA val (all are the same)
                       tmp[is.na(tmp[['DataFileName']]), 'DataFileName'] <- x

                       # select rows for lipids missing from sample and in list of std lipids
                       tmp[is.na(tmp[['TotalPeakArea1']]), ]

                     }
  )

  # combine list to df showing all sample/lipid combos missing
  missing_std_df <- do.call('rbind', tmp_list)

  # print warning if standards were missing from any samps (user can check)
  if(nrow(missing_std_df) > 0){
    message(paste0('Warning: Standards missing from at least one sample.\n',
                   'Check data file and/or chromatograms before proceeding\n\n')
    )
    #cat(missing_std_df[which(missing_std_df[['Count']] > 1),][['DataFileName']])
  }

  return(missing_std_df[c('Batch_', 'DataFileName', 'FAME')])

}

#' Calculate distribution of lipids
#'
#' \code{count_lipids} tabulates the number of occurrences and frequency
#' of detection of each lipid and returns as a tibble.
#'
#' @param df Dataframe or tibble containing the following columns: Batch,
#' DataFileName, RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1,
#' Name.
#'
#' @return A tibble in long format containing the name of each lipid and the
#' frequency with which it is detected in a sample.
#'
#' @examples
#'
#' @export
#'
count_lipids <- function(df){

  n_samples <- length(unique(df[['DataFileName']]))

  tmp_df <- aggregate(TotalPeakArea1~BatchDataFileName+Name, data = df,
            FUN = length, drop = TRUE)
  lipid_freq_df <- aggregate(TotalPeakArea1~Name, data = tmp_df,
            FUN = function(x){sum(x)/n_samples}, drop = TRUE)

  names(lipid_freq_df) <- c('Name', 'Freq')

  return(lipid_freq_df)
}

#' Inspect PLFA peak list quality.
#'
#' \code{quality_check} returns a list of dataframes containing the following
#' information related to the quality of the data in the peak list:
#' Samples with duplicate peak names
#' Samples that are missing one or more of the standards peaks 13:0, 16:0, 19:0
#' The distribution of lipids among samples.
#' Function will also notifiy you of batches with these issues in the console
#' whether or not you store the return values
#'
#' @param df Dataframe or tibble containing the following columns: Batch,
#' DataFileName, RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1,
#' Name. df can also be a list of dataframes, with each element following the
#' above format.
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
#' @export
#'
quality_check <- function(df){

  qa_funs <- list(duplicate_lipids = find_duplicates,
                  missing_lipids = find_missing,
                  lipid_frequency = count_lipids)

  dtype <- check_format(df)

  if (dtype == 'data.frame'){

    ls <- split(df, f = df[['Batch']])

    qa_stats_list <- lapply(ls, function(x){
      cat(paste0('Batch: ', unique(x[['Batch']]), '\n------\n'))
      lapply(qa_funs, function(f){f(x)})
    })

    batch_names <- vapply(ls, function(x){unique(as.character(x[['Batch']]))},
                          FUN.VALUE = 'character')

  } else if (dtype == 'list'){
    qa_stats_list <- lapply(df, function(x){
      cat(paste0('Batch: ', unique(x[['Batch']]), '\n------\n'))
      lapply(qa_funs, function(f){f(x)})
    })

    batch_names <- vapply(df, function(x){unique(as.character(x[['Batch']]))},
                          FUN.VALUE = 'character')
  }

  # TO DO: Make this message show which samples had no associated batch info
  #   This might require us to move this up into the lapply() loops above
  if(NA %in% batch_names){message(paste0('Warning: Some samples in this batch',
                                         'are missing batch data'))}

  names(qa_stats_list) <- paste0('Batch_',
                                 batch_names[which(!is.na(batch_names))])

  invisible(qa_stats_list) # Create list of for assignemnt w/o printing

}

# function for summarizing the fractional difference between batches

# function for summarizing standards missing expected peaks, average peak area for standards (or indiv classes of standards)

# Some stats on peak heights

# could also add an option to filter peak heights below a certain threshold

# Function to determine intra-batch variation in RT of major/reliable lipids
