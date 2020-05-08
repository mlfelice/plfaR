find_duplicates_base <- function(df){

  tmp_df <- aggregate(TotalPeakArea1~BatchDataFileName+Name, data = df,
                      FUN = length, drop = TRUE)
  names(tmp_df)[names(tmp_df) == 'TotalPeakArea1'] <- 'Count'
  duplicates_df <- tmp_df[tmp_df[['Count']] > 1, ]

  if(any(duplicates_df[['Count']] > 1)){
    cat('Warning: Duplicate peaks detected. Check naming before proceeding\n\n')

  }

  return(duplicates_df)
}

find_missing_base <- function(df, lipids = c('13:0', '16:0', '19:0')){

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

count_lipids_base <- function(df){

  n_samples <- length(unique(df[['DataFileName']]))

  tmp_df <- aggregate(TotalPeakArea1~BatchDataFileName+Name, data = df,
            FUN = length, drop = TRUE)
  lipid_freq_df <- aggregate(TotalPeakArea1~Name, data = tmp_df,
            FUN = function(x){sum(x)/n_samples}, drop = TRUE)

  names(lipid_freq_df) <- c('Name', 'Freq')

  return(lipid_freq_df)
}

quality_check_base <- function(df){

  qa_funs <- list(duplicate_lipids = find_duplicates_base,
                  missing_lipids = find_missing_base,
                  lipid_frequency = count_lipids_base)

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
