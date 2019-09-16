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

  data_files <- unique(df$DataFileName)
  all_lipids <- lipid_reference['fame']

  tmp_list <- lapply(data_files,
    function(x){
      tmp <- merge(all_lipids, df[df[['DataFileName']] == x, ], by.x = 'fame',
            by.y = 'Name', all.x = TRUE)
      tmp[is.na(tmp[['DataFileName']]), 'DataFileName'] <-  #Would like to find a more elegant solution for filling missing values
        tmp[!is.na(tmp[['DataFileName']]), 'DataFileName'][[1]]
      tmp[is.na(tmp[['Batch']]), 'Batch'] <-
        tmp[!is.na(tmp[['Batch']]), 'Batch'][[1]]
      tmp[is.na(tmp[['TotalPeakArea1']]) & tmp[['fame']] %in% lipids, ]
    }
  )
  missing_std_df <- do.call('rbind', tmp_list)

  if(nrow(missing_std_df) > 0){
    cat('Warning: Standards missing from at least one sample.\n')
    cat('Check data file and/or chromatograms before proceeding\n\n')
    #cat(missing_std_df[which(missing_std_df[['Count']] > 1),][['DataFileName']])
  }

  return(missing_std_df[c('Batch', 'DataFileName', 'fame')])

}

count_lipids_base <- function(df){

  n_samples <- length(unique(df[['DataFileName']]))

  tmp_df <- aggregate(TotalPeakArea1~BatchDataFileName+Name, data = df,
            FUN = length, drop = TRUE)
  lipid_freq_df <- aggregate(TotalPeakArea1~Name, data = tmp_df,
            FUN = function(x){sum(x)/n_samples}, drop = TRUE)

  return(lipid_freq_df)
}

quality_check_base <- function(df){

  qa_funs <- list(duplicate_lipids = find_duplicates_base,
                  missing_lipids = find_missing_base,
                  lipid_frequency = count_lipids_base)

  dtype <- check_format(df)

  if (dtype == 'data.frame'){

    qa_stats_list <- lapply(split(df, f = df[['Batch']]), function(x){
      cat(paste0('Batch: ', unique(x[['Batch']]), '\n------\n'))
      lapply(qa_funs, function(f){f(x)})
    })
  } else if (dtype == 'list'){
      qa_stats_list <- lapply(df, function(x){
        cat(paste0('Batch: ', unique(x[['Batch']]), '\n------\n'))
        lapply(qa_funs, function(f){f(x)})
      })
  }

  names(qa_stats_list) <- unique(as.character(df[['Batch']]))

  invisible(qa_stats_list) # Create list of for assignemnt w/o printing

}

# function for summarizing the fractional difference between batches

# function for summarizing standards missing expected peaks, average peak area for standards (or indiv classes of standards)

# Some stats on peak heights

# could also add an option to filter peak heights below a certain threshold
