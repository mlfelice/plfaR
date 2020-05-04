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
# Need to make decisions about what should go in ref lipids df, and change code accordingly
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
    message(
            paste0('Warning: Standards missing from at least one sample.\n',
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
relative_peak_area <- function(list){
  mean_list <- lapply(list,
                      function(x){
                        aggregate(TotalPeakArea1 ~ Batch, data = x, FUN = mean, drop = TRUE)
                      })

  mean_df <- do.call('rbind', mean_list)

  mean_df['fractional_mean_area'] <- sapply(mean_df[2], function(x){x/max(x)})

  sum_list <- lapply(list,
                     function(x){
                       aggregate(TotalPeakArea1 ~ Batch, data = x, FUN = sum, drop = TRUE)
                     })

  sum_df <- do.call('rbind', sum_list)

  sum_df['fractional_sum_area'] <- sapply(sum_df[2], function(x){x/max(x)})

  summary_df <- merge(mean_df, sum_df, by = c('Batch'))

  colnames(summary_df) <- c('Batch', 'MeanPeakArea', 'RelMeanPeakArea', 'SumPeakArea', 'RelSumPeakArea')
  return(summary_df)
}

# function for summarizing standards missing expected peaks, average peak area for standards (or indiv classes of standards)

# Some stats on peak heights

# could also add an option to filter peak heights below a certain threshold

# Function to determine intra-batch variation in RT of major/reliable lipids
rt_var <- function(batch_df){
  rt_sd_vec <- tapply(batch_df[['RetTimeSecs']],
                      INDEX = batch_df[['Name']],
                      FUN = sd, na.rm = TRUE)
  data.frame(Batch = unique(batch_df[['Batch']]), Name = names(rt_sd_vec), StDev = rt_sd_vec)
}

rt_var_list <- function(list){
  summary_list <- lapply(list, rt_var)
  summary_df <- do.call('rbind', summary_list)
  summary_df <- reshape(data = summary_df, idvar = 'Name', timevar = 'Batch', direction = 'wide')
  colnames(summary_df) <- gsub(pattern = 'StDev.', replacement = 'Batch', x = colnames(summary_df))
  rownames(summary_df) <- summary_df[['Name']]
  return(summary_df)
}


# Function for counting fraction of peaks/samples that make it through quality control steps

#Function for checking for lipids not matching the database (this will avoid issues with typos)
