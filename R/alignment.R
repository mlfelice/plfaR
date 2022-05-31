# TODO: replace tidyverse with base R
# TODO: include sample datasets that would allow running of examples
# TODO: add export functions for saving processed data
#### Use openxlsx to create multi-sheet files that document the processing choices made along with the data

#' Ensure proper formatting of sample names
#'
#' Eliminate special characters and ensure proper formatting of sample names to
#' comply with GCalignR formatting requirements. This function is used
#' internally in the \code{df_to_gcinput()} function
#'
#' Input must be a named list of dataframes, where each list element is named
#' with corresponding sample name, with no special characters in sample names.
#' Each dataframe must have at minimum one column for  retention time, but can
#' also include columns for additional data. Retention time column must be
#' named 'rt'
#'
#' @param \code{vector} Character vector containing sample names.
#'
#' @return Returns a vector of sample names with special characters removed.
#' ' ' and '-' are replaced with '_' and '_to_', respectively, and removes file
#'  '.raw' file extension. The special characters replaced are currently
#'  limited, but future versions will be more comprehensive.
#'
#' @examples
#'
#'
#'
format_names_gcalignr <- function(vector){
  # Function replaces ' ' and '-' with '_' and '_to_', respectively, and
  # removes file '.raw' file extension. Utlimately, should replace all
  # special characters, but these are the only ones found in SPRUCE data.
  # This is used in df_to_gc_input()

  # strips special characters and file extension
  samp_name <- vector %>%
    str_replace_all(pattern = '\\s', replacement = '_') %>%
    str_replace_all(pattern = '-', replacement = '_to_') %>%
    str_remove('.raw')

  samp_name
}
# Example
#format_names_gcalignr(batch_fame_df[['BatchDataFileName']])



#' Format GC data for GCalignR
#'
#' Convert GC peak lists into the format required for alignment using the
#' GCalignR package. This function formats data to meet the folowing
#' requirements:
#'
#' Input must be a named list of dataframes, where each list element is named
#' with corresponding sample name, with no special characters in sample names.
#' Each dataframe must have at minimum one column for  retention time, but can
#' also include columns for additional data. Retention time column must be
#' named 'rt'
#'
#' @param \code{batch_df} peak list dataframe or tibble from GC run in
#' IonVantage. The recommended column are as follows:
#' Batch, DataFileName, RetTimeSecs, MajorHeightnA, TotalPeakArea1,
#' DisplayDelta1, Name. This can be the output from \code{import_batch()} or
#' \code{import_batch_multi()}.
#'
#' @return Returns a nanmed list of dataframes, with each dataframe
#' reepresenting the peak list for a given sample.
#'
#' @examples
#'
#' @export
#'
df_to_gcinput <- function(batch_df){ # same as above, but avoids for loop
  #
  # Converts dataframe of imported GC batch to suitable input for GCalignR
  #
  # Args:
  #   batch_df: dataframe created from Isoprime batch log file (cols = ID_With batch,
  #       Name, RT (Sec),	Height (nA),	Corrected 13C	Peak Name	Sum Peak Area)
  #     - This should be direct import of 'F1CO2-log' sheet of results
  #
  # Returns:
  #   input_list: list of dataframes formatted for direct input to align_chromatograms()
  #     - each element will be sample name (whitespace and special char replaced with "_")
  #     - dataframe cols = "rt" and "area"

  # Convert input dataframe to list of dataframes w/proper colnames
  input_list <- lapply(unique(batch_df$DataFileName),
                       function(x){df <- batch_df %>%
                         filter(DataFileName == x) %>%
                         select(RetTimeSecs, TotalPeakArea1, MajorHeightnA) %>%
                         arrange(RetTimeSecs) %>% # Need to sort on RT or interferes w/alignment
                         rename(rt = RetTimeSecs, area = TotalPeakArea1,
                                height = MajorHeightnA)
                       })


  # Extract sample names from dataframe and replace special char w/_
  sample_names <- vapply(unique(batch_df$DataFileName),
                         function(x){
                           format_names_gcalignr(x)
                         },
                         FUN.VALUE = 'character')

  # Replace column names
  names(input_list) <- unique(sample_names)

  input_list
}
# Example
#check_input(df_to_gcinput(batch_fame_df))


#' Apply lipid names to aligned chromatograms
#'
#' Applies lipid names to all samples Using the specified named reference
#' chromatogram
#'

#'
#' @param \code{gc_align} GCalignR object containing peak alignment.
#'
#' @param \code{ref_samp_df} Named peak list dataframe in same format as output
#' \code{from import_batch()} or \code{import_batch_multi()}.
#'
#' @param \code{ref_samp} Character string of the name of the sample used as
#' the reference chromatogram for naming. THe sample name must be in proper
#' GCalignR format.
#'
#' @return This version returnss a long-form dataframe while maintaining peak
#' height and area data. Output includes following columns: BatchDataFileName,
#' RetTimeSecs, TotalPeakArea1, MajorHeightnA, Peak_Name
#'
#' @examples
#'
#' @export
#'
name_alignment <- function(gcalign, ref_samp_df, ref_samp){
  # This function attaches names to peak list from GCalignR. This version
  # outputs long-form dataframe while maintaining peak height and area
  # data. Output includes following columns: DataFileName,
  # RetTimeSecs, TotalPeakArea1, MajorHeightnA, Peak_Name

  # extract aligned peak list from GCalignR obj
  # sort = F should prevent resorting df, but still sorts on the merge col
  # Identical samp names will get prefixed w/name of list item (eg. rt., area.)
  gcalign[['aligned']][['rt']] <- merge(gcalign[['aligned']][['rt']],
                                        ref_samp_df,
                                        # could make this more automatic if desired
                                        # Is reference stored anywhere?
                                        by.x = c(ref_samp), # Formatted name of sample used as reference
                                        by.y = c('RetTimeSecs'),
                                        all.x = T,
                                        sort = F) %>%
    # Need to re-sort on mean_RT, which is how all alignments are initially sorted
    arrange(mean_RT)
  # also sort the area and height df in the aligned list to ensure matching
  # Before I had both sorted, the sometimmes failed for unknown reason
  gcalign[['aligned']][['area']] <- arrange(gcalign[['aligned']][['area']], mean_RT)
  gcalign[['aligned']][['height']] <- arrange(gcalign[['aligned']][['height']], mean_RT)

  # Convert to long format
  # names_sep arg specifies new separate val cols created for each prefix
  # portion after prefix will become entry in DataFileName col
  tmp_df <- do.call('cbind', gcalign[['aligned']]) %>%

    pivot_longer(cols = -rt.Name,
                 names_to = c('.value', 'DataFileName'),
                 names_sep = '\\.') %>%
    rename(Name = rt.Name, RetTimeSecs = rt,
           TotalPeakArea1 = area, MajorHeightnA = height)

  # If aligned rt, area, height somehow get out of order, data will get
  # associated with wrong peak. This adds a check for correct binding
  if(nrow(tmp_df %>% filter(RetTimeSecs == 0 & TotalPeakArea1 != 0,
                            RetTimeSecs == 0 & MajorHeightnA != 0)) > 0){

    warning('Peak data improperly aligned with retention time')

  }

  # Filter RetTimeSecs == 0, since these represent no peak
  out_df <- tmp_df %>% filter(RetTimeSecs != 0) %>%
    select(DataFileName, Name, RetTimeSecs,
           MajorHeightnA, TotalPeakArea1)

  out_df

}
# Example
#print(name_alignment(al, fame_ref, format_names_gcalignr(ref_samp)))

#' Run peak alignment and assign peak names
#'
#' Uses a single reference chromatogram to align a batch of chromatograms and
#' assign peak names based on the reference. This function wraps all stages of
#' data prep.
#'
#'
#' @param \code{df} peak list dataframe or tibble from GC run in
#' IonVantage. The recommended column are as follows:
#' Batch, DataFileName, RetTimeSecs, MajorHeightnA, TotalPeakArea1,
#' DisplayDelta1, Name. This can be the output from \code{import_batch()} or
#' \code{import_batch_multi()}.
#'
#' @param \code{reference} Character string of the name of file used as your
#' reference chromatogram for alignment and naming. This should be in the same
#' format as in the peak list dataframe (i.e. not formatted for GCalignR)
#'
#' @return
#'
#' @examples
#'
#' @export
#'
align_name_chrom <- function(df, #ref_samp,
                             rt_col_name = 'rt', # retention time variable name
                             rt_cutoff_low = 400, # remove peaks below 15 Minutes
                             rt_cutoff_high = 3600, # remove peaks exceeding 45 Minutes
                             reference = NULL, # choose automatically NOTE: Needs to be in format matching source df DataFileName
                             max_linear_shift = 25, # max. shift for linear corrections
                             max_diff_peak2mean = 20, # max. distance of a peak to the mean across samples
                             min_diff_peak2peak = 20, # min. expected distance between peaks
                             blanks = NULL, # negative control
                             delete_single_peak = FALSE, # delete peaks that are present in just one sample
                             write_output = NULL,# add variable names to write aligned data to text files)
                             remove_empty = F){ # removes row if no peak

  # This is just a giant wrapper function created specifically for the purpose
  # of quickly testing automated naming of PLFA batches using GCalignR. Input a
  # PLFA peak list, the filename of a reference sample, and any
  # GCalignR::align_chromatograms() args. This will return a dataframe showing
  # retention time of each peak for each file in the batch along with the
  # automatically generated name and the originally assigned name. Must use
  # named peak lists for this to work properly.
  # Additionally, this will print GCalignR::gc_heatmap() and summary tibbles
  # showing the number of orginal to auto-assigned name mismatches for each
  # and overa the entire batch.
  # (user-defined fun) NOTE: spaces (and special chars) converted to underscore
  input_list <- df_to_gcinput(df)

  # Check the input for formatting errors
  GCalignR::check_input(input_list)

  # summarize distribution of dist between peaks - use to set alignment params
  GCalignR::peak_interspace(data = input_list, rt_col_name = 'rt')

  ### NEW ###
  # Extract named peak list for reference sample
  # This will be used to assign names to aligned peak lists
  # use a dataframe to preserve original names so we can restore in the end
  sample_name_df <- df %>%
    mutate(GCalignName = format_names_gcalignr(DataFileName)) %>%
    select(OriginDataFileName = DataFileName, GCalignName) %>%
    distinct()

  ### NEW ###

  # store a dataframe of the original names and corresponding GCalignR names in order
  # to restore at the end
  #origin_name_df <- data.frame()

  name_ref <- df %>%
    filter(DataFileName == reference) %>%
    select(RetTimeSecs, Name)  #####switched from Peak_Name to Name to be consistent with plfaR format

  # Remove special characters and make sure supplied sample name matches
  # GCalignR format
  ref_samp <- format_names_gcalignr(reference)

  # Function just wraps a few steps that I'm using to test how alignment works
  # Mostly just modifying the input list and re-running


  align <- GCalignR::align_chromatograms(data = input_list, # input data
                                         rt_col_name = rt_col_name, # retention time variable name
                                         rt_cutoff_low = rt_cutoff_low, # remove peaks below 15 Minutes
                                         rt_cutoff_high = rt_cutoff_high, # remove peaks exceeding 45 Minutes
                                         reference = ref_samp, # choose automatically
                                         max_linear_shift = max_linear_shift, # max. shift for linear corrections
                                         max_diff_peak2mean = max_diff_peak2mean, # max. distance of a peak to the mean across samples
                                         min_diff_peak2peak = min_diff_peak2peak, # min. expected distance between peaks
                                         blanks = blanks, # negative control
                                         delete_single_peak = delete_single_peak, # delete peaks that are present in just one sample
                                         write_output = write_output,
                                         remove_empty = remove_empty) # add variable names to write aligned data to text files

  # Display heatmap figure of alignment
  #print(GCalignR::gc_heatmap(align))

  named <- name_alignment(align, name_ref, ref_samp)

  final_named_df <- named %>%
    left_join(sample_name_df, by = c('DataFileName' = 'GCalignName')) %>%
    select(-DataFileName, DataFileName = OriginDataFileName)


}
