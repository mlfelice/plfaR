#' Import single batch PLFA peak list from Excel
#'
#' This is the first function in the plfaR workflow, loading named PLFA peak
#' list from a single batch for use in downstream analysis.
#'
#' @param file_path Full file path, including file extension for the Excel file
#' containing PLFA data. The input Excel file must contain Batch, DataFileName,
#' RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1, and Name columns.
#'
#' @return Tibble suitable for downstream analyses in plfaR. Contains the
#' input columns as well as BatchDataFileName.
#'
#' @examples
#'

# Maybe instead, we just have a display function that converts to wide
import_batch_tidy <- function(file_path){
  # This function imports named peak list for downstream processing/analysis
  #
  # Args:
  #   file_path: This is the full path+filename of source data
  #     - Input file should be in Excel format (xlsx or xls ok)
  #     - If more than one worksheet, relevant data must be
  #       included in tabs with following names:
  #       - 'named_peaks'
  #     - File should have the following columns with headers:
  #       - Batch
  #       - DataFileName
  #       - RetTimeSecs
  #       - MajorHeightnA
  #       - TotalPeakArea1
  #       - DisplayDelta1
  #       - Name
  #
  # To do: include data checks:
  #   - file format is correct
  #   - if more than one tab, make sure tab name is present
  #   - Make check peak names against the reference peak list
  #     - ensures that no mistakes in typing peaks (ie. extra spaces, etc)
  #   - We could alternatively make a function for using regex to replace
  #     the common errors (eg. trailing space, gamma instead of w, spaces
  #     after comma)
  readxl::read_excel(file_path, sheet = 'named_peaks', na = 'NA') %>%
    select(-BiomarkerRTBased, -Notes) %>%
    filter(!is.na(Name)) %>% # Remove unnamed peaks
             # add if statement so that if there is no batch column, but batch name is in
             #filename, it can create a batch column
             mutate(DisplayDelta1 = as.numeric(DisplayDelta1),
             #Batch = str_extract(string = file_path, pattern = '[Bb]atch ?[0-9]+'),
              BatchDataFileName = paste(Batch, DataFileName, sep = '_')) %>% # ? looks for 0 or 1)
             select(Batch, DataFileName, BatchDataFileName, everything()) # Ensure order
}

# one issue is that the quality control functions were designed to loop
# over a list of dataframes, not a single dataframe with all samples
# Need to modify so that both work on single dataframe (or list of df)
# Also want to modify this so that you can either supply a directory or
# a list of file paths
import_batch_multi_tidy <- function(dir, keyword){
  batch_files <- list.files(path = source_dir, pattern = keyword,
                            full.names = TRUE) %>%
    lapply(import_batch_tidy)

  all_batch_df <- bind_rows(batch_files)
  invisible(all_batch_df)
}


# Add something for importing metadata here
# Make it a requirement to have SampleID match between metadata and peak list

