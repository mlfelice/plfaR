###############################################################################
# Functions for importing PLFA data and formatting for use in processing
# functions.
#
# This file contains the most current implemented version of the functions
###############################################################################

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
#' @export

# Maybe instead, we just have a display function that converts to wide
import_batch <- function(file_path){
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
  column_names <- c('Batch', 'DataFileName', 'RetTimeSecs', 'MajorHeightnA',
                    'TotalPeakArea1', 'DisplayDelta1', 'Name')
  column_types <- c('text', 'text', 'numeric', 'numeric', 'numeric', 'guess',
                    'text')

  tmp_df <- readxl::read_excel(file_path, sheet = 'named_peaks',
                               range = readxl::cell_cols('A:G'),
                               col_types = column_types, na = 'NA')
  tmp_df <- tmp_df[column_names] #This won't be necessary if we import the correct format
  tmp_df[['DisplayDelta1']] <- as.character(tmp_df[['DisplayDelta1']])  # Sometimes this col has "#TypeMismatch", which throws an error for readxl if you specify numeric

  tmp_df <- tmp_df[!is.na(tmp_df[['Name']]), ]  # Remove unnamed peaks
  tmp_df[['BatchDataFileName']] <- paste(tmp_df[['Batch']], tmp_df[['DataFileName']], sep = '_')

  # Convert area and d13C columns to numeric if not already
  tmp_df['TotalPeakArea1'] <- as.numeric(tmp_df[['TotalPeakArea1']])
  tmp_df['DisplayDelta1'] <- as.numeric(tmp_df[['DisplayDelta1']])
    # add if statement so that if there is no batch column, but batch name is in
    #filename, it can create a batch column
#    mutate(DisplayDelta1 = as.numeric(DisplayDelta1),
#           Batch = str_extract(string = file_path, pattern = '[Bb]atch ?[0-9]+'),
#           BatchDataFileName = paste(Batch, DataFileName, sep = '_')) %>% # ? looks for 0 or 1)
  imported_df <- tmp_df[c('BatchDataFileName', column_names)] # Ensure order
}

#' Import multiple PLFA peak lists (Excel files) in a directory
#'
#' This is the first function in the plfaR workflow, loading named PLFA peak
#' list from mutliple batches for use in downstream analysis.
#'
#' @param dir Full path to the directory containing the PLFA data Excel files
#' to be imported. The input Excel files must contain Batch, DataFileName,
#' RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1, and Name columns.
#'
#' @return A named list. Each element is a tibble suitable for downstream
#' analyses in plfaR. Tibbles contain the input columns as well as
#' BatchDataFileName.
#'
#' @examples
#'
#' @export

import_batch_multi <- function(dir, keyword){
  # This function imports all named peak lists in a directory for downstream
  # processing/analysis
  #
  # Args:
  #   dir: This is full path to the directory containing PLFA data
  #     - Input files should be in Excel format (xlsx or xls ok)
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
  batch_files <- list.files(path = dir, pattern = keyword,
                            full.names = TRUE)
  batch_list <- lapply(batch_files, import_batch_base)

  #all_batch_df <- bind_rows(batch_files)
  invisible(batch_list)
}

# would be cool to have a merge emthod like in phyloseq object
