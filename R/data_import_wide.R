#source_path <- file.choose()

##############
# Make function for importing wide gc data
import_batch_wide <- function(file_path){
  column_names <- c('Batch', 'DataFileName', 'RetTimeSecs', 'MajorHeightnA',
                    'TotalPeakArea1', 'DisplayDelta1', 'Name')
  column_types <- c('text', 'text', 'numeric', 'numeric', 'numeric', 'guess',
                    'text')

  tmp_df <- readxl::read_excel(file_path, sheet = 'named_peaks',
                               range = readxl::cell_cols('A:G'),
                               col_types = column_types, na = 'NA')
  tmp_df <- tmp_df[column_names] #This won't be necessary if we import the correct format
  tmp_df[['DisplayDelta1']] <- as.numeric(tmp_df[['DisplayDelta1']])  # Sometimes this col has "#TypeMismatch", which throws an error for readxl if you specify numeric

  tmp_df <- tmp_df[!is.na(tmp_df[['Name']]), ]  # Remove unnamed peaks
  tmp_df[['BatchDataFileName']] <- paste(tmp_df[['Batch']], tmp_df[['DataFileName']], sep = '_')

  Check_duplicate_lipids(tmp_df)

  import_funs <- list(
    import_area <- function(batch_df){
      reshape::cast(batch_df, DataFileName ~ Name,  # make wide
                    value = 'TotalPeakArea1', fun.aggregate = sum) #%>%  # vals don't change
    },

    import_height <- function(batch_df){
      reshape::cast(batch_df, DataFileName ~ Name,  # make wide
                    value = 'MajorHeightnA', fun.aggregate = sum) #%>%  # vals don't change
    },

    import_delta <- function(batch_df){
      reshape::cast(batch_df, DataFileName ~ Name,  # make wide
                    value = 'DisplayDelta1', fun.aggregate = sum) #%>%  # vals don't change
    }

  )

  batch_list <- lapply(import_funs, function(f){f(tmp_df)})
  names(batch_list) <- c('peak_area', 'peak_height', 'd13c')
  return(batch_list)
}

import_batch_multi_wide <- function(dir, keyword){
  batch_files <- list.files(path = source_dir, pattern = keyword,
                            full.names = TRUE)
  batch_list <- lapply(batch_files, import_batch_wide)

  #all_batch_df <- bind_rows(batch_files)
  invisible(batch_list)
}

#batch_df <- readxl::read_excel(source_path, sheet = 'named_peaks', na = 'NA') %>%
#  select(DataFileName, MajorHeightnA, TotalPeakArea1, DisplayDelta1, Name)

#batch_df <- batch_df[!is.na(batch_df[['Name']]),]

#Check_duplicate_lipids(batch_df)

# During reshape, NA's become 0, I think
#batch_area_df <-
#  reshape::cast(batch_df, DataFileName ~ Name,  # make wide
#       value = 'TotalPeakArea1', fun.aggregate = sum) #%>%  # vals don't change
#batch_height_df <-
#  reshape::cast(batch_df, DataFileName ~ Name,  # make wide
#                value = 'MajorHeightnA', fun.aggregate = sum) #%>%  # vals don't change
#batch_delta_df <-
#  reshape::cast(batch_df, DataFileName ~ Name,  # make wide
#                value = 'DisplayDelta1', fun.aggregate = sum) #%>%  # vals don't change
#import_funs <- list(
#  import_area <- function(filepath){
#    reshape::cast(batch_df, DataFileName ~ Name,  # make wide
#                  value = 'TotalPeakArea1', fun.aggregate = sum) #%>%  # vals don't change
#  },
#
#  import_height <- function(filepath){
#    reshape::cast(batch_df, DataFileName ~ Name,  # make wide
#                  value = 'MajorHeightnA', fun.aggregate = sum) #%>%  # vals don't change
#  },

# import_delta <- function(filepath){
#    reshape::cast(batch_df, DataFileName ~ Name,  # make wide
#                  value = 'DisplayDelta1', fun.aggregate = sum) #%>%  # vals don't change
#  }

#)

#batch6_list <- lapply(import_funs, function(f){f(df)})
#names(batch6_list) <- c('peak_area', 'peak_height', 'd13c')


########### Function for importing metadata


# Cols in the df used to develop this script:
# ID_With batch
# Name
# RT (Sec)
# Height (nA)
# Corrected 13C
# Peak Name
# Sum Peak Area

#areaw_df <- dcast(batch_narm_df, DataFileName ~ BiomarkerFinal,  # make wide
#                  value.var = 'TotalPeakArea1', fun.aggregate = sum) %>%  # vals don't change
#  mutate(batch = as.character(7)) %>% #str_extract(`ID_With batch`, '[a-zA-Z0-9]\\d+(?=_)')) # Make column for batch ID; may need to get batch some other way, as not all include bath name
#  select(-everything(), batch, everything())

#frac_diff <- function(x){
#  max(x)
#}

#c_names <- names(areaw_df)[3:nrow(areaw_df)]
