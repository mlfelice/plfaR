# Current development:
# Write code to import data from multiple batches, append batch #, and merge into single dataframe for
# remaining analysis

######
# 2nd interation integrating the quality checks - not the neatest, but seems to work

# 1. Load packages ####
#----------------------
library(tidyverse)
library(readxl)
library(writexl)
library(reshape2)
library(stringr)

# 2. Define functions ####
# ------------------------
id_dups <- function(df, start_col = 1){
  #
  # function checks and outputs biomarkers with duplicate entries within sample
  #
  # Args:
  #   df: dataframe of count of peaks
  #       - wide format (cols = biomarkers, rows = samples)
  #       - sample names can either be rownames or a column named 'biomarker'
  #   start_col: the first column containing biomarkers
  #
  # Returns:
  #   dataframe of biomarkers and samples with duplicate biomarkers (if any)
  #-------------------------------------------------------------------------

  bm_col <- grepl('[Bb]iom', colnames(df))

  if(any(bm_col) == TRUE){  # converts col to rowname if biomarks are a col
    df <- column_to_rownames(df, colnames(df[bm_col]))
  }

  # store indices of col and row w/relevant data
  i <- start_col
  j <- ncol(df)

  if(any(df[, i:j] > 1) == TRUE){

    # get rownames/colnames of entries > 1 (duplicate name x sample)
    files_dup <- colnames(df[, i:j])[col(df[, i:j])[which(df[, i:j] > 1)]]
    bms_dup <- rownames(df[, i:j])[row(df[, i:j])[which(df[, i:j] > 1)]]
    # dataframe of sample x biomarkers w/duplicate
    dup_df <- tibble(DataFileName = files_dup, Biomarker = bms_dup)

  } else{print('No duplicate peak names identified')}

  return(dup_df)

}
# old version 8/6/2019
#id_dups <- function(df){
  #
  # function checks and outputs biomarkers with duplicate entries within sample
  #
  # Args:
  #   df: dataframe of count of peaks
  #       - wide format (cols = biomarkers, rows = samples)
  #       - sample names can either be rownames or a column named 'biomarker'
  #   start_col: the first column containing biomarkers
  #
  # Returns:
  #   dataframe of biomarkers and samples with duplicate biomarkers (if any)
  #-------------------------------------------------------------------------
#  files_col <- grepl('DataFileName', colnames(df))

#  if(any(files_col) == TRUE){  # converts col to rowname if biomarks are a col
#    df <- column_to_rownames(df, colnames(df[files_col]))
#  }


#  if(any(df > 1) == TRUE){

    # get rownames/colnames of entries > 1 (duplicate name x sample)
#    bms_dup <- colnames(df)[col(df)[which(df > 1)]]
#    files_dup <- rownames(df)[row(df)[which(df > 1)]]
    # dataframe of sample x biomarkers w/duplicate
#    dup_df <- tibble(DataFileName = files_dup, Biomarker = bms_dup)

#  } else{print('No duplicate peak names identified')}

#  return(dup_df)

#}
id_dups <- function(df){
  #
  # function checks and outputs biomarkers with duplicate entries within sample
  #
  # Args:
  #   df: dataframe of count of peaks
  #       - wide format (cols = biomarkers, rows = samples)
  #       - sample names can either be rownames or a column named 'biomarker'
  #   start_col: the first column containing biomarkers
  #
  # Returns:
  #   dataframe of biomarkers and samples with duplicate biomarkers (if any)
  #-------------------------------------------------------------------------
  files_col <- grepl('DataFileName', colnames(df))

  if(any(files_col) == TRUE){  # converts col to rowname if biomarks are a col
    df <- column_to_rownames(df, colnames(df[files_col]))
  }


  if(any(df > 1) == TRUE){

    # get rownames/colnames of entries > 1 (duplicate name x sample)
    bms_dup <- colnames(df)[col(df)[which(df > 1)]]
    files_dup <- rownames(df)[row(df)[which(df > 1)]]
    # dataframe of sample x biomarkers w/duplicate
    dup_df <- tibble(DataFileName = files_dup, Biomarker = bms_dup)

  } else{print('No duplicate peak names identified')}

  return(dup_df)

}

#find_miss <- function(df, lipids = c('13:0', '16:0', '19:0')){
  #
  # Returns lists for each lipid of samples without that lipid
  #
  # Args:
  #   df:
  #   lipids:
  #
  # Returns:
  #   returns a list of standard? lipids containing vector of samples
  #      not containing it
  #
  #------------------------------------------------------------------
#  no_lipid <- lapply(lipids, function(x){df[which(df[x] == 0), 1]})
#  names(no_lipid) <- lipids

#  return(no_lipid)

#}
#find_miss <- function(df, lipids = c('13:0', '16:0', '19:0')){
  #
  # function checks and outputs biomarkers with duplicate entries within sample
  #
  # Args:
  #   df: dataframe of count of peaks
  #       - wide format (cols = biomarkers, rows = samples)
  #       - sample names can either be rownames or a column named 'biomarker'
  #   start_col: the first column containing biomarkers
  #
  # Returns:
  #   dataframe of biomarkers and samples with duplicate biomarkers (if any)
  #-------------------------------------------------------------------------

#  files_col <- grepl('DataFileName', colnames(df))

#  if(any(files_col) == TRUE){  # converts col to rowname if biomarks are a col
#    df <- column_to_rownames(df, colnames(df[files_col]))
#  }

#  df <- df[, lipids]


#  if(any(df == 0) == TRUE){

    # get rownames/colnames of entries > 1 (duplicate name x sample)
#    bms_miss <- colnames(df)[col(df)[which(df == 0)]]
#    files_miss <- rownames(df)[row(df)[which(df == 0)]]
    # dataframe of sample x biomarkers w/duplicate
#    miss_df <- tibble(DataFileName = files_miss, Biomarker = bms_miss)

#  } else{print('No missing standards identified')}

#  return(miss_df)

#}
find_miss <- function(df, lipids = c('13:0', '16:0', '19:0')){
  #
  # function checks and outputs biomarkers with duplicate entries within sample
  #
  # Args:
  #   df: dataframe of count of peaks
  #       - wide format (cols = biomarkers, rows = samples)
  #       - sample names can either be rownames or a column named 'biomarker'
  #   start_col: the first column containing biomarkers
  #
  # Returns:
  #   dataframe of biomarkers and samples with duplicate biomarkers (if any)
  #-------------------------------------------------------------------------

  files_col <- grepl('DataFileName', colnames(df))

  if(any(files_col) == TRUE){  # converts col to rowname if biomarks are a col
    df <- column_to_rownames(df, colnames(df[files_col]))
  }

  df <- df[, lipids]


  if(any(df == 0) == TRUE){

    # get rownames/colnames of entries > 1 (duplicate name x sample)
    bms_dup <- colnames(df)[col(df)[which(df == 0)]]
    files_dup <- rownames(df)[row(df)[which(df == 0)]]
    # dataframe of sample x biomarkers w/duplicate
    miss_df <- tibble(DataFileName = files_dup, Biomarker = bms_dup)

  } else{print('No missing standards identified')}

  return(miss_df)

}

#count_lips <- function(df, start_col = 2){
  #
  # determines detection freqeuncy for each lipid in input dataframe
  #
  # Args:
  #   df: dataframe with # of peaks for each biomarker (cols) and sample (rows)
  #   start_col: first column containing biomarker
  #
  # Returns:
  #   dataframe listing the detection frequency for each biomarker (rows)

#  df_trim <- df[,start_col:ncol(df)] # Subset dataframe (remove sample name column)
#  column_names <- vector(length = ncol(df_trim)) # Empty vector to house sample names
#  lipid_freq <- vector(length = ncol(df_trim)) # Empty vector to house freq vals

#  for(i in 1:ncol(df_trim)){ # Loop through all cols to calc frequency of detection
#    column_names[i] <- names(df_trim)[i]
#    lipid_freq[i] <- sum(df_trim[[i]])/length(df_trim[[i]])
#  }

#  lipid_freq_df <- tibble(Lipid = column_names, Detection_freq = lipid_freq) # Make dataframe from vecs

#  print(lipid_freq_df)
#}
count_lips <- function(df, start_col = 2){
  #
  # determines detection freqeuncy for each lipid in input dataframe
  #
  # Args:
  #   df: dataframe with # of peaks for each biomarker (cols) and sample (rows)
  #       - should only contain biomarker frequencies and filenames
  #
  # Returns:
  #   dataframe listing the detection frequency for each biomarker (rows)

  files_col <- grepl('DataFileName', colnames(df))

  if(any(files_col) == TRUE){  # converts col to rowname if biomarks are a col
    df <- column_to_rownames(df, colnames(df[files_col]))
  }

  freq <- sapply(df, function(x){sum(x)/length(x)})

  lipid_freq_df <- tibble(Lipid = names(df), Detection_freq = freq) # Make dataframe from vecs

  return(lipid_freq_df)
}

# Transform values in the batch dataframe (a) by multiplying by the correction factor/fractional difference from max
### Deprecated version
normalize_area <- function(df, frac_area_df, batches, start_col = 2){

  i = start_col

  for(b in batches){
    cf <- frac_area_df[!is.na(frac_area_df$Batch) & frac_area_df$Batch == b,
                       ][[i]]

    normalized_df <- df

    normalized_df[!is.na(normalized_df$Batch) & normalized_df$Batch == b,
                  ][[i]] <- normalized_df[!is.na(normalized_df$Batch) &
                                            normalized_df$Batch == b,
                                          ][[i]] * cf
  }
  return(normalized_df)
}
###
normalize_area <- function(df, biomarkers){

  cf_long_df <- df %>% group_by(Batch) %>%
    summarise_at(biomarkers, funs(mean), na.rm = TRUE) %>%
    mutate_at(vars(biomarkers), function(x) {max(x)/x}) %>%
    gather(key = 'Biomarker', value = 'ConversionFactor',
           c('??':'8:0??'))

  area_long_df <- gather(data = df, key = 'Biomarker',
                         value = 'TotalPeakArea1', c('??':'8:0??'))

  normalized_df <- area_long_df %>%
    left_join(cf_long_df, by = c('Batch', 'Biomarker')) %>%
    mutate(NormalizedArea = ConversionFactor * TotalPeakArea1) %>%
    select(-TotalPeakArea1, -ConversionFactor) %>%
    spread(key = Biomarker, value = NormalizedArea)

  return(normalized_df)
}

subt_stand <- function(df, stand_df){  # deprecated version

  for(s in unique(stand_df[['biomarker']])){
    #search_string <- paste0(str_extract(s, '[0-9]+'), '.[Ss]tandard')
    fnames_df <- subset(stand_df, biomarker == s)
    stand_vec <- df[which(df$DataFileName %in% fnames_df[['name']]), s]
    #stand_vec <- df[which(str_detect(df[ ,2], ls1[ls1$biomarker == s, 2])),
    #                which(names(df) == s)]
    stand_val <- mean(stand_vec)

    df[, which(names(df) == s)] <-
      df[, which(names(df) == s)] - stand_val
    df[df[[s]] < 0, which(names(df) == s)] <- 0
  }
  return(df)
}

subt_stand <- function(df, standard_fnames){  # Latest version

  df2 <- df[which(df[['DataFileName']] %in% standard_fnames &
                    !is.na(df[['BiomarkerFinal']])),] %>%
    group_by(Batch, BiomarkerFinal) %>%
    summarise_at(vars(TotalPeakArea1), mean, na.rm = TRUE) %>%
    rename(StandardArea = TotalPeakArea1) %>%
    ungroup() %>%
    right_join(df, by = c('Batch', 'BiomarkerFinal')) %>%
    mutate(TotalPeakArea1 - StandardArea)

  df2[, df2[['TotalPeakArea']] < 0] <- 0

  return(df2)
}

get_mw <- function(df, marker){
  # Returns the molecular weight of specified biomarker
  #
  # Args:
  #   df: reference dataframe containing biomarker names and molecular weight
  #   marker: the biomarker whose molecular weight you want to look up
  #
  # Returns:
  #   molecular weight of specified biomarker (numeric)
  #
  df[which(df$`FAME ID` == marker), 2][[1]]
}

apply_k <- function(df, stand_df, mwt_df, standard_conc = 250, inj_vol = 2,
                    standard = '13:0', soil_wt_df, vial_vol = 20,
                    start_col = 4){
  # Converts the peak area to lipid concentration
  #
  # Args:
  #   df: dataframe of normalized peak areas with 13:0 and 19:0 subtracted
  #   stand_df: dataframe identifying filenames for standards
  #   mwt_df: reference dataframe containing molecular weights of biomarkers
  #   standard_conc: concentration (nmol/uL) of specified standard used
  #   inj_vol: volume of specified standard injected
  #   standard: lipid standard used for area to conc calculation
  #   soil_wt_df: dataframe containing the recorde weights of peat
  #   vial_vol: total volume in GC vial
  #   start_col: index of column with first biomarker
  #
  # Returns:
  #   dataframe containing total biomass
  temp_df <- df[, start_col:ncol(df)]  # Remove sample names from dataframe

  fnames_df <- subset(stand_df, biomarker == standard)
  stand_vec <- df[which(df$DataFileName %in% fnames_df[['name']]), standard]
  stand_val <- mean(stand_vec)
  kval <- stand_val / standard_conc / inj_vol

  for(i in mwt_df[[1]]){  # First column of reference dataframe (biomarker names)
    mw <- get_mw(mwt_df, i)
    df[i] <- (df[[i]] / kval) * (vial_vol / 2) / (mw)
  }

  df_long <- df %>% gather(Biomarker, Concentration,
                    names(df)[start_col]:names(df)[length(names(df))]) %>%
    select(-Batch) %>%
    left_join(soil_wt_df) %>%
    mutate(nmol_per_g_soil = Concentration / SampleWt)

  df <- df_long %>%
    select(-Concentration) %>%
    spread(Biomarker, nmol_per_g_soil)


  cat('No reference molecular weight for: \n\n')
  print(
    names(df[,start_col:ncol(df)])[which(!names(df[, start_col:ncol(df)]) %in%
                                   mwt_df[[1]])]
  )

  df <- df %>% mutate(
    #total_biomass = select_(.dots = match(mwt_df[[1]], names(df))) %>% rowSums) %>%
    total_biomass = rowSums(select(df, which(names(.) %in% mwt_df[[1]])),
                            na.rm = TRUE)) #%>%
  #select(batch, DataFileName, total_biomass)

  return(df)

}

o_wd <- getwd() #store original wd for restore

setwd('C:/Users/Mark/Dropbox/umn_gutknecht_postdoc/spruce_project/plfa_13c/')

# Load sample wght data sheet
samp_wt_df <- read_xlsx(path = '20190705_spruce_plfa_sample-prep-wksht.xlsx',
                        sheet = 'SPRUCE batches ', na = 'NA') %>%  # sht name incl trailing space
  fill(`GC Batch #`) %>%
  #select(`Sample ID`, `Renamed on day 4`, `GC Batch #`) %>%
  mutate(SampleID = paste(`Sample ID`, '2016',
                          str_extract(string = `GC Batch #`,
                                      pattern = '[0-9]+$'), sep = '_')) %>%
  mutate(DataFileName = ifelse((!is.na(`Renamed on day 4`) &
                                  grepl('^[0-9]+', `Renamed on day 4`) == TRUE),
                               paste0(`Renamed on day 4`, '.raw'), NA),
         Batch = str_extract(string = `GC Batch #`, pattern = '[0-9]$+')) %>%
  rename(SampleWt = `Sample weight`) %>%
  select(Batch, DataFileName, SampleID, SampleWt)


# 3. Load raw data for multiple batch files in single directory ####
#-------------------------------------------------------------------

source_dir <- choose.dir()

batch_files <- list.files(path = source_dir, full.names = TRUE)

batch_list <- vector('list', length = 5)
names(batch_list) <- c('named_peaks', 'qc_df', 'dup_lipids', 'missing_stds',
                       'lipid_freq')
for(n in names(batch_list)){
  batch_list[[n]] <- vector('list', length = length(batch_files))
}

i = 1
for(bf in batch_files){
  batch_name <- str_extract(string = bf, pattern = '[Bb]atch ?[0-9]+')  # ? looks for 0 or 1

  batch_list[['named_peaks']][[i]] <- read_xlsx(bf, sheet = 'named_peaks', na = 'NA') %>%
    select(-BiomarkerRTBased, -Notes) %>%
    filter(!is.na(BiomarkerFinal) & BiomarkerFinal != 'Check chromatogram for 18 peaks' & BiomarkerFinal != 'nothing')
  batch_list[['named_peaks']][[i]]['Batch'] <- batch_name
  batch_list[['named_peaks']][[i]] <- batch_list[['named_peaks']][[i]] %>%
    mutate(BatchDataFileName = paste(Batch, DataFileName, sep = '_')) %>%
    select(Batch, DataFileName, BatchDataFileName, everything())
  names(batch_list[['named_peaks']])[[i]] <- batch_name
  batch_list[['named_peaks']][[i]]['DisplayDelta1'] <- as.character(batch_list[['named_peaks']][[i]]['DisplayDelta1'])
###
  batch_list[['qc_df']][[i]] <- dcast(batch_list[['named_peaks']][[i]], DataFileName ~ BiomarkerFinal,
                 value.var = 'TotalPeakArea1',
                 fun.aggregate = length)  # Make wide and aggregate
  names(batch_list[['qc_df']])[[i]] <- batch_name

  #### Aggregate data - counts of biomarker (rows) in ID'd in each sample (cols) ####
#  qc_df2 <- dcast(batch_list[['named_peaks']][[i]], BiomarkerFinal ~ DataFileName,
#                  value.var = 'TotalPeakArea1',
#                  fun.aggregate = length)  # Make wide and aggregate
  qc_df <- dcast(batch_list[['named_peaks']][[i]], DataFileName ~ BiomarkerFinal,
                 value.var = 'TotalPeakArea1',
                 fun.aggregate = length)  # Make wide and aggregate

  batch_list[['dup_lipids']][[i]] <- id_dups(qc_df)
  names(batch_list[['dup_lipids']])[[i]] <- batch_name



  # check to see that 16:O, 13:O, & 19:O in each sample
  batch_list[['missing_stds']][[i]] <- find_miss(qc_df)
  names(batch_list[['missing_stds']])[[i]] <- batch_name

  # sum each column and divide by nrow to determine the percentage of samples with a named peak
  batch_list[['lipid_freq']][[i]] <- count_lips(qc_df)
  names(batch_list[['lipid_freq']])[[i]] <- batch_name

  i = i + 1
}

all_batch_df <- bind_rows(batch_list[['named_peaks']])

areaw_df <- dcast(all_batch_df, BatchDataFileName ~ BiomarkerFinal,  # make wide
                  value.var = 'TotalPeakArea1', fun.aggregate = sum) %>%  # vals don't change
  mutate(Batch = str_extract(BatchDataFileName, '[Bb]atch ?[0-9]+'),
         DataFileName = gsub('[Bb]atch ?[0-9]+_', '', BatchDataFileName)) %>%
  select(Batch, DataFileName, everything())

c_names <- names(areaw_df)[4:length(names(areaw_df))]

# Make dataframe with batch mean peak area for all biomarkers
area_frac_df <- areaw_df %>% group_by(Batch) %>%
  summarise_at(c_names, funs(mean), na.rm = TRUE)# %>%
#summarise_at(c_names, funs(frac_diff))

# Calculate fractional difference (relative to max)
area_frac_df <- bind_cols(area_frac_df[,1], # Create dataframe of correction factors
                          as_tibble(lapply(area_frac_df[,2:ncol(area_frac_df)],
                                           function(x){max(x)/x})))




normalized_area_df = normalize_area(df = areaw_df, frac_area_df = area_frac_df,
                                    batches = c('Batch4', 'Batch5', 'Batch6', 'Batch7'), start_col = 4)


#batches <- c('B14')
#i = 2
#for(b in batches){
#  cf <- area_frac_df[!is.na(area_frac_df$batch) & area_frac_df$batch == b, ][[i]]
#
#  a[!is.na(a$batch) & a$batch == b, ][[i]] <- a[!is.na(a$batch) & a$batch == b, ][[i]] * cf

# Subtract values of blanks 13:0 and 19:0 if difference > 0 ---> Not totally sure I understand this step


# Make dataframe containing biomarkers and associated file names
ls1 <- data.frame(biomarker = c('13:0', '13:0'),
                  name = c('Internal std 1.raw', 'Internal std 2.raw'),
                  stringsAsFactors = FALSE)  # Make df w/samp names of markers

stands_sub_df <- subt_stand(normalized_area_df, stand_df = ls1)
#Asdflkj;sldakf;sdlaf;lkasdjf;dlkjas;lkj;dlkja;slkdjfasdfkdsljlsadk-80characters
# Calculate K-value and apply it to peak areas (this converts to nmol/__, etc)
mwt_ref_path <- 'C:/Users/Mark/Dropbox/umn_gutknecht_postdoc/spruce_project/plfa_13c/20190521_spruce_plfa_biomarker-molecular-wt.xlsx'
mwt_ref_df <- read_xlsx(mwt_ref_path)

mwt_ref_df <- mwt_ref_df[2:nrow(mwt_ref_df), ]

nmol_df <- apply_k(df = normalized_area_df, stand_df = ls1, , mwt_df = mwt_ref_df,
                   soil_wt_df = samp_wt_df, start_col = 4)

# Calculate % of total biomass and append as final column
#nmol_df <- nmol_df %>% mutate(
#  total_biomass = rowSums(select(., -batch, -DataFileName), na.rm = TRUE)) %>%
#  select(batch, DataFileName, total_biomass, `10:0`:`8:0`)

# Convert biomarker concentrations to % of total biomass
perc_df <- bind_cols(nmol_df['Batch'], nmol_df['DataFileName'],
                     lapply(nmol_df[, 6:(ncol(nmol_df))],
                            function(x){x/nmol_df[['total_biomass']] * 100}))

# Calculate fungal:bacterial ratio  -- This needs work
f_lipids <- c('16:1 w5c', '18:1 w9c', '18:2 w6,9c')
b_lipids <- c('13:0 iso', '13:0 anteiso', '14:0 3OH', '15:0 iso', '15:0 anteiso',
              '16:0 iso', '16:1 w7c', '16:0 10me', '17:0 iso', '17:0 anteiso',
              '18:1 w9t', '18:1 w7c', '18:0 10me')

# Calculate and add percent fungi,bacteria, and fungal:bacterial ratio
perc_df['fungi_conc'] <- rowSums(nmol_df[, f_lipids], na.rm = TRUE)
perc_df['bact_conc'] <- rowSums(nmol_df[,b_lipids[b_lipids %in% names(nmol_df)]],
                                na.rm = TRUE)
perc_df['fb'] <- perc_df['fungi_conc']/perc_df['bact_conc']
perc_df['fungi_perc'] <- perc_df['fungi_conc']/nmol_df['total_biomass'] * 100
perc_df['bact_perc'] <- perc_df['bact_conc']/nmol_df['total_biomass'] * 100


#####


# messing around, this is my own version of str_extract()
find_string <- function(pattern, text){
  match_result <- regexpr(pattern = pattern, text = text)
  start_match <- match_result[1]
  end_match <- match_result[1] + attributes(match_result)$match.length
  string_match <- substring(text, first = start_match, last = end_match - 1)
  return(string_match)
}
find_string(pattern = '[Bb]atch.?[0-9]+', text = 'SPRUCE_plfa_batch9.raw')

?select


######
######
##### first iteration combining multiple batches into list/same dataframe
library(stringr)
source_dir <- choose.dir()

batch_files <- list.files(path = source_dir, full.names = TRUE)

batch_list <- vector('list', length = length(batch_files))
i = 1
for(bf in batch_files){
  batch_list[[i]] <- read_xlsx(bf, sheet = 'named_peaks', na = 'NA') %>%
    select(-BiomarkerRTBased, -Notes)

  batch_name <- str_extract(string = bf, pattern = '[Bb]atch ?[0-9]+')  # ? looks for 0 or 1

  batch_list[[i]]['Batch'] <- batch_name

  names(batch_list)[[i]] <- batch_name

  i = i + 1
}

all_batch_df <- bind_rows(batch_list)

######
######

#### Load Jess's template to check our calcs ####
source_path <- file.choose()

test_plfa_df <- batch_raw_df <- read_xls(source_path,
                                         sheet = 'All Named Data', na = 'NA',
                                         range = cell_cols('A:F'),
                                         col_types = c('text', 'numeric',
                                                       'numeric', 'numeric',
                                                       'text', 'text')) %>%
  rename(BatchDataFileName = RowLabels, BiomarkerFinal = Names) %>%
  mutate(Batch = str_extract(string = BatchDataFileName, pattern = '[Bb]atch ?[0-9]+'),
         DataFileName = gsub('[Bb]atch ?[0-9]+_', '', BatchDataFileName),
         PN = gsub('[Bb]atch ?[0-9]+_', '', BatchDataFileName)) %>%
  select(Batch, DataFileName, BatchDataFileName, everything())

test_plfa_df[test_plfa_df$DisplayDelta1 == '#Type mismatch', 'DisplayDelta1'] <- NA
test_plfa_df$DisplayDelta1 <- as.numeric(test_plfa_df$DisplayDelta1)

filter_df <- c('NB[0-9]', '_[0-9]', 'FAME', 'Hex', 'M1', '[A-Za-z][0-9]\\.[0-9]',
               '_C2', 'intern', 'PB1')

test_p_df <- test_plfa_df %>%
  filter(grepl(paste(filter_df, collapse = '|'), DataFileName))

test_p_df[test_p_df$DataFileName == 'Batch2_N2.2N.raw',]
test_plfa_df[grepl('[A-Za-z][0-9]\\.[0-9]', test_plfa_df$DataFileName),]

grepl('_[A-Za-z][0-9]\\.[0-9]', 'Batch2_N2.2N.raw')

####
# Need to modify apply_k() to use actual soil weights associated with each sample
# Best solution is probably to do whole nmol calc (minus the soil weight component)
# Then join soil weights df to the result and use mutate to apply the soil weight


o_wd <- getwd() #store original wd for restore

setwd('C:/Users/Mark/Dropbox/umn_gutknecht_postdoc/spruce_project/plfa_13c/')


df_long <- nmol_df %>% gather(Biomarker, Concentration, names(nmol_df)[4]:names(nmol_df)[length(names(nmol_df))]) %>%
  select(-Batch) %>%
  left_join(samp_wt_df, by = 'DataFileName') %>%
  mutate(nmol_per_g_soil = Concentration / SampleWt)
df <- df_long %>%
  select(-Concentration) %>%
  spread(Biomarker, nmol_per_g_soil)
