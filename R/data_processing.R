#' Subtract surrogate and methylation standard peak areas
#'
#' Gutknecht lab PLFA runs all include a 19;0 surrogate standard spiked
#' into soil samples before extraction, and extracted lipids are suspended
#' in 13:0 control solution. We need to subtract the peak area due to these
#' samples, so that we only measure lipids from the soil sample. This function
#' subtracts 13:0 and 19:0 peak areas to account for this. User must specify
#' the DataFileNames associated with the controls used for peak subtraction.
#' User may also select peaks for subtraction (default is 13:0 and 19:0).
#'
#' Might want to integrate id of blanks etc. in the metadata.
#'
#'@param \code{df} Dataframe or tibble with NormalizedPeakArea. Use the output
#'from \code{normalize_area_tidy()} function.
#'
#'@param \code{blanks}
#'
#'@param \code{lipids}
#'
#'@return
#'
#'@examples
#'
#'
subtract_blanks_tidy <- function(df, blanks = c('100.raw', '101.raw'),
                            lipids = c('13:0', '19:0')){  # version that subtracts avg of blanks only

  blanks_df <- df %>%
    filter(DataFileName %in% blanks) %>%
    group_by(Name) %>%
    mutate(BlankArea = if_else(Name %in% lipids, mean(TotalPeakArea1), 0)) %>%
    select(Name, BlankArea)

  df %>% ungroup() %>%
    left_join(blanks_df, by = 'Name') %>%
    mutate(BlankArea = if_else(is.na(BlankArea), 0, BlankArea)) %>%
    group_by(Batch, Name) %>%
    mutate(AreaMinusBlank = if_else(TotalPeakArea1 - BlankArea < 0, 0,
                                    TotalPeakArea1 - BlankArea))
}

#'
#'
#'
#'
#'
#calc_kval <- function(df, standard_fnames, standard_conc = 250, inj_vol = 2, #standard_fnames should be blanks, not 13:0 standard, I think
#                      standard = '13:0'){
#  kval_df <- df[which(df[['DataFileName']] %in% standard_fnames &
#                        df[['Name']] == standard),] %>%
#    group_by(Name) %>% # took this away, b/c Jess uses single K-value/ whole set
#    mutate_at(vars(TotalPeakArea1), funs(StandardArea = mean(., na.rm = TRUE))) %>%  # intra-batch mean of blanks
#    ungroup() %>%
#    #rename(StandardArea = TotalPeakArea1, StandardBiomarker = Name) %>%  # rename so distinct after join
#    rename(StandardBiomarker = Name) %>%  # rename so distinct after join
#    #group_by(Batch, Biomarker) %>%
#    # currently makes one kval/batch, Jess version of Excel uses single kval for all
#    mutate(kval = StandardArea / !!standard_conc / !!inj_vol) %>%  # calc kval - one/batch
#    #filter(Biomarker == standard & !is.na(Batch)) %>%
#    select(StandardBiomarker, StandardArea, kval)
#
#  kval <- kval_df[['kval']][[1]]
#
#}

#' Convert peak areas to concentration
#'
#' This function converts peak areas to lipid concentration in units of nmol/g
#' dry soil.Peak areas should be normalized and have standards subtracted prior
#' to running this function.
#'
#' @param \code{df}
#' @param \code{standard_fnames}
#' @param \code{mw_df}
#' @param \code{standard_conc}
#' @param \code{innj_vol}
#' @param \code{standard}
#' @param \code{soil_wt_df}
#' @param \code{vial_vol}
#'
#' @return
#'
#' @examples
#'
#'

# Think about what we want the final output to look like - keeping everything is good for troubleshooting,
# but probably not the best in terms of presenting clear results
# Now that the biomarker calcs returns a dataframe, for lipid concentrations, can just stop analysis after
# this point
# An advantage of doing all in one/appending kvals to the df is you can keep df as output and use pipes
calc_concentration_tidy <- function(df, standard_fnames, mw_df = lipid_reference, standard_conc = 250, inj_vol = 2, #standard_fnames should be blanks, not 13:0 standard, I think
                       standard = '13:0', soil_wt_df, vial_vol = 20){
  kval_df <- df[which(df[['DataFileName']] %in% standard_fnames &
                        df[['Name']] == standard),] %>%
    group_by(Batch, Name) %>%
    summarise_at(vars(TotalPeakArea1), mean, na.rm = TRUE) %>%  # intra-batch mean of blanks
    ungroup() %>%
    rename(StandardArea = TotalPeakArea1, StandardBiomarker = Name) %>%  # rename so distinct after join
    #group_by(Batch, Biomarker) %>%
    mutate(kval = (StandardArea / !!standard_conc / !!inj_vol)) %>%  # calc kval - one/batch
    #filter(Biomarker == standard & !is.na(Batch)) %>%
    right_join(df, by = c('Batch')) %>%
    left_join(mw_df, by = c('Name' = 'fame')) %>%
    select(-Batch) %>%
    left_join(samp_wt_df, by = 'DataFileName') %>%
    mutate(nmol_g = (TotalPeakArea1 / kval) * (!!vial_vol / 2) /
             (molecular_weight_g_per_mol * SampleWt))
  return(kval_df)

}

#' Calculate mol fractions of biomarkers and microbial groups
#'
#' This function calculates mole fractions of fungal and bacterial lipids. In
#' the future, this will also calculate moel fractions of Gram (+) and Gram(-)
#' bacteria, actinomycetes, and ___?
#'
#' @param \code{df} dataframe or tibble containing lipid concentrations. Use
#' output from \code{apply_kval()} function.
#'
#' @return
#'
#' @examples
#'

# consider returning a separate summary dataframe for the summary values
# keeping the lipids makes for long df with lots of NA
calc_biomarkers_tidy <- function(df){
  f_lipids <- c('16:1 w5c', '18:1 w9c', '18:2 w6,9c') # consider loading these with the package
  b_lipids <- c('13:0 iso', '13:0 anteiso', '14:0 3OH', '15:0 iso', '15:0 anteiso',
                '16:0 iso', '16:1 w7c', '16:0 10me', '17:0 iso', '17:0 anteiso',
                '18:1 w9t', '18:1 w7c', '18:0 10me')
  gram_pos_lipids <- c('13:0 iso', '13:0 anteiso', '15:0 iso', '15:0 anteiso',
                       '17:iso', '17:0 anteiso')
  gram_neg_lipids <- c('14:0 3OH', '16:1 w7c', '16:1 w9c', '18:1 w7c',
                       '18:1 w9t')
  actino_lipids <- c('16:0 10me', '18:0 10me')
  amf_lipids <- c('16:1 w5c')
  s_fungi_lipids <- c('18:1 w9c', '18:1 w6,9c')
  anaerobe_lipids <- c('19:0 cyclo')

  # Try to find a more efficient way of doing this: perhaps apply, or mutate_at?
  # maybe put indicator lipids in a list of vectors we can loop over
  df <- df %>% filter(!is.na(BatchDataFileName)) %>% # This removes cells added with complete() - no longer needed b/c we're not finding averages; it actually screws up the summary
    group_by(BatchDataFileName, DataFileName, SampleID) %>%
    summarise(TotalBiomass = sum(nmol_g[!is.na(indicates)], na.rm = TRUE),
           TotalFungi = sum(nmol_g[Name %in% !!f_lipids], na.rm = TRUE),
           TotalBact = sum(nmol_g[Name %in% !!b_lipids], na.rm = TRUE),
           GramNeg = sum(nmol_g[Name %in% !!gram_neg_lipids], na.rm = TRUE),
           GramPos = sum(nmol_g[Name %in% !!gram_pos_lipids], na.rm = TRUE),
           Actinomycetes = sum(nmol_g[Name %in% !!actino_lipids], na.rm = TRUE),
           AMF = sum(nmol_g[Name %in% !!amf_lipids], na.rm = TRUE),
           SFungi = sum(nmol_g[Name %in% !!s_fungi_lipids], na.rm = TRUE),
           Anaerobes = sum(nmol_g[Name %in% !!anaerobe_lipids], na.rm = TRUE),

           FungiToBact = TotalFungi / TotalBact,
           PercentFungi = TotalFungi / TotalBiomass * 100,
           PercentBact = TotalBact / TotalBiomass * 100,
           PercentGramNeg = GramNeg / TotalBiomass * 100,
           PercentGramPos = GramPos / TotalBiomass * 100,
           PercentActinomycetes = Actinomycetes / TotalBiomass * 100,
           PercentAMF = AMF / TotalBiomass * 100,
           PercentSFungi = SFungi / TotalBiomass * 100,
           PercentAnaerobes = Anaerobes / TotalBiomass * 100

    )
  return(df)

}

clean_nmol_df_tidy <- function(df){
  clean_df <- df %>%
    mutate(Replicate = str_extract(SampleID, '(?<=[0-9])[A-Za-z]*(?=_)'),
           Year = str_extract(SampleID, '(?<=_)[0-9]+(?=_)'),
           Plot = str_extract(SampleID, '^[0-9]+(?=-)'),
           DepthNum = str_extract(SampleID, '(?<=-)[0-9](?=[A-Za-z])')) %>%
    select(-c(kval, TotalPeakArea1, `Molecular weight (g/mol)`, SampleWt))
  return(clean_df)
}

#' Normalize peak areas across batches
#'
#' GC response varies varies across batches, producing different peak areas
#' for same lipid cocncentration. This function normalizes the peak area
#' across batches for consistent results within an analysis. Current
#' implementation normalizes to the batch with the smallest mean peak area
#' (non-detects not treated as 0 and not included in the average). Averages
#' used for the correction factor only includes experimental samples, blanks
#' (B), and methylation controls (C) (ie. internal standard, M1M2, and FAME).
#' This is consistent with the Excel example calculations provided by Jess.
#'
#' This should only be used VERY carefully and under extreme circumstances, as
#' it introduces some dicey assumptions and may not be totally statistically
#' valid
#'
#' @param \code{df} peak list as dataframe or tibble, loaded using
#' \code{load_batch()} function.
#'
#' @return Returns with same fields as input dataframe as well as
#' Normalized CFactor (conversion factor; used for normalization),
#' NormalizedArea and other intermeidates used in normalization
#' calculations.
#'
#' @examples
#'
#'
# Is there any advantage to keeping the old cols (eg. TotalPeakArea1)
# Need to refactor/optimize this function - got it working, but was focused
# on function over form.
normalize_area_tidy <- function(df){
  cf_numerator_df <- df %>%
    group_by(Batch, DataFileName, BatchDataFileName, Name) %>% #in case dups here on purpose, this merges them
    summarise(TotalPeakArea1 = sum(TotalPeakArea1)) %>%
    ungroup() %>%
    filter(!stringr::str_detect(DataFileName, 'Internal|M1M2|FAME')) %>%
    group_by(BatchDataFileName) %>%
    mutate(AllLipidArea = sum(TotalPeakArea1)) %>%
    ungroup %>%
    group_by(Batch) %>%
    mutate(MeanBatchArea = mean(AllLipidArea)) %>%
    ungroup() %>%
    filter(MeanBatchArea == min(MeanBatchArea)) %>%
    group_by(Name) %>%
    summarise(CFNumerator = mean(TotalPeakArea1))

  df %>% #group_by(Batch, Name) %>%
    #summarise(BatchMeanPeakArea = mean(TotalPeakArea1, na.rm = TRUE)) %>%
    group_by(Batch, DataFileName, BatchDataFileName, Name) %>% #in case dups here on purpose, this merges them
    summarise(TotalPeakArea1 = sum(TotalPeakArea1)) %>%
    ungroup() %>%
    left_join(cf_numerator_df, by = 'Name') %>%
    group_by(Batch, Name) %>%
    mutate(CFactor = CFNumerator / mean(TotalPeakArea1[!stringr::str_detect(DataFileName, 'Internal|M1M2|FAME')])) %>%
    #right_join(df, by = c('Batch', 'Name')) %>%
    mutate(NormalizedArea = TotalPeakArea1 * CFactor)  ### changed 8/15/2019

  #should also add something to make CF values = 1 if there are no samples for denominator
  #had to add some lines to remove teh standards from these calculations
}


#Add a function for displaying data, which will reshape to wide format for easier viewing
# might have to add an argument for which measurement to show

# Add a function for batch importing of batches --> basically the lapply stuff from
# 20190811_spruce_plfa_data-processing.R


# confirmed Blanks have 13:0 only, controls should have 13:0 and 19:0
