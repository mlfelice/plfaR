subtract_blanks_base <- function(df, blanks = c('1.raw', '2.raw'),
                            lipids = c('13:0', '19:0')){  # version that subtracts avg of blanks only

  blanks_df <- df[df[['DataFileName']] %in% blanks, ]
  blanks_df <- aggregate(TotalPeakArea1 ~ Batch + Name,
                         FUN = function(x){
                           blank_area <- mean(x, na.rm = TRUE)
                         }, data = blanks_df
                      )
  blanks_df <- blanks_df[blanks_df[['Name']] %in% lipids, ]

  names(blanks_df)[names(blanks_df) == 'TotalPeakArea1'] <- 'StdArea'

  full_df <- merge(df, blanks_df, by = c('Batch', 'Name'), all.x = TRUE)
  #full_df[is.na(full_df[['StdArea']]), 'StdArea'] <- 0
  full_df <- transform(full_df,
                       StdArea = ifelse(is.na(StdArea), 0, StdArea))
  full_df <- transform(full_df,
                       AreaMinusBlanks = ifelse(TotalPeakArea1 - StdArea > 0,
                                                TotalPeakArea1 - StdArea, 0))

  return(full_df)
}

area_to_concentration_base <- function(df, standard_fnames, soil_wt_df,
                                       mw_df = lipid_reference,
                                       standard_conc = standard_conc, inj_vol = inj_vol, #standard_fnames should be blanks, not 13:0 standard, I think
                                       standard = standard, vial_vol = vial_vol){



  standard_vec <- unlist(df[df[['DataFileName']] %in% standard_fnames &
                             df[['Name']] == standard, 'TotalPeakArea1'])

  kval <- mean(standard_vec) / standard_conc / inj_vol

  tmp_df <- merge(df, soil_wt_df, by = c('Batch', 'DataFileName'), all.x = TRUE)
  tmp_df <- merge(tmp_df, mw_df, by.x = 'Name', by.y = 'fame')
  tmp_df <- transform(tmp_df, nmol_g = (AreaMinusBlanks / kval) *
                         (vial_vol / 2) /
                         (molecular_weight_g_per_mol * SampleWt))

  return(tmp_df)

}

process_peak_area_base <- function(dat, standard_fnames, mw_df = lipid_reference,
                              standard_conc = 250, inj_vol = 2, #standard_fnames should be blanks, not 13:0 standard, I think
                              standard = '13:0', soil_wt_df, vial_vol = 50,
                              lipids = c('13:0', '19:0')){

  # could also structure so that user can input 'blanks' arg that would
  # supersede the SampleType from metadata
  blanks <- unlist(soil_wt_df[grepl('[Bb]lank', soil_wt_df[['SampleType']]),
                              'DataFileName'])
  dtype <- check_format(dat)

  if (dtype == 'data.frame') {
    tmp_df <- subtract_blanks_base(dat, blanks = blanks,
                    lipids = c('13:0', '19:0'))
    area_to_concentration_base(tmp_df, standard_fnames, mw_df = lipid_reference,
                      standard_conc = standard_conc, inj_vol = inj_vol, #standard_fnames should be blanks, not 13:0 standard, I think
                      standard = standard, soil_wt_df, vial_vol = vial_vol)
  }
}

###########
# TO DO: output info on the warnings that are thrown if there are duplicates
# during reshape
# An alternative to this approach is to calculate mol% for each lipid and just
# attach the name of the group it represents, so that final dataframe has both
# individual lipid data and the data associated with a microbial group. This
# looks like the approach that Jess/Cameron took. The advantage would be
# keeping all data in one df, and simplified calculations.

calculate_indicators_base <- function(df, soil_wt_df){

  indicator_list <- list(f_lipids = c('16:1 w5c', '18:1 w9c', '18:2 w6,9c'),
                         b_lipids = c('13:0 iso', '13:0 anteiso', '14:0 3OH',
                                      '15:0 iso', '15:0 anteiso', '16:0 iso',
                                      '16:1 w7c', '16:0 10me', '17:0 iso',
                                      '17:0 anteiso', '18:1 w9t', '18:1 w7c',
                                      '18:0 10me'),
                         gram_pos_lipids = c('13:0 iso', '13:0 anteiso',
                                             '15:0 iso', '15:0 anteiso',
                                             '17:iso', '17:0 anteiso'),
                         gram_neg_lipids = c('14:0 3OH', '16:1 w7c',
                                             '16:1 w9c', '18:1 w7c',
                                             '18:1 w9t'),
                         actino_lipids = c('16:0 10me', '18:0 10me'),
                         amf_lipids = c('16:1 w5c'),
                         s_fungi_lipids = c('18:1 w9c', '18:1 w6,9c'),
                         anaerobe_lipids = c('19:0 cyclo'),
                         # May need to redefine the total_biomass. Right now, it's just all of the
                         # indicator lipids. Not sure if there should be others as well
                         total_biomass = c("16:1 w5c", "18:1 w9c", "18:2 w6,9c",
                                           "13:0 iso", "13:0 anteiso",
                                           "14:0 3OH", "15:0 iso",
                                           "15:0 anteiso", "16:0 iso",
                                           "16:1 w7c", "16:0 10me", "17:0 iso",
                                           "17:0 anteiso", "18:1 w9t",
                                           "18:1 w7c", "18:0 10me", "17:iso",
                                           "16:1 w9c", "18:1 w6,9c",
                                           "19:0 cyclo")
  )

  df_wide <- reshape(data = df[c('DataFileName', 'Name', 'nmol_g')],
                     timevar = 'Name', idvar = 'DataFileName',
                     direction = 'wide')
  names(df_wide) <- gsub('nmol_g.', '', names(df_wide))
  # The reshape step will produce NA for any lipids that were not detected.
  # Some entries will have 0 due to subtracting the standards (only)

  df_wide[is.na(df_wide)] <- 0 # replace NA with 0 for calculations

  # without drop = FALSE, subsetting returns vector and rowSums() errors
  total_biomass <- rowSums(df_wide[,
                                   names(df_wide) %in% indicator_list[['total_biomass']],
                                   drop = FALSE],
                           na.rm = TRUE) # need to limit this to only microbial lipids

  # Convert biomarker concentrations to % of total biomass
  # can use this code chunk if we go the route of calcs for indiv biomarkers
  #perc_df <- cbind(df_wide['DataFileName'], total_biomass,
  #                 lapply(df_wide[, 2:(ncol(df_wide))],
  #                        function(x){x/total_biomass * 100}))

  # aggregate lipid nmol_g by microbial group (sum)
  nm_df <- cbind(df_wide['DataFileName'],
                 lapply(indicator_list,
                        function(x){
                          tmp_df <- df_wide[,
                                            names(df_wide) %in% x,
                                            drop = FALSE]
                          rowSums(tmp_df, na.rm = TRUE)
                        }
                 )
  )

  perc_df <- cbind(nm_df['DataFileName'],
                   lapply(names(indicator_list),
                          function(x){ nm_df[x] / total_biomass * 100}
                   )
  )

  # NaN's produced from 0/0 (if there was no total_biomass and no detected lipids for an indicator)

  perc_df['f_to_b'] <- nm_df[['f_lipids']] / nm_df[['b_lipids']] * 100
  #f_to_b is also appended as a column now, will be treated like other indicators if we reshape

  ncp <- ncol(perc_df)  # define this based on ncol so doesn't break if we change
  # biomarkers
  perc_df_long <- reshape(perc_df, varying = names(perc_df)[c(2:ncp)],
                          v.names = 'Percent',
                          timevar = 'Indicator', idvar = 'DataFileName',
                          times = names(perc_df)[c(2:ncp)], direction = 'long')

  ncn <- ncol(nm_df)  # define this based on ncol so doesn't break if we change
  # biomarkers
  nm_df_long <- reshape(nm_df, varying = names(nm_df)[c(2:ncn)],
                        v.names = 'nmol_g',
                        timevar = 'Indicator', idvar = 'DataFileName',
                        times = names(nm_df)[c(2:ncn)], direction = 'long')

  # combine the mol%, absolute abundance, and metadata
  perc_df_long <- merge(perc_df_long, soil_wt_df,
                        by = 'DataFileName', all.x = TRUE)
  perc_df_long <- merge(perc_df_long, nm_df_long,
                        by = c('DataFileName', 'Indicator'), all.x = TRUE)

}

#################################################
# This is a little rough and not yet functioning
# Try to get it right later
#calc_biomarkers_base <- function(df){
#  indicator_list <- list(f_lipids <- c('16:1 w5c', '18:1 w9c', '18:2 w6,9c'),

#                         b_lipids <- c('13:0 iso', '13:0 anteiso', '14:0 3OH', '15:0 iso', '15:0 anteiso',
#                                       '16:0 iso', '16:1 w7c', '16:0 10me', '17:0 iso', '17:0 anteiso',
#                                       '18:1 w9t', '18:1 w7c', '18:0 10me'),

#                         gram_pos_lipids <- c('13:0 iso', '13:0 anteiso', '15:0 iso', '15:0 anteiso',
#                                              '17:iso', '17:0 anteiso'),
#
#                         gram_neg_lipids <- c('14:0 3OH', '16:1 w7c', '16:1 w9c', '18:1 w7c',
#                                              '18:1 w9t'),

#                         actino_lipids <- c('16:0 10me', '18:0 10me'),

                         #amf_lipids <- c('16:1 w5c'),

#                         #s_fungi_lipids <- c('18:1 w9c', '18:1 w6,9c'),
#
#                         anaerobe_lipids <- c('19:0 cyclo'))

#  # Try to find a more efficient way of doing this: perhaps apply, or mutate_at?
#  df <- df[!is.na(df[['BatchDataFileName']]), ]  # This removes cells added with complete() - no longer needed b/c we're not finding averages; it actually screws up the summary
#
#  tot_biomass_df <- aggregate(nmol_g ~ Batch + DataFileName,
#                              FUN = function(x){
#                                TotalBiomass <- sum(x, na.rm = TRUE)
#                              }, data = df[!is.na(df[['indicates']]), ]
#  )


#  tmp_df <- lapply(indicator_list,
#         function(x){print(x)
#           aggregate(nmol_g ~ Batch + DataFileName,
#                     FUN = function(y){
#                       sum(y, na.rm = TRUE)
#                     },
#                     data = df[df[['Name']] %in% x, ])
#         }
#  )



#              FungiToBact = TotalFungi / TotalBact
#              PercentFungi = TotalFungi / TotalBiomass * 100
#              PercentBact = TotalBact / TotalBiomass * 100
#              PercentGramNeg = GramNeg / TotalBiomass * 100
#              PercentGramPos = GramPos / TotalBiomass * 100
#              PercentActinomycetes = Actinomycetes / TotalBiomass * 100
#              PercentAMF = AMF / TotalBiomass * 100
#              PercentSFungi = SFungi / TotalBiomass * 100
#              PercentAnaerobes = Anaerobes / TotalBiomass * 100
#

#  return(df)
#
#}
