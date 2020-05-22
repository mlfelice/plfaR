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
# This alternate approach could also be used for the d13C

calculate_indicators_base <- function(df, soil_wt_df){

  indicator_list <- list(f_lipids <- c('16:1 w5c', '18:1 w9c', '18:2 w6,9c'),
                         b_lipids <- c('13:0 iso', '13:0 anteiso', '14:0 3OH', '15:0 iso', '15:0 anteiso',
                                       '16:0 iso', '16:1 w7c', '16:0 10me', '17:0 iso', '17:0 anteiso',
                                       '18:1 w9t', '18:1 w7c', '18:0 10me'),
                         gram_pos_lipids <- c('13:0 iso', '13:0 anteiso', '15:0 iso', '15:0 anteiso',
                                              '17:iso', '17:0 anteiso'),
                         gram_neg_lipids <- c('14:0 3OH', '16:1 w7c', '16:1 w9c', '18:1 w7c',
                                              '18:1 w9t'),
                         actino_lipids <- c('16:0 10me', '18:0 10me'),
                         amf_lipids <- c('16:1 w5c'),
                         s_fungi_lipids <- c('18:1 w9c', '18:1 w6,9c'),
                         anaerobe_lipids <- c('19:0 cyclo')
  )
  # Naming list facilitates naming the output of calcs
  names(indicator_list) <- c('f_lipids', 'b_lipids', 'gram_pos_lipids', 'gram_neg_lipids',
                             'actino_lipids', 'amf_lipids', 's_fungi_lipids',
                             'anaerobe_lipids')


  nmol_df <- reshape(data = df[c('DataFileName', 'Name', 'nmol_g')],
                   timevar = 'Name', idvar = 'DataFileName',
                   direction = 'wide')
  names(nmol_df) <- gsub('nmol_g.', '', names(nmol_df))
  nmol_df[is.na(nmol_df)] <- 0

  # without drop = FALSE, subsetting returns vector and rowSums() errors
  nmol_df['total_biomass'] <- rowSums(nmol_df[, 2:ncol(nmol_df), drop = FALSE],
                                      na.rm = TRUE) # need to limit this to only microbial lipids

  # Convert biomarker concentrations to % of total biomass -- I think this is redundant
  #perc_df <- cbind(nmol_df['DataFileName'], nmol_df['total_biomass'],
  #                 lapply(nmol_df[, 2:(ncol(nmol_df)-1)],
  #                        function(x){x/nmol_df[['total_biomass']] * 100}))


  perc_df <- cbind(nmol_df['DataFileName'], nmol_df['total_biomass'],
                   lapply(indicator_list,
                          function(x){
                           tmp_df <- nmol_df[,
                                             names(nmol_df) %in% x,
                                             drop = FALSE]
                           rowSums(tmp_df, na.rm = TRUE) /
                             nmol_df[, 'total_biomass'] * 100
                          }
    )
  )

  perc_df['fb'] <- rowSums(nmol_df[,
                                   names(nmol_df) %in% indicator_list[['f_lipids']],
                                   drop = FALSE]) /
                   rowSums(nmol_df[,
                                   names(nmol_df) %in% indicator_list[['b_lipids']],
                                   drop = FALSE])

  # calculate total mass/abundance of lipids belonging to a microbial group
  nm_df <- cbind(nmol_df['DataFileName'],
                   lapply(indicator_list,
                          function(x){
                            tmp_df <- nmol_df[,
                                              names(nmol_df) %in% x,
                                              drop = FALSE]
                            rowSums(tmp_df, na.rm = TRUE)
                          }
                   )
  )

  perc_df_long <- reshape(perc_df, varying = names(perc_df)[c(2:10)],
                          v.names = 'Percent_or_fraction',
                          timevar = 'Indicator', idvar = 'DataFileName',
                          times = names(perc_df)[c(2:10)], direction = 'long')

  nm_df_long <- reshape(nm_df, varying = names(nm_df)[c(2:9)],
                          v.names = 'nmol_g',
                          timevar = 'Indicator', idvar = 'DataFileName',
                          times = names(perc_df)[c(2:9)], direction = 'long')

#  b <- as.data.frame(sapply(a,'['))

#  b <- lapply(indicator_list, function(x){
#    tmp_df <- nmol_df[, colnames(nmol_df)[colnames(nmol_df) %in% x]]
#    rowSums(as.data.frame(tmp_df), na.rm = TRUE)
#  }
#  )


#  b <- cbind(nmol_df[1], b)

  perc_df_long <- merge(perc_df_long, soil_wt_df, by = 'DataFileName', all.x = TRUE)

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
