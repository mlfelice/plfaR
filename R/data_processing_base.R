subtract_blanks <- function(df, blanks = c('1.raw', '2.raw'),
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

area_to_concentration <- function(df, standard_fnames, mw_df = lipid_reference,
                               standard_conc = 250, inj_vol = 2, #standard_fnames should be blanks, not 13:0 standard, I think
                               standard = '13:0', soil_wt_df, vial_vol = 20){



  standard_vec <- df[df[['DataFileName']] %in% standard_fnames &
                             df[['Name']] == standard, 'TotalPeakArea1']
  kval <- mean(standard_vec) / standard_conc / inj_vol

  tmp_df <- merge(df, samp_wt_df, by = c('Batch', 'DataFileName'), all.x = TRUE)
  tmp_df <- merge(tmp_df, mw_df, by.x = 'Name', by.y = 'fame')
  tmp_df <- transform(tmp_df, nmol_g = (TotalPeakArea1 / kval) *
                         (vial_vol / 2) /
                         (molecular_weight_g_per_mol * SampleWt))

  return(tmp_df)

}

process_peak_area <- function(dat, standard_fnames, mw_df = lipid_reference,
                              standard_conc = 250, inj_vol = 2, #standard_fnames should be blanks, not 13:0 standard, I think
                              standard = '13:0', soil_wt_df, vial_vol = 20,
                              blanks = c('1.raw', '2.raw'),
                              lipids = c('13:0', '19:0')){

  dtype <- check_dtype(dat)

  if (dtype == 'data.frame') {
    tmp_df <- subtract_blanks(dat, blanks = c('1.raw', '2.raw'),
                    lipids = c('13:0', '19:0'))
    area_to_concentration(tmp_df, standard_fnames, mw_df = lipid_reference,
                      standard_conc = 250, inj_vol = 2, #standard_fnames should be blanks, not 13:0 standard, I think
                      standard = '13:0', soil_wt_df, vial_vol = 20)
  }
}
