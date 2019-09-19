
subt_stand <- function(df, stand_df){

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
                    standard = '13:0', soil_mass = 1, vial_vol = 20,
                    start_col = 2){
  # Converts the peak area to lipid concentration
  #
  # Args:
  #   df: dataframe of normalized peak areas with 13:0 and 19:0 subtracted
  #   stand_df: dataframe identifying filenames for standards
  #   mwt_df: reference dataframe containing molecular weights of biomarkers
  #   standard_conc: concentration (nmol/uL) of specified standard used
  #   inj_vol: volume of specified standard injected
  #   standard: lipid standard used for area to conc calculation
  #   vial_vol: total volume in GC vial
  #   start_col: index of column with first biomarker
  #
  # Returns:
  #   dataframe containing total biomass
  temp_df <- df[, 3:ncol(df)]  # Remove sample names from dataframe

  fnames_df <- subset(stand_df, biomarker == standard)
  stand_vec <- df[which(df$DataFileName %in% fnames_df[['name']]), standard]
  stand_val <- mean(stand_vec)
  kval <- stand_val / standard_conc / inj_vol

  for(i in mwt_df[[1]]){  # First column of reference dataframe (biomarker names)
    mw <- get_mw(mwt_df, i)
    df[i] <- (df[[i]] / kval) * (vial_vol / 2) / (mw * soil_mass)
  }


  cat('No reference molecular weight for: \n\n')
  print(
    names(df[,4:ncol(df)])[which(!names(df[,4:ncol(df)]) %in% mwt_df[[1]])]
  )

  df <- df %>% mutate(
    #total_biomass = select_(.dots = match(mwt_df[[1]], names(df))) %>% rowSums) %>%
    total_biomass = rowSums(select(df, which(names(.) %in% mwt_df[[1]])),
                            na.rm = TRUE)) #%>%
  #select(batch, DataFileName, total_biomass)

  return(df)

}
