###############################################################################
# This file includes the current implementation of functions for processing raw
# PLFA data into concentrations, total biomass, and indicator group abundance/
# concentration.
###############################################################################

#' Subtract peak area from control and standard lipids
#'
#' All samples are typically run with a 19:0 extraction control and 13:0
#' standard. This function subtracts peak area contributed by these additions
#' prior to converting peak area to concentration.
#'
#' @param df Dataframe or tibble with named PLFA peak data following format
#' specified by \code{import_batch()}.
#'
#' @param blanks Vector of the names of samples run as blanks. This should
#' include all samples you want to use for peak subtraction (eg. 13:0 and 19:0
#' in most cases).
#'
#' @return Tibble similar to input with the addition of AreaMinusBlanks column
#' showing peak areas minus area contributed by standard/control lipids.
#'
#' @examples
#'
#' @export
#'
subtract_blanks <- function(df, blanks){  # version that subtracts avg of blanks only
  # Inputs:
  #   df (dataframe): peak list dataframe
  #   blanks (list): named list with list item names corresponding to lipid,
  #     and each item consisting of a character vector containing the
  #     DataFileName corresponding to the samples to use for blanks
  #     subtraction. Blanks are non-experimental samples containing surrogate
  #     or internal standards also present in experimental samples that need
  #     to be subtracted from experimental samples. In the Gutknecht lab, this
  #     is typically 13:0 and 19:0.
  #
  #
  # Would it make more sense to just add info on the role of sample in the
  # metadata?

  # TO DO: add ability to check list names against lipid names in the reference
  # so that it throws an error if you've chosen a non-existant lipid

  # Filter peak list dataframe for blank samples (defined by func arg)
  # The [[1]] is necessary, as subsetting list by name returns a list, even
  # with double brackets
  b_lipids <- names(blanks)
  blanks_df <- do.call('rbind',
                       lapply(b_lipids,
                              function(x){df[df[['DataFileName']] %in%
                                               blanks[[x]] &
                                             df[['Name']] == x, ]}))

  # Average the peak area of each lipid within each batch
  # Resulting df has 3 cols: Batch, Name, and TotalPeakArea1
  blanks_df <- aggregate(TotalPeakArea1 ~ Batch + Name,
                         FUN = function(x){
                           mean(x, na.rm = TRUE)  # modified
                         }, data = blanks_df
  )

  # Rename TotalPeakArea1 to StdArea so that resulting df shows value used for
  # blank/peak subtraction when we merge
  names(blanks_df)[names(blanks_df) == 'TotalPeakArea1'] <- 'BlnkArea'

  full_df <- merge(df, blanks_df, by = c('Batch', 'Name'), all.x = TRUE)
  #full_df[is.na(full_df[['StdArea']]), 'StdArea'] <- 0
  # After merging, lipids other than the blanks you are subtracting will be NA
  # This converts to 0, so that the arithmetic produces a number rather than
  # NA
  full_df <- transform(full_df,
                       BlnkArea = ifelse(is.na(BlnkArea), 0, BlnkArea))
  full_df <- transform(full_df,
                       AreaMinusBlanks = ifelse(TotalPeakArea1 - BlnkArea > 0,
                                                TotalPeakArea1 - BlnkArea, 0))

  return(full_df)
}

#' Convert peak areas to concentration (nmol/g dry soil)
#'
#' This function takes a formatted peak list and calculates kval in order to
#' convert peak areas to concentration in units of nmol of lipid/g dry soil.
#' This function does not include blank subtraction. If blanks need to be
#' subtracted, use \code{process_peak_area()}.
#'
#' @param df Dataframe or tibble with named PLFA peak data following format
#' specified by \code{import_batch()}.
#'
#' @param standard_fnames Character vector with names of samples to be used as
#' standards for calculating the kval. Names should match names in DataFileName
#' column of input dataframe.
#'
#' @param mw_df Dataframe including all lipids (fame column) of interest and
#' their associated molecular weights (molecular_weight_g_per_mol column). By
#' default, this function uses the included \code{lipid_reference} dataframe.
#'
#' @param standard_conc Numeric value indicating the concentraion of standard
#' in ng/uL. Default is 250 ng/uL.
#'
#' @param inj_vol Volume of sample injected into GC in uL. Default is 2uL.
#'
#' @param standard String indicating the lipid used as the standard. Default is
#' 13:0.
#'
#' @param soil_wt_df Dataframe or tibble containing the dry weight (g) of soil
#' used for lipid extraction for each sample.
#'
#' @param vial_vol Volume of solution run on GC. Default is 50 uL
#'
#' @param blanks Vector of the names of samples run as blanks. This should
#' include all samples you want to use for peak subtraction (eg. 13:0 and 19:0)
#' in most cases.
#'
#' @return Tibble similar to input with the addition of concentration column.
#'
#' @examples
#'
#' @export
#'
area_to_concentration <- function(df, avg_std_area, soil_wt_df,
                                        mw_df = lipid_reference,
                                        standard_conc, inj_vol, standard, #standard_fnames should be blanks, not 13:0 standard, I think
                                        vial_vol){
  # standard fnames seems kind of redundant, since we're just pulling the
  # Inputs:
  #   df (dataframe):
  #   standard_fnames (???): replaced by avg_std_area. This was the vector of
  #     DataFileName s associated with standards we wanted to use
  #   soil_wt_df (numeric):
  #   mw_df: (dataframe):
  #   standard_conc (numeric):
  #   inj_vol (numeric):
  #   standard (string):
  #   vial_vol (numeric):
  #   standard_vec (numeric vector): This is a vector of peak areas for
  #     the standard biomarker (usu 13:0) used for kval calc. These are
  #     associated with the DataFileNames defined by standard_fnames.
  #     I think it makes more sense to provide only the standard vector or avg
  #     peak area as an arg. #### working on changing this now ###
  #     ### This is now replaced by avg_std_area
  #

  # moving this to process_peak_area2()
  # standard_vec <- unlist(df[df[['DataFileName']] %in% standard_fnames &
  #                            df[['Name']] == standard, 'TotalPeakArea1'])

  kval <- avg_std_area / standard_conc / inj_vol

  # Merge the metadata and lipid reference dataframs to the input peak list
  # These have parameters needed for calculations
  # If/else statement ensures Batch is same type (char/num) between df,
  # otherwise merge would fail silently
  if(class(df[['Batch']]) == class(soil_wt_df[['Batch']])){
    tmp_df <- merge(df, soil_wt_df, by = c('Batch', 'DataFileName'),
                    all.x = TRUE)
  } else{stop('Batch column type does not match between df and soil_wt_df')}

  tmp_df <- merge(tmp_df, mw_df, by.x = 'Name', by.y = 'fame')

  # Run actual calculation applying kval to all peak areas
  tmp_df <- transform(tmp_df, nmol_g = (AreaMinusBlanks / kval) *
                        (vial_vol / 2) /
                        (molecular_weight_g_per_mol * SampleWt))

  return(tmp_df)

}

#' Convert peak areas to concentration (nmol/g dry soil)
#'
#' This function takes a formatted peak list and carries out blank subtraction
#' and calculates kval in order to convert peak area for each peak to
#' concentration.
#'
#' @param dat Dataframe or tibble with named PLFA peak data following format
#' specified by \code{import_batch()}.
#'
#' @param standard_fnames Character vector with names of samples to be used as
#' standards for calculating the kval. Names should match names in DataFileName
#' column of input dataframe.
#'
#' @param mw_df Dataframe including all lipids (fame column) of interest and
#' their associated molecular weights (molecular_weight_g_per_mol column). By
#' default, this function uses the included \code{lipid_reference} dataframe.
#'
#' @param standard_conc Numeric value indicating the concentraion of standard
#' in ng/uL. Default is 250 ng/uL.
#'
#' @param inj_vol Volume of sample injected into GC in uL. Default is 2uL.
#'
#' @param standard String indicating the lipid used as the standard. Default is
#' 13:0.
#'
#' @param soil_wt_df Dataframe or tibble containing the dry weight (g) of soil
#' used for lipid extraction for each sample.
#'
#' @param vial_vol Volume solution run on GC. Default is 50 uL
#'
#' @param blanks Vector of the names of samples run as blanks. This should
#' include all samples you want to use for peak subtraction (eg. 13:0 and 19:0)
#' in most cases.
#'
#' @return Tibble similar to input with the addition of concentration column.
#'
#' @examples
#'
#' @export
#'
process_peak_area <- function(dat, standard_fnames, mw_df = lipid_reference,
                                   standard_conc = 250, inj_vol = 2, #standard_fnames should be blanks, not 13:0 standard, I think
                                   standard = '13:0', soil_wt_df, vial_vol = 50,
                                   blanks){

  # Inputs:
  #   dat (dataframe or list?):
  #   standard_fnames (char vector):
  #   mw_df (dataframe):
  #   standard_conc (numeric):
  #   inj_vol (numeric):
  #   standard (string):
  #   vial_vol (numeric):
  #   blanks (list): named list with list item names corresponding to lipid,
  #     and each item consisting of a character vector containing the
  #     DataFileName corresponding to the samples to use for blanks
  #     subtraction. Blanks are non-experimental samples containing surrogate
  #     or internal standards also present in experimental samples that need
  #     to be subtracted from experimental samples. In the Gutknecht lab, this
  #     is typically 13:0 and 19:0.


  # This creates a numeric vector of the standard (usu 13:0) values
  # associated with supplied DataFileName
  standard_vec <- unlist(dat[dat[['DataFileName']] %in% standard_fnames &
                               dat[['Name']] == standard, 'TotalPeakArea1'])
  avg_std_area <- mean(standard_vec)
  # This is what gets used for applying
  #kval, values from subtract_blanks() only used for subtracting

  ####
  # could also structure so that user can input 'blanks' arg that would
  # supersede the SampleType from metadata
  # This is the old way. Probably best to do something where the metadata
  # tells which to use
  #blanks <- unlist(soil_wt_df[grepl('[Bb]lank', soil_wt_df[['SampleType']]),
  #                            'DataFileName'])
  ###
  dtype <- check_format(dat)

  if (dtype == 'data.frame') {
    tmp_df <- subtract_blanks(dat, blanks)
    #tmp_df <- subtract_blanks2(dat, blanks = blanks,
    #                               lipids = c('13:0', '19:0'))
    area_to_concentration(tmp_df, avg_std_area, mw_df = lipid_reference,
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

#' Calculate concentration and relative abundance of indicator lipids
#'
#' This function aggregates lipid concentration data into indicator groups
#' (eg. fungal lipids, bacterial lipids) and calculates their relative
#' abundance. Indicator groups reported as the sum of lipid concentration /
#' abundance of all lipids in that group. This also calculates total biomass
#' and fungal:bacterial lipid ratio.
#'
#' @param dat Dataframe or tibble with named PLFA peak data following format
#' specified by \code{import_batch()}.
#'
#' @param df Dataframe or tibble generated from \code{process_peak_area()} or
#' a dataframe/tibble with the same fields and format
#'
#' @param soil_wt_df Dataframe or tibble containing the dry weight (g) of soil
#' used for lipid extraction for each sample.
#'
#' @return Tibble similar to input with the addition of columns for indicator
#' groups (including total biomass), relative abundance, concentration, and
#' fungal:bacterial lipid ratio.
#'
#' @examples
#'
#' @export
#'

calculate_indicators <- function(df, soil_wt_df){

  indicator_list <- list(f_lipids = c('18:1 w9c', '18:2 w6,9c'),
                         # Lipids representing total bacteria for f to b ratio not shown in Cameron's data, so for now, I will use lipids for Gram +/-, actino, and anaerobe
                         b_lipids = c('15:0 iso', '15:0 anteiso',
                                      '16:1 w7c', '16:0 10me',
                                      '18:1 w9c', '18:1 w9t',
                                      '18:2 w6,9c', '18:0 10me', '19:0 cyclo'),
                         gram_pos_lipids = c('15:0 iso', '15:0 anteiso'),
                         gram_neg_lipids = c('16:1 w7c', '18:1 w9t'),
                         actino_lipids = c('16:0 10me', '18:0 10me'),
                         anaerobe_lipids = c('19:0 cyclo'),
                         protozoa = c('20:4 w6,9,12,15'),
                         # May need to redefine the total_biomass. Right now, it's just all of the
                         # indicator lipids. Not sure if there should be others as well
                         total_biomass = c('8:0', '10:0', '10:0 2OH', '11:0',
                                           '12:0', '12:0 2OH', '12:0 3OH',
                                           '13:0', '14:0', '14:0 2OH',
                                           '14:0 3OH', '14:1', '15:0',
                                           '15:0 anteiso', '15:0 iso',
                                           '16:0 10me', '16:0 2OH', '16:0 iso',
                                           '16:1 w5c', '16:1 w7c', '16:1 w9c',
                                           '17:0', '17:0 iso', '17:0 anteiso',
                                           '17:0 cyclo', '17:1', '17:1 iso',
                                           '18:0', '18:0 10me', '18:1 w9c',
                                           '18:1 w9t', '18:2 w6,9c', '19:0',
                                           '19:0 cyclo', '19:1', '20:0')
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
  if(class(df[['Batch']]) == class(soil_wt_df[['Batch']])){
    perc_df_long <- merge(perc_df_long, soil_wt_df,
                          by = 'DataFileName', all.x = TRUE)
  } else{stop('Batch column type does not match between df and soil_wt_df')}

  perc_df_long <- merge(perc_df_long, nm_df_long,
                        by = c('DataFileName', 'Indicator'), all.x = TRUE)

}


#' Convert peak areas to concentration (nmol/g dry soil)
#'
#' ## This is the most up to date peak area processing function (5/21/2022) ##
#' This function takes a formatted peak list and carries out blank subtraction
#' and calculates kval in order to convert peak area for each peak to
#' concentration. This new version is modified from \code{process_peak_area()}
#' in order to simplify specifying the standards and blanks. Instead of
#' inputting a list of lipids and blanks or lipids and standards, you simply
#' include a SampleType column in your weight sheet/sample denoting your blanks
#' and standards. Default arguments specify the most common lipids used as
#' standards and for blank subtraction and common SampleType names for them.
#' However, these can be explicitly entered to fit any analysis/naming scheme.
#'
#' @param dat Dataframe or tibble with named PLFA peak data following format
#' specified by \code{import_batch()}.
#'
#' @param standard_fnames Character vector with names of samples to be used as
#' standards for calculating the kval. Names should match names in DataFileName
#' column of input dataframe.
#'
#' @param mw_df Dataframe including all lipids (fame column) of interest and
#' their associated molecular weights (molecular_weight_g_per_mol column). By
#' default, this function uses the included \code{lipid_reference} dataframe.
#'
#' @param standard_conc Numeric value indicating the concentraion of standard
#' in ng/uL. Default is 250 ng/uL.
#'
#' @param inj_vol Volume of sample injected into GC in uL. Default is 2uL.
#'
#' @param standard String indicating the lipid used as the standard. Default is
#' 13:0.
#'
#' @param soil_wt_df Dataframe or tibble containing the dry weight (g) of soil
#' used for lipid extraction for each sample.
#'
#' @param vial_vol Volume solution run on GC. Default is 50 uL
#'
#' @param blanks named vector specifying names of samples to use for blanks
#' subtraction and which lipids need subtraction. The name of each vector item
#' should denote the lipid and each actual items in the vector should denote
#' the associated SampleType. Include all lipids you want to use for peak
#' subtraction. As an examnple, the default argument is
#' \code{blanks = c('19:0' = Blank, '13:0' = Standard)}.
#'
#' @return Tibble similar to input with the addition of concentration column.
#'
#' @examples
#'
#' @export
#'
process_peak_area_ <- function(dat,
                               #standard_fnames,
                               mw_df = lipid_reference,
                               standard_conc = 250, inj_vol = 2, #standard_fnames should be blanks, not 13:0 standard, I think
                               standard = '13:0', soil_wt_df, vial_vol = 50,
                               blanks = c('19:0' = 'Control', '13:0' = 'Standard')){
                                 # Inputs:
                                 #   dat (dataframe or list?):
                                 #   standard_fnames (char vector):
                                 #   mw_df (dataframe):
                                 #   standard_conc (numeric):
                                 #   inj_vol (numeric):
                                 #   standard (string):
                                 #   vial_vol (numeric):
                                 #   blanks (list): named list with list item names corresponding to lipid,
                                 #     and each item consisting of a character vector containing the
                                 #     DataFileName corresponding to the samples to use for blanks
                                 #     subtraction. Blanks are non-experimental samples containing surrogate
                                 #     or internal standards also present in experimental samples that need
                                 #     to be subtracted from experimental samples. In the Gutknecht lab, this
                                 #     is typically 13:0 and 19:0.


                                 # This creates a numeric vector of the standard (usu 13:0) values
                                 # associated with supplied DataFileName
                                 #####

                                # Pull the names of files used as standards and add to vector
                                 standard_fnames <- unlist(soil_wt_df[soil_wt_df[['SampleType']] == standard,
                                                                      'DataFileName'])

                                 # Pull the TotalPeakArea1 for the specified lipid associated with each file
                                 # There should only be one lipid used as a standard. Otherwise, will throw an
                                 # error.
                                 # TODO: update to allow multiple standards, or build in more informative
                                 # error handling.
                                 standard_vec <- unlist(dat[dat[['DataFileName']] %in% standard_fnames &
                                                              dat[['Name']] == names(standard), 'TotalPeakArea1'])


                                 #####
                                 #standard_vec <- unlist(dat[dat[['DataFileName']] %in% standard_fnames &
                                 #                             dat[['Name']] == standard, 'TotalPeakArea1'])
                                 avg_std_area <- mean(standard_vec)

                                 # This is what gets used for applying
                                 #kval, values from subtract_blanks() only used for subtracting

                                 ####
                                 # could also structure so that user can input 'blanks' arg that would
                                 # supersede the SampleType from metadata
                                 # This is the old way. Probably best to do something where the metadata
                                 # tells which to use
                                 #blanks <- unlist(soil_wt_df[grepl('[Bb]lank', soil_wt_df[['SampleType']]),
                                 #                            'DataFileName'])
                                 ###

                                 #####
                                 # This is the main differing code chunk. Creates a list of lipids for blank subtraction and associated file names
                                 # This creates a named list in same format as input for process_peak_area(), so we don't have to change
                                 # the subtract_blanks() function
                                 # Loops over the named vector and finds that data file names that matches each SampleType.
                                 # The output list items take on the name of the input vector item name
                                 blank_ls <- lapply(c('19:0' = 'Control', '13:0' = 'Standard'),
                                                    function(x){unlist(md[which(md[['SampleType']] ==
                                                                                  x), 'DataFileName'])})
                                 #####
                                 dtype <- check_format(dat)

                                 if (dtype == 'data.frame') {
                                   tmp_df <- subtract_blanks(dat, blank_ls)
                                   #tmp_df <- subtract_blanks2(dat, blanks = blanks,
                                   #                               lipids = c('13:0', '19:0'))
                                   area_to_concentration(tmp_df, avg_std_area, mw_df = lipid_reference,
                                                         standard_conc = standard_conc, inj_vol = inj_vol, #standard_fnames should be blanks, not 13:0 standard, I think
                                                         standard = standard, soil_wt_df, vial_vol = vial_vol)
                                 }
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

############################################
# This is the more comprehensive set of indicators I used for original analysis
# Eventually, we may want the option of using different pre-defined indicator
# sets

#indicator_list <- list(f_lipids = c('16:1 w5c', '18:1 w9c', '18:2 w6,9c'),
#                       b_lipids = c('13:0 iso', '13:0 anteiso', '14:0 3OH',
#                                    '15:0 iso', '15:0 anteiso', '16:0 iso',
#                                    '16:1 w7c', '16:0 10me', '17:0 iso',
#                                    '17:0 anteiso', '18:1 w9t', '18:1 w7c',
#                                    '18:0 10me'),
#                       gram_pos_lipids = c('13:0 iso', '13:0 anteiso',
#                                           '15:0 iso', '15:0 anteiso',
#                                           '17:iso', '17:0 anteiso'),
#                       gram_neg_lipids = c('14:0 3OH', '16:1 w7c',
#                                           '16:1 w9c', '18:1 w7c',
#                                           '18:1 w9t'),
#                       actino_lipids = c('16:0 10me', '18:0 10me'),
#                       amf_lipids = c('16:1 w5c'),
#                       s_fungi_lipids = c('18:1 w9c', '18:1 w6,9c'),
#                       anaerobe_lipids = c('19:0 cyclo'),
#                       # May need to redefine the total_biomass. Right now, it's just all of the
#                       # indicator lipids. Not sure if there should be others as well
#                       total_biomass = c("16:1 w5c", "18:1 w9c", "18:2 w6,9c",
#                                         "13:0 iso", "13:0 anteiso",
#                                         "14:0 3OH", "15:0 iso",
#                                         "15:0 anteiso", "16:0 iso",
#                                         "16:1 w7c", "16:0 10me", "17:0 iso",
#                                         "17:0 anteiso", "18:1 w9t",
#                                         "18:1 w7c", "18:0 10me", "17:iso",
#                                         "16:1 w9c", "18:1 w6,9c",
#                                         "19:0 cyclo")
#)



