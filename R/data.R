#' Molecular weights of microbially-derived phospholipids
#'
#' COmpiled from Jess Gutknecht and literature review. Weights are for
#' methylated lipids and are rounded to the nearest whole number.
#'
#' @format A dataframe with  rows and  variables:
#' \describe{
#'   \item{fame}{Lipid name, following __ convention}
#'   \item{molecular_weight_g_per_mol}{methylated lipid molecular weight}
#'   \item{indicates}{biomarker which lipid corresponds to}
#'   \item{C_PLFA}{number of carbon atoms in the lipid}
'lipid_reference'

# To update package data to reflect source edits, run code below:

#lipid_reference <- readxl::read_xlsx(paste0('C:/Users/Mark/Dropbox/',
#                                    'umn_gutknecht_postdoc/spruce_project/',
#                                    'plfa_13c/20190521_spruce_plfa_biomarker',
#                                    '-molecular-wt.xlsx')
#)
#
#devtools::use_data(lipid_reference, overwrite = TRUE)
