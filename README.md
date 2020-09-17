plfaR Readme
================

This is a package for running QC on and processing named PLFA data. The
three major components are calculating lipid concentrations, biomarker
concentrations and mol %, and d13C isotope data, with output in a format
that can easily be further analyzed using tidyverse and statistical
packages in R.

# Getting started

## Installation

If the devtools package isn’t already installed, download and install
from CRAN.

``` r
install.packages("devtools")
```

Install plfaR package.

``` r
devtools::install_github(repo = 'mlfelice/plfaR')
```

The package is now ready for processing PLFA data.

## Basic Workflow

### Summary

This vignette shows an example workflow for basic importing and
processing of PLFA data from raw named peak lists to lipid and biomarker
concentrations and mol %.

### Processing

Import packages

``` r
library(plfaR)  # PLFA processing
library(tidyverse)  # Basic data wrangling
library(readxl)  # Import Excel files
```

Import sample metadata. At minimum, you should have the following
columns: \* Batch \* DataFileName \* SampleType \* SampleWt You can also
include any additional data you will want for analysis.

``` r
stm_1415 <- read_xls(path = paste0('C:\\Users\\Mark\\Dropbox\\',                                                                          'umn_gutknecht_postdoc\\',
                                   'spruce_project\\', 
                                   'plfa_13c\\data\\',
                                   'SPRUCE DATA Nov 2016_all the ',
                                   'final calculations.xls'), 
                     sheet = 'Stats Ready Data Values', na = 'NA',
                     range = cell_cols('A:I')) %>%
  rename(BatchDataFileName = `Sample ID`) %>%
  mutate(Month = factor(Month, levels = c('June', 'July', 'August', 
                                          'Sept', 'October'))) %>%
  select(c(BatchDataFileName, Plot, DepthInterval = `Depth Increment`, 
           CO2EventualTreatment = Treatment, Temp = `Soil Probe Temp`, 
           GWCPerc = `% Gravimetric Water Content`)) # Excluded Year and Month to avoid dups in merge

c_types1 <- c('text', 'text', 'text', 'text', 'numeric', 'numeric', 'text', 
              'numeric', 'text', 'text', 'numeric', 'numeric', 'numeric',
              'text', 'numeric')
md_1415 <- read_xls(path = paste0('C:\\Users\\Mark\\Dropbox\\',
                                  'umn_gutknecht_postdoc\\',
                                  'spruce_project\\', 
                                  'plfa_13c\\data\\',
                                  'SPRUCE DATA Nov 2016_all the ',
                                  'final calculations.xls'), 
                    sheet = 'Data for pivot tables', na = 'NA',
                    col_types = c_types1) %>%
  mutate(TopCoreDepthcm = as.numeric(if_else(str_extract(`Depth Increment (cm)`, 
                                                         '(?<=\\()[0-9]+') == '1', 
                                             '0', 
                                             paste0('-', str_extract(`Depth Increment (cm)`, 
                                                                     '(?<=\\()[0-9]+')
                                             )
  )
  ),
  BottomCoreDepthcm = as.numeric(paste0('-', 
                                        str_extract(`Depth Increment (cm)`, 
                                                    '(?<=, )[0-9]+'))),
  MidCoreDepthcm = TopCoreDepthcm + (BottomCoreDepthcm - TopCoreDepthcm)/2) %>%
  rename(BatchDataFileName = `ID_with Batch ID`) %>%
  distinct(BatchDataFileName, .keep_all = TRUE) %>% # cut metadata repeated for each biomarker
  left_join(stm_1415, by = c('BatchDataFileName')) %>%
  select(c(DataFileName = Name, Batch = `Batch ID`, SampleType = Category, 
           Month, Year, Plot, DepthInterval, MidCoreDepthcm, 
           SampleWt = `Weight (g)`, Temp, GWCPerc))
```

    ## Warning in read_fun(path = enc2native(normalizePath(path)), sheet_i = sheet, :
    ## Expecting numeric in M2261 / R2261C13: got '.'

    ## Warning in mask$eval_all_mutate(dots[[i]]): NAs introduced by coercion

``` r
# Cell M2261 is blank, gives warning about expecting numeric

# Clean up workspace
rm(c_types1, stm_1415)
```

Import your PLFA data. This should be your named peak list, and should
be an Excel file with the following formatting: \* Either single tab OR
peak list data on a worksheet named ‘named\_peaks’ \* File should have
the following columns with headers: + Batch + DataFileName + RetTimeSecs
+ MajorHeightnA + TotalPeakArea1 + DisplayDelta1 + Name

``` r
# Define file path
source_dir <- 'C:\\Users\\Mark\\Desktop\\temp plfa practice - can delete whenever'  
source_file <- 'spruce_2015_b14_example.xlsx'

# Imported data is dataframe with cols listed above
peak_list <- import_batch(paste0(source_dir, '/', source_file))
```

Run a QC tests on the data. If there are any issues with the data, an
error message will be displayed indicating the issue. This also
generates a list containing summary dataframes showing 1) samples with
duplicate lipids (and which lipids), 2) samples missing standard peaks
(13:0, 16:0\*, or 19:0), and 3) the frequency (fraction) with which each
lipid is present in a sample.

\*NOTE: Although 16:0 isn’t technically a standard, it is present in
Sphagnum, so should be present in any peat samples. This can likely be
ignored if working with mineral soils.

``` r
qc_df <- quality_check(peak_list)
```

    ## Batch: B14
    ## ------

    ## Warning: Standards missing from at least one sample.
    ## Check data file and/or chromatograms before proceeding

``` r
# QC summary dataframes can be accessed individually
qc_df[[1]][['duplicate_lipids']]
qc_df[[1]][['missing_lipids']]
qc_df[[1]][['lipid_frequency']]
```

To process the data, we first convert peak areas to nmol/g soil \* In
addition to supplying the peak list, we need to supply a vector
containing the sample names for the blanks (should only have 13:0 and
19:0) and standards (ie. the samples on which you want to base kval
calculations.

NOTE: process\_peak\_area() also has several other parameters that can
be changed, but which are constant in our lab.

``` r
# Create vector of filenames of standards (can also manually supply)
stds <- c('13 Standard 1.raw', '4x concentraed 19 standard.raw')

# Blanks represent peak area subtracted for 13:0 and 19:0, so input a named 
# list with a vector for each lipid to subtract
blanks <- list( '13:0' = stds,
                '19:0' = c('Control 1 for plot 16 17 and 7.raw',
                           'Control 1 for plots 19 20  21.raw',
                           'Control 2 for plot 16 17 and 7.raw',
                           'Control 2 for plots 19 20  21.raw)'))

# Convert the peak areas to concentration (nmol/g dry soil)
conc_df <- process_peak_area(peak_list, 
                                  blanks = blanks, 
                                  standard_fnames = stds,
                                  soil_wt_df = md_1415)
```

Next, we can calculate the concentrations of groups of indicator lipids
as well as their mol %.

``` r
# I'm removing 16:0 and 19:0 peaks for this analysis
# use the output from concentration calcs as the input
indicators_df <- filter(conc_df, !Name %in% c('16:0', '19:0')) %>%
                 calculate_indicators(soil_wt_df = md_1415)
```

These dataframes can now be analyzed with your preferred stats pipeline.

### Other information

The following are the groupings of lipids used for each biomarker group.

    ## $f_lipids
    ## [1] "18:1 w9c"   "18:2 w6,9c"
    ## 
    ## $b_lipids
    ## [1] "15:0 iso"     "15:0 anteiso" "16:1 w7c"     "16:0 10me"    "18:1 w9c"    
    ## [6] "18:1 w9t"     "18:2 w6,9c"   "18:0 10me"    "19:0 cyclo"  
    ## 
    ## $gram_pos_lipids
    ## [1] "15:0 iso"     "15:0 anteiso"
    ## 
    ## $gram_neg_lipids
    ## [1] "16:1 w7c" "18:1 w9t"
    ## 
    ## $actino_lipids
    ## [1] "16:0 10me" "18:0 10me"
    ## 
    ## $anaerobe_lipids
    ## [1] "19:0 cyclo"
    ## 
    ## $protozoa
    ## [1] "20:4 w6,9,12,15"
    ## 
    ## $total_biomass
    ##  [1] "8:0"          "10:0"         "10:0 2OH"     "11:0"         "12:0"        
    ##  [6] "12:0 2OH"     "12:0 3OH"     "13:0"         "14:0"         "14:0 2OH"    
    ## [11] "14:0 3OH"     "14:1"         "15:0"         "15:0 anteiso" "15:0 iso"    
    ## [16] "16:0 10me"    "16:0 2OH"     "16:0 iso"     "16:1 w5c"     "16:1 w7c"    
    ## [21] "16:1 w9c"     "17:0"         "17:0 iso"     "17:0 anteiso" "17:0 cyclo"  
    ## [26] "17:1"         "17:1 iso"     "18:0"         "18:0 10me"    "18:1 w9c"    
    ## [31] "18:1 w9t"     "18:2 w6,9c"   "19:0"         "19:0 cyclo"   "19:1"        
    ## [36] "20:0"
