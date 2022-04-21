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
# Import metadata Excel file
md <- read_xlsx('plfaR_example-metadata.xlsx')

# Preview metadata
head(md)
```

    ## # A tibble: 6 x 4
    ##   Batch DataFileName                       SampleType SampleWt
    ##   <chr> <chr>                              <chr>         <dbl>
    ## 1 B14   13 Standard 1.raw                  Standard         NA
    ## 2 B14   Blank 1 for plot 16 17 and 7.raw   Blank            NA
    ## 3 B14   Blank 1 for plot 19 20  21.raw     Blank            NA
    ## 4 B14   Blank 2 for plot 16 17 and 7.raw   Blank            NA
    ## 5 B14   Blank 2 for plot 19 20  21.raw     Blank            NA
    ## 6 B14   Control 1 for plot 16 17 and 7.raw Control          NA

Import your PLFA data. This should be your named peak list, and should
be an Excel file with the following formatting: \* Either single tab OR
peak list data on a worksheet named ‘named\_peaks’ \* File should have
the following columns with headers: + Batch + DataFileName + RetTimeSecs
+ MajorHeightnA + TotalPeakArea1 + DisplayDelta1 + Name

``` r
peak_list <- import_batch('plfaR_example-peak-list.xlsx')

# Preview imported data
head(peak_list)
```

    ## # A tibble: 6 x 8
    ##   BatchDataFileNa~ Batch DataFileName RetTimeSecs MajorHeightnA TotalPeakArea1
    ##   <chr>            <chr> <chr>              <dbl>         <dbl>          <dbl>
    ## 1 B14_Plot 16 100~ B14   Plot 16 100~        526.          1.88        443498.
    ## 2 B14_Plot 17 0-2~ B14   Plot 17 0-2~        526.          1.96        458006.
    ## 3 B14_Plot 17 20-~ B14   Plot 17 20-~        526.          1.32        311325.
    ## 4 B14_Plot 17 50-~ B14   Plot 17 50-~        526.          2.60        695964.
    ## 5 B14_Plot 17 100~ B14   Plot 17 100~        526.          1.82        490706.
    ## 6 B14_Plot 17 150~ B14   Plot 17 150~        526.          2.57        718811.
    ## # ... with 2 more variables: DisplayDelta1 <dbl>, Name <chr>

Run QC tests on the data. If there are any issues with the data, an
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

To process the data, we first calculate lipid concentrations in nmol/g
soil \* In addition to supplying the peak list, we need to supply a
vector containing the sample names for the blanks (should only have 13:0
and 19:0) and standards (ie. the samples on which you want to base kval
calculations).

NOTE: process\_peak\_area() also has several other parameters that can
be changed, but which are constant in our lab.

``` r
# Create vector of filenames of standards (can also manually supply)
stds <- c('13 Standard 1.raw')

# Blanks represent peak area subtracted for 13:0 and 19:0, so input a named 
# list with a vector for each lipid to subtract
blanks <- list( '13:0' = c('13 Standard 1.raw'),
                '19:0' = c('Control 1 for plot 16 17 and 7.raw',
                           'Control 1 for plots 19 20  21.raw',
                           'Control 2 for plot 16 17 and 7.raw',
                           'Control 2 for plots 19 20  21.raw)'))


# Convert the peak areas to concentration (nmol/g dry soil)
conc_df <- process_peak_area(peak_list, 
                             blanks = blanks, 
                             standard_fnames = stds,
                             soil_wt_df = md)

head(conc_df)
```

    ##       Name Batch             DataFileName            BatchDataFileName
    ## 1 10:0 2OH   B14 Plot 17 150-200cm 8g.raw B14_Plot 17 150-200cm 8g.raw
    ## 2 10:0 2OH   B14  Plot 17 50-100cm 8g.raw  B14_Plot 17 50-100cm 8g.raw
    ## 3 10:0 2OH   B14 Plot 17 100-150cm 8g.raw B14_Plot 17 100-150cm 8g.raw
    ## 4 10:0 2OH   B14 Plot 21 50-100cm  8g.raw B14_Plot 21 50-100cm  8g.raw
    ## 5 10:0 2OH   B14 Plot 20 20-50 cm  6g.raw B14_Plot 20 20-50 cm  6g.raw
    ## 6 10:0 2OH   B14 Plot 19 50-100cm  8g.raw B14_Plot 19 50-100cm  8g.raw
    ##   RetTimeSecs MajorHeightnA TotalPeakArea1 DisplayDelta1 BlnkArea
    ## 1       714.1     1.9397842       746308.6    -17.226493        0
    ## 2       714.2     1.6354058       821001.4     -8.610726        0
    ## 3       714.1     1.6381963       646790.2    -34.335439        0
    ## 4       714.4     2.2091562      1003103.8    -26.853440        0
    ## 5       714.9     0.5327858       194985.3    -30.641190        0
    ## 6       714.6     1.9758513       971162.6    -26.955573        0
    ##   AreaMinusBlanks SampleType SampleWt molecular_weight_g_per_mol indicates
    ## 1        746308.6    Unknown        8                        202      <NA>
    ## 2        821001.4    Unknown        8                        202      <NA>
    ## 3        646790.2    Unknown        8                        202      <NA>
    ## 4       1003103.8    Unknown        8                        202      <NA>
    ## 5        194985.3    Unknown        6                        202      <NA>
    ## 6        971162.6    Unknown        8                        202      <NA>
    ##   C_PLFA    nmol_g
    ## 1     10 0.9448399
    ## 2     10 1.0394022
    ## 3     10 0.8188478
    ## 4     10 1.2699471
    ## 5     10 0.3291398
    ## 6     10 1.2295089

Next, we can calculate the concentrations of groups of indicator lipids
as well as their mol %.

``` r
indicators_df <- calculate_indicators(conc_df, soil_wt_df = md)

head(indicators_df)
```

    ##        DataFileName       Indicator Percent Batch SampleType SampleWt nmol_g
    ## 1 13 Standard 1.raw   actino_lipids     NaN   B14   Standard       NA      0
    ## 2 13 Standard 1.raw anaerobe_lipids     NaN   B14   Standard       NA      0
    ## 3 13 Standard 1.raw        b_lipids     NaN   B14   Standard       NA      0
    ## 4 13 Standard 1.raw        f_lipids     NaN   B14   Standard       NA      0
    ## 5 13 Standard 1.raw          f_to_b     NaN   B14   Standard       NA     NA
    ## 6 13 Standard 1.raw gram_neg_lipids     NaN   B14   Standard       NA      0

Next, we can calculate the corrected isotope values for individual
lipds. This will calculate 2 corrections: 1) normalize measured isotope
values to the international reference scale, and 2) correct values for
the d13C of methanol added during methylation. The lipid d13C values can
be used to calculate the d13C of indicator groups.

Indicator group d13C are calculated from the unweighted average of
individual lipid d13C.

``` r
# Process lipid d13C to normalize to internationl refernce scale and account 
# for the contribution of methanol to d13C
d13c_lipid_df <- correct_iso_base(df = peak_list, 
                                  d13c_correction = 0, 
                                  methanol_13c = -55.84, 
                                  min_height = 1)
```

    ## 144 of 364 peaks removed from batch B14

``` r
# Calculate d13C of indicator groups
d13c_indicator_df <- indic_iso_base(d13c_lipid_df, md)

head(d13c_lipid_df)
```

    ##       Name             BatchDataFileName Batch              DataFileName
    ## 1 10:0 2OH  B14_Plot 17 150-200cm 8g.raw   B14  Plot 17 150-200cm 8g.raw
    ## 2 10:0 2OH  B14_Plot 17 100-150cm 8g.raw   B14  Plot 17 100-150cm 8g.raw
    ## 5 10:0 2OH  B14_Plot 21 50-100cm  8g.raw   B14  Plot 21 50-100cm  8g.raw
    ## 6 10:0 2OH   B14_Plot 17 50-100cm 8g.raw   B14   Plot 17 50-100cm 8g.raw
    ## 7 10:0 2OH  B14_Plot 16 150-200cm 8g.raw   B14  Plot 16 150-200cm 8g.raw
    ## 8 10:0 2OH B14_Plot 20 150-200cm  7g.raw   B14 Plot 20 150-200cm  7g.raw
    ##   RetTimeSecs MajorHeightnA TotalPeakArea1 DisplayDelta1
    ## 1       714.1      1.939784       746308.6    -17.226493
    ## 2       714.1      1.638196       646790.2    -34.335439
    ## 5       714.4      2.209156      1003103.8    -26.853440
    ## 6       714.2      1.635406       821001.4     -8.610726
    ## 7       714.4      1.076172       460225.0    -17.875433
    ## 8       714.4      1.575394       618243.9    -24.667249
    ##   molecular_weight_g_per_mol indicates C_PLFA d13C_corrected
    ## 1                        202      <NA>     10     -13.365143
    ## 2                        202      <NA>     10     -32.184982
    ## 5                        202      <NA>     10     -23.954784
    ## 6                        202      <NA>     10      -3.887799
    ## 7                        202      <NA>     10     -14.078976
    ## 8                        202      <NA>     10     -21.549974

``` r
head(d13c_indicator_df)
```

    ##        DataFileName     BatchDataFileName       Indicator avg_d13C_corrected
    ## 1 13 Standard 1.raw B14_13 Standard 1.raw gram_pos_lipids                 NA
    ## 2 13 Standard 1.raw B14_13 Standard 1.raw        protozoa                 NA
    ## 3 13 Standard 1.raw B14_13 Standard 1.raw   total_biomass           -28.8483
    ## 4 13 Standard 1.raw B14_13 Standard 1.raw   actino_lipids                 NA
    ## 5 13 Standard 1.raw B14_13 Standard 1.raw        b_lipids                 NA
    ## 6 13 Standard 1.raw B14_13 Standard 1.raw anaerobe_lipids                 NA
    ##   Batch SampleType SampleWt
    ## 1   B14   Standard       NA
    ## 2   B14   Standard       NA
    ## 3   B14   Standard       NA
    ## 4   B14   Standard       NA
    ## 5   B14   Standard       NA
    ## 6   B14   Standard       NA

These dataframes can now be analyzed with your preferred stats pipeline.

### Other information

The following are the groups of lipids used for each biomarker group.

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
