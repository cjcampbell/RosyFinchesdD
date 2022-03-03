# RosyFinchesdD

Code and output of an analysis investigating the feather and claw origins of two species of rosy-finch using stable hydrogen isotope analysis.

## This repository contains:
### /out/
Species-specific raster files with probability-of-origin maps.
| Files                                   | Description |
| -----------                            | ----------- |
| BLRF_normalisedProbabilityMaps.grd,<br>BLRF_normalisedProbabilityMaps.gri  | Probability-of-origin maps from BLRF samples, normalized so probability of origin in focal region sums to 1.
| GCRF_normalisedProbabilityMaps.grd,<br>GCRF_normalisedProbabilityMaps.gri  | Probability-of-origin maps from GCRF samples, normalized so probability of origin in focal region sums to 1.
| BLRF_oddsProbabilityMaps.grd,<br>BLRF_oddsProbabilityMaps.gri  | Probability-of-origin maps from BLRF samples, transformed to relativized odds-of-origin values.
| GCRF_oddsProbabilityMaps.grd,<br>GCRF_oddsProbabilityMaps.gri  | Probability-of-origin maps from GCRF samples, transformed to relativized odds-of-origin values.

### /R/
Scripts used to read, tidy, and standardize data and run analyses.
| Script                                 | Description |
| -----------                            | ----------- |
| 00_Setup.R                             | Set up workspace, define extents, CRS |
| 01_loadIsotopeData.R                   | Import data, combine with metadata, tidy |
| 02_LoadSpatialData.R                   | Load administrative area, create  isoscapes, and create range map datasets |
| 03_FitApplyTransferFunctions.R         | Partially reproducing the results of Hobson 2012 to generate a bird keratin transfer function for migratory birds. Apply to create keratin isoscapes |
| 04_MakeProbabilityOfOriginMaps.R       | Create probability-of-origin maps for each sample, generate plots for each sample origin |

## Session Info
All analyses were conducted with the following configuration:

R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] tidyr_1.2.0        stringr_1.4.0      smatr_3.4-8        rmapshaper_0.4.5   rgdal_1.5-28      
 [6] readxl_1.3.1       purrr_0.3.4        measurements_1.4.0 lwgeom_0.2-8       lubridate_1.8.0   
[11] ggstar_1.0.3       ggmap_3.0.0        ggplot2_3.3.5      geosphere_1.5-14   chron_2.3-56      
[16] assignR_2.1.1      sf_1.0-6           dplyr_1.0.7        isocat_0.2.6       raster_3.5-15     
[21] sp_1.4-6          

loaded via a namespace (and not attached):
  [1] minqa_1.2.4          colorspace_2.0-2     rjson_0.2.21         deldir_1.0-6        
  [5] ellipsis_0.3.2       class_7.3-20         sjlabelled_1.1.8     estimability_1.3    
  [9] parameters_0.16.0    rstudioapi_0.13      httpcode_0.3.0       proxy_0.4-26        
 [13] fansi_1.0.2          mvtnorm_1.1-3        codetools_0.2-18     splines_4.1.2       
 [17] knitr_1.37           sjmisc_2.8.9         geojsonlint_0.4.0    ade4_1.7-18         
 [21] jsonlite_1.7.3       nloptr_2.0.0         ggeffects_1.1.1      broom_0.7.12        
 [25] png_0.1-7            rgeos_0.5-9          effectsize_0.6.0.1   compiler_4.1.2      
 [29] httr_1.4.2           sjstats_0.18.1       emmeans_1.7.2        backports_1.4.1     
 [33] assertthat_0.2.1     Matrix_1.4-0         cli_3.1.1            adehabitatHR_0.4.19 
 [37] tools_4.1.2          coda_0.19-4          gtable_0.3.0         glue_1.6.1          
 [41] V8_4.0.0             Rcpp_1.0.8           carData_3.0-5        cellranger_1.1.0    
 [45] vctrs_0.3.8          crul_1.2.0           sjPlot_2.8.10        ebirdst_0.3.3       
 [49] nlme_3.1-155         iterators_1.0.13     insight_0.15.0       xfun_0.29           
 [53] lme4_1.1-27.1        lifecycle_1.0.1      terra_1.5-17         MASS_7.3-55         
 [57] scales_1.1.1         parallel_4.1.2       mvnfast_0.2.7        curl_4.3.2          
 [61] gridExtra_2.3        ggspatial_1.1.5      stringi_1.7.6        jsonvalidate_1.3.2  
 [65] bayestestR_0.11.5    maptools_1.1-2       SDMetrics_0.0.0.9000 foreach_1.5.2       
 [69] e1071_1.7-9          boot_1.3-28          RgoogleMaps_1.4.5.3  rlang_1.0.0         
 [73] pkgconfig_2.0.3      bitops_1.0-7         lattice_0.20-45      tidyselect_1.1.1    
 [77] plyr_1.8.6           magrittr_2.0.2       R6_2.5.1             generics_0.1.2      
 [81] DBI_1.1.2            pillar_1.7.0         foreign_0.8-82       withr_2.4.3         
 [85] units_0.8-0          stars_0.5-5          datawizard_0.2.3     abind_1.4-5         
 [89] tibble_3.1.6         performance_0.8.0    modelr_0.1.8         crayon_1.4.2        
 [93] car_3.0-12           KernSmooth_2.23-20   utf8_1.2.2           jpeg_0.1-9          
 [97] grid_4.1.2           CircStats_0.2-6      classInt_0.4-3       pbmcapply_1.5.0     
[101] xtable_1.8-4         adehabitatMA_0.3.14  munsell_0.5.0        adehabitatLT_0.3.25 
