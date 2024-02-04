## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: compare species overlap between our study and other study (https://doi.org/10.1111/1365-2435.13127)
# (discussion section)

## Contacts ----
# Chongyi Jiang (chongyi.jiang@uni-jena.de)
## Create dates: 2024-01-27 ----
## Licence: CC BY 4.0 ----

## Information of related software and package versions ----
# R version 4.2.0 (2022-04-22)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.3
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
# 
# locale:
#     [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] boot_1.3-28.1   geiger_2.0.11   ggpubr_0.6.0    nortest_1.0-4   phylolm_2.6.2  
# [6] phytools_2.0-3  maps_3.4.2      ape_5.7-1       lubridate_1.9.3 forcats_1.0.0  
# [11] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0    
# [16] tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] subplex_1.8             Rcpp_1.0.11             mvtnorm_1.2-4          
# [4] lattice_0.22-5          listenv_0.9.0           digest_0.6.33          
# [7] foreach_1.5.2           utf8_1.2.4              parallelly_1.36.0      
# [10] R6_2.5.1                backports_1.4.1         coda_0.19-4            
# [13] pillar_1.9.0            rlang_1.1.2             rstudioapi_0.15.0      
# [16] car_3.1-2               phangorn_2.11.1         Matrix_1.6-4           
# [19] combinat_0.0-8          igraph_1.6.0            munsell_0.5.0          
# [22] broom_1.0.5             compiler_4.2.0          numDeriv_2016.8-1.1    
# [25] pkgconfig_2.0.3         mnormt_2.1.1            optimParallel_1.0-2    
# [28] globals_0.16.2          tidyselect_1.2.0        expm_0.999-8           
# [31] quadprog_1.5-8          codetools_0.2-19        fansi_1.0.6            
# [34] future_1.33.1           tzdb_0.4.0              withr_2.5.2            
# [37] MASS_7.3-60             grid_4.2.0              nlme_3.1-164           
# [40] gtable_0.3.4            lifecycle_1.0.4         magrittr_2.0.3         
# [43] scales_1.3.0            carData_3.0-5           future.apply_1.11.1    
# [46] cli_3.6.2               stringi_1.8.3           ggsignif_0.6.4         
# [49] scatterplot3d_0.3-44    doParallel_1.0.17       generics_0.1.3         
# [52] vctrs_0.6.5             fastmatch_1.1-4         deSolve_1.40           
# [55] iterators_1.0.14        tools_4.2.0             glue_1.6.2             
# [58] hms_1.1.3               abind_1.4-5             parallel_4.2.0         
# [61] timechange_0.2.0        colorspace_2.1-0        rstatix_0.7.2          
# [64] clusterGeneration_1.3.8


library(tidyverse)

fragmented_sp <- readxl::read_excel("data/species_list_ploidy_fragmentation.xlsx") %>% 
    pull(SP)

our_sp <- read_delim("results/02_traits/regions_traits_resident_sp.txt") %>% 
    pull(species) %>% unique()

length(intersect(fragmented_sp, our_sp)) / length(fragmented_sp)
