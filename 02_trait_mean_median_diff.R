#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: compare mean and median of traits ----
# this will produce the Figure S1

## Contacts ----
# Chongyi Jiang (chongyi.jiang@uni-jena.de)
## Create dates: 2024-09-24 ----
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
# [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
# [6] readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0
# [11] Taxonstand_2.4  pbapply_1.7-2  
# 
# loaded via a namespace (and not attached):
# [1] rstudioapi_0.15.0 magrittr_2.0.3    hms_1.1.3         tidyselect_1.2.0 
# [5] munsell_0.5.0     timechange_0.2.0  colorspace_2.1-0  R6_2.5.1         
# [9] rlang_1.1.2       fansi_1.0.6       tools_4.2.0       parallel_4.2.0   
# [13] grid_4.2.0        gtable_0.3.4      utf8_1.2.4        cli_3.6.2        
# [17] withr_2.5.2       lifecycle_1.0.4   tzdb_0.4.0        vctrs_0.6.5      
# [21] glue_1.6.2        stringi_1.8.3     compiler_4.2.0    pillar_1.9.0     
# [25] generics_0.1.3    scales_1.3.0      pkgconfig_2.0.3  

library(tidyverse)

rm(list = ls())
trait_mean <- read_delim("results/02_traits/BE_species_traits_both.txt") %>% 
    dplyr::select(BEname, usedPloidy, usedGS, multiPloidy, multiGS)
trait_median <- read_delim("results/02_traits/BE_species_traits_both_median.txt") %>% 
    dplyr::select(BEname, usedPloidy, usedGS, multiPloidy, multiGS)
names(trait_mean) <- c("BEname", "meanPloidy", "meanGS", "multiPloidy", "multiGS")
names(trait_median) <- c("BEname", "medianPloidy", "medianGS", "multiPloidy", "multiGS")

GS_diff <- trait_mean %>% 
    left_join(trait_median) %>% 
    filter(multiGS == "Y") %>% 
    mutate(diff = meanGS - medianGS) %>% 
    dplyr::select(BEname, meanGS, medianGS, diff) %>% 
    rename(mean = meanGS, median = medianGS) %>% 
    mutate(traits = "GS_diff")
mean(GS_diff$diff==0) # percentage of GS difference equal to 0: 66%

PL_diff <- trait_mean %>% 
    left_join(trait_median) %>% 
    filter(multiPloidy == "Y") %>% 
    mutate(diff = meanPloidy - medianPloidy) %>% 
    dplyr::select(BEname, meanPloidy, medianPloidy, diff) %>% 
    rename(mean = meanPloidy, median = medianPloidy) %>% 
    mutate(traits = "PL_diff")
mean(PL_diff$diff==0) # percentage of PL difference equal to 0: 89%

trait_all <- rbind(GS_diff, PL_diff)

trait_all$traits <- factor(trait_all$traits, 
                           levels = c("GS_diff", "PL_diff"),
                           labels = c("Genome Size (2C/pg)", "Ploidy Level"))

## the big value of difference in gneome size are those species with large genome size
## the relative small value of difference in ploidy is because most species are just 
# diploidy, and the mean and median is very small difference
ggplot(trait_all, aes(x = diff)) +
    geom_histogram() +  
    facet_wrap(~traits, scales = "free") + 
    theme_bw() + 
    xlab("Difference between mean and median") +
    ylab("# Species") + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"))

ggsave("results/06_manuscript/FigureS1_mean_median_different.png", 
       width = 5.8, height = 3.2)
