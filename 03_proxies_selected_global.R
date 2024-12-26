## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: put all selected stress indicators into one table
# this dataset is for analysis in the discussion section

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
## LUI
# LUI overall (species level)

LUI_cond <- function(end){
    # end <- 2020
    year_end <- paste0("Year", end)
    LUI_file <- "results/03_proxies/LUI_global/LUI_overall_years.txt"
    LUI_df <- read.delim(LUI_file, stringsAsFactors = FALSE) %>% 
        filter(endYEAR == year_end) %>%
        dplyr::select(PlotID, endYEAR, G_STD, M_STD, F_STD, LUI)
}
LUI_df_sp <- LUI_cond(2020) %>% dplyr::select(-endYEAR)
LUI_df_community <- purrr::map(c(2008:2020), LUI_cond) %>% bind_rows()

# LUI interannual variation
LUI_sd_cond <- function(end){
    year_end <- paste0("Year", end)
    LUI_file2 <- "results/03_proxies/LUI_global/LUI_sds.txt"
    LUI_sd_df_sp <- read.delim(LUI_file2, stringsAsFactors = FALSE) %>%
        filter(endYEAR == year_end) %>%
        dplyr::select(PlotID, endYEAR, sd_G_STD, sd_M_STD, sd_F_STD, sd_LUI)
}
LUI_sd_sp <- LUI_sd_cond(2020) %>% dplyr::select(-endYEAR)
LUI_sd_community <- purrr::map(c(2008:2020), LUI_sd_cond) %>% bind_rows()

## Landscape
# patch
filter_cond <- function(){
    year_resolution <- "2008_LUagg2"
    type <- "grass.dry"
    proxy <- c("Area")
    return (list(year_resolution, type, proxy))
}
patch_cond <- filter_cond()
patch_file <- "results/03_proxies/landscape/EP_all.txt"
patch_df <- read.delim(patch_file, stringsAsFactors = FALSE) %>% 
    filter(grepl(patch_cond[[1]], region_year_resolution)) %>% 
    filter(type == patch_cond[[2]]) %>% 
    filter(proxy %in% patch_cond[[3]]) %>% 
    dplyr::select(PlotID, value) %>% 
    rename("Patch.Area" = value)

# class
filter_cond <- function(){
    year <- "_2008_LUagg2"
    # buffer <- "2000"
    buffer <- "1000"
    proxies <- c("C.patch.area.w.Gd", "C.PX.100.Gd", 
                 "C.Hanski.Gd", "C.NN.dist.Gd", 
                 "C.patch.no.Gd", "C.edge.dens.Gd",
                 "C.shape.ind.Gd")
    return (list(year, buffer, proxies))
}
class_cond <- filter_cond()
class_file <- "results/03_proxies/landscape/LUagg2_all.txt"
class_df <- read.delim(class_file, stringsAsFactors = FALSE) %>% 
    filter(grepl(class_cond[[1]], region_year_resolution)) %>% 
    filter(buffer == class_cond[[2]]) %>% 
    filter(proxy %in% class_cond[[3]]) %>% 
    dplyr::select(PlotID, proxy, value) %>% 
    pivot_wider(names_from = proxy, values_from = value)

proxies <- list(LUI_df_sp, LUI_sd_sp,
                patch_df, class_df) %>% reduce(inner_join, by = 'PlotID')

proxies %>%  # none of the indicator contains NA
    summarise(across(-PlotID, ~sum(is.na(.)))) %>% 
    pivot_longer(cols = everything())
names(proxies)[2:9] <- c("Grazing", "Mowing", "Fertilization", "LUI", 
                         "Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd")
output <- "results/03_proxies/global_all_proxies_sp.txt"
write.table(proxies, sep = "\t",
            col.names = TRUE, row.names = FALSE,
            file = output)

############### For community level analysis
LUI_community <- LUI_df_community %>% left_join(LUI_sd_community, by = c("PlotID", "endYEAR"))
names(LUI_community)[3:10] <- c("Grazing", "Mowing", "Fertilization", "LUI", 
                                "Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd")
fragmentation <- patch_df %>% inner_join(class_df, by = 'PlotID')
proxies_community <- LUI_community %>% left_join(fragmentation, by = 'PlotID') 
output <- "results/03_proxies/global_all_proxies_community.txt"
write.table(proxies_community, sep = "\t",
            col.names = TRUE, row.names = FALSE,
            file = output)

######### calculate difference across regions (Welch's F test) species level df
df_long <- proxies %>% 
    pivot_longer(cols = -PlotID, names_to = "env_proxies") %>% 
    mutate(region = case_when(str_detect(PlotID, "AEG") ~ "Alb", 
                              str_detect(PlotID, "HEG") ~ "Hai",
                              str_detect(PlotID, "SEG") ~ "Sch"))
df_long$env_proxies <- factor(df_long$env_proxies, 
                              levels = unique(df_long$env_proxies))
Welch_F_test  <- df_long %>% 
    group_by(env_proxies) %>% 
    summarise(
        broom::tidy(oneway.test(value ~ region, data = pick(everything())))) %>% 
    arrange(env_proxies) %>% 
    rename(env_proxy = "env_proxies")
output <- "results/06_manuscript/TableS1_F_test_global.txt"
write.table(Welch_F_test, sep = "\t",
            col.names = TRUE, row.names = FALSE,
            file = output)



