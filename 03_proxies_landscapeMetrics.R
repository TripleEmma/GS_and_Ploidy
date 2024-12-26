## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: tidy up the landscape metrics 

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
rm(list=ls())

output_dir <- "results/03_proxies/landscape/"

# landscape data was downloaded from Biodiversity Exploratories online database
# (www.biodiversity-exploratories.de/en/public-data-bexis/, dataset ID: 25747)

# EP (patch level)
EP_path <- "data/BE/Landscape/EP/"
(EP_pattern <- list.files(path = EP_path, pattern = "\\.txt"))
list_EP <- vector(mode = "list")
for (i in EP_pattern){
    # print(i)
    EP_file <- paste0(EP_path, i)
    region <- unlist(strsplit(i, '_'))[1]
    year <- unlist(strsplit(i, '_'))[2]
    resolution <- unlist(strsplit(i, '_'))[3]
    list_name <- paste0(region, "_", year, "_", resolution)
    tmp_df <- read.delim(EP_file, sep = " ", stringsAsFactors = FALSE)
    names(tmp_df)[1] <- "PlotID"
    names(tmp_df)[2] <- "type"
    row_num <- nrow(tmp_df)
    for(i in seq(row_num)){
        if(!(grepl("grass.dry", tmp_df[i, 2]) | grepl("grassland", tmp_df[i, 2]))){
            tmp_df[i, c(4:13)] <- NA
        }
    }
    list_EP[[list_name]] <- tmp_df
}

EP_dfs <- bind_rows(list_EP, .id = "region_year_resolution")
EP_dfs2 <- EP_dfs %>% 
    mutate(ID = NULL) %>% 
    pivot_longer(-c(region_year_resolution, PlotID, type),
                 names_to = "proxy", values_to = "value")
output <- paste0(output_dir, "EP_all.txt") # plots other than grassland or grassland.dry are assigned 'NA' for each metric
write.table(EP_dfs2, sep = "\t",
            col.names = TRUE, row.names = FALSE, 
            file = output)

# buffer (radius) class level 
buffer_path <- "data/BE/Landscape/buffer/"
(buffer_pattern <- list.files(path = buffer_path, pattern = "\\.txt"))
list_buffer <- vector(mode = "list")
for (i in buffer_pattern){
    # print(i)
    # i <- "Alb_1820_LUagg1_100.txt"
    buffer_file <- paste0(buffer_path, i)
    preffix <- unlist(strsplit(i, '\\.'))[1]
    region <- unlist(strsplit(preffix, '_'))[1]
    year <- unlist(strsplit(preffix, '_'))[2]
    resolution <- unlist(strsplit(preffix, '_'))[3]
    buffer <- unlist(strsplit(preffix, '_'))[4]
    list_name <- paste0(region, "_", year, "_", resolution)
    tmp_df <- read.delim(buffer_file, sep = " ", stringsAsFactors = FALSE)
    region_label <- str_extract(tmp_df[, 1], ".EG")
    region_num <- str_pad(str_extract(tmp_df[, 1], "\\d+"), 2, pad = "0")
    tmp_df[, 1] <- paste0(region_label, region_num)
    names(tmp_df)[1] <- "PlotID"
    if(resolution == "LUagg1"){
        a <- tmp_df %>% dplyr::select(1, starts_with("L."), starts_with("Shannon."), ends_with("Gr"))
    } else {
        a <- tmp_df %>% dplyr::select(1, starts_with("L."), starts_with("Shannon."), ends_with("Gd"))
    }
    a$buffer <- buffer
    if(is.null(list_buffer[[list_name]])){
        list_buffer[[list_name]] <- a
    }else{ # different radius
        list_buffer[[list_name]] <- rbind(list_buffer[[list_name]], a)
    }
}
LUagg1_list <- names(list_buffer) %>% str_detect("_LUagg1") %>% keep(list_buffer, .)
LUagg2_list <- names(list_buffer) %>% str_detect("_LUagg2") %>% keep(list_buffer, .)
LUagg1_dfs <- bind_rows(LUagg1_list, .id = "region_year_resolution")
LUagg2_dfs <- bind_rows(LUagg2_list, .id = "region_year_resolution")


buffer_LUagg1_dfs <- LUagg1_dfs %>%
    pivot_longer(-c(region_year_resolution, PlotID, buffer),
                 names_to = "proxy", values_to = "value")

buffer_LUagg2_dfs <- LUagg2_dfs %>%
    pivot_longer(-c(region_year_resolution, PlotID, buffer),
                 names_to = "proxy", values_to = "value")

output <- paste0(output_dir, "LUagg1_all.txt")
write.table(buffer_LUagg1_dfs, sep = "\t",
            col.names = TRUE, row.names = FALSE, 
            file = output)

output <- paste0(output_dir, "LUagg2_all.txt")
write.table(buffer_LUagg2_dfs, sep = "\t",
            col.names = TRUE, row.names = FALSE, 
            file = output)
