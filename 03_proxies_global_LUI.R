#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: tidy up the land-use intensity data 
# the dataset under processed is standardized across regions (not within region)
# this dataset is for the analysis in the discussion section

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


rm(list=ls())
library(tidyverse)

LUI <- data.frame(PLOTID = character(),
                  YEAR = character(),
                  EXPLO = character(),
                  G_STD = numeric(),
                  M_STD = numeric(),
                  F_STD = numeric(),
                  LUI = numeric())

# LUI dataset was downloaded from BExIS (www.bexis.uni-jena.de/lui/LUICalculation/index)
# data was standardized across regions

folders <- list.files(path = "data/BE/LUI/", pattern = "global_")
for (folder in folders){
    # folder <- "global_2006-2020"
    file_path <- paste0("data/BE/LUI/", folder, "/")
    txts <- list.files(path = file_path, pattern = "\\.txt")
    txt <- txts[grep("^LUI", basename(txts))]
    input_file <- paste0("data/BE/LUI/", folder, "/", txt)
    a <- read.delim(input_file, sep = ",", stringsAsFactors = FALSE)
    LUI <- rbind(LUI, a) %>% 
        filter(grepl("global", EXPLO))
}

names(LUI)[1] <- "PlotID"
region_label <- str_extract(LUI$PlotID, ".EG")
region_num <- str_pad(str_extract(LUI$PlotID, "\\d+"), 2, pad = "0")
LUI$PlotID <- paste0(region_label, region_num)
LUIsorted <- arrange(LUI, YEAR, PlotID)
output_dir <- "results/03_proxies/LUI_global/"
output <- paste0(output_dir, "LUI_global.txt")
write.table(LUIsorted, sep = "\t",
            col.names = TRUE, row.names = FALSE, 
            file = output)

## get overall LUI
filter_cond <- function(endYear){
    year_start <- 2006
    year_end <- as.integer(endYear)
    # year_end <- 2020
    a <- paste(seq(year_start, year_end), collapse = ", ")
    b <- paste0("overall(", a, ")")
    LUI_file <- "results/03_proxies/LUI_global/LUI_global.txt"
    LUI_df <- read.delim(LUI_file, stringsAsFactors = FALSE) %>% 
        filter(YEAR == b) %>%
        dplyr::select(PlotID, G_STD, M_STD, F_STD, LUI)
    print(paste(b, "total", nrow(LUI_df), "plots"))
    LUI_df$endYEAR <- paste0("Year", year_end)
    # LUI_df$years <- as.character(year_end)
    return(LUI_df)
}
overall_years_LUI <- purrr::map(c(2008:2020), filter_cond) %>% bind_rows()

output_dir <- "results/03_proxies/LUI_global/"
output <- paste0(output_dir, "LUI_overall_years.txt")
write.table(overall_years_LUI, sep = "\t",
            col.names = TRUE, row.names = FALSE,
            file = output)


## calculate sd of LUI
sd_LUI <- function(end){
    start_year <- 2006
    end_year <- end
    # Year <<- paste0("separately(",seq(start_year, end_year),")") #?
    Year <- paste0("separately(",seq(start_year, end_year),")")
    print(Year)
    sds <- LUIsorted %>% filter(YEAR %in% Year) %>%
        group_by(PlotID) %>%
        summarise(sd_G_STD = sd(G_STD), sd_M_STD = sd(M_STD),
                  sd_F_STD = sd(F_STD), sd_LUI = sd(LUI))
    sds$endYEAR <- paste0("Year", end)
    a <- sds %>% mutate(EXPLO = case_when(grepl("AEG", PlotID) ~ "global(ALB)",
                                          grepl("HEG", PlotID) ~ "global(HAI)",
                                          grepl("SEG", PlotID) ~ "global(SCH)")) %>%
        dplyr::select(PlotID, endYEAR, EXPLO, everything())
}

all_LUI_sd <- purrr::map(c(2008:2020), sd_LUI) %>% bind_rows()

output <- paste0(output_dir, "LUI_sds.txt")
write.table(all_LUI_sd, sep = "\t",
            col.names = TRUE, row.names = FALSE,
            file = output)
