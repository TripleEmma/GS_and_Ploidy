#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: retrieves traits value from kew databases if not available from BiolFlor ----
# 1. cleaning up the Kew dataset 
# 2. assign one estimate for each trait for each species 
# 3. combine trait data retrieved from BiolFlor and Kew, and produce plot to display it
# (if multiple records for a single species, use the median of the records)
# Median is an alternative to average, we could have used the median to represent species with multiple trait records. 

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
# [1] Taxonstand_2.4  pbapply_1.7-2   ggpubr_0.6.0    lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1  
# [7] dplyr_1.1.4     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4  
# [13] tidyverse_2.0.0
#
# 
# loaded via a namespace (and not attached):
# [1] pillar_1.9.0      compiler_4.2.0    tools_4.2.0       lifecycle_1.0.4   gtable_0.3.4     
# [6] timechange_0.2.0  pkgconfig_2.0.3   rlang_1.1.2       cli_3.6.2         rstudioapi_0.15.0
# [11] parallel_4.2.0    withr_2.5.2       generics_0.1.3    vctrs_0.6.5       hms_1.1.3        
# [16] grid_4.2.0        tidyselect_1.2.0  glue_1.6.2        R6_2.5.1          rstatix_0.7.2    
# [21] fansi_1.0.6       carData_3.0-5     car_3.1-2         tzdb_0.4.0        magrittr_2.0.3   
# [26] scales_1.3.0      backports_1.4.1   abind_1.4-5       colorspace_2.1-0  ggsignif_0.6.4   
# [31] utf8_1.2.4        stringi_1.8.3     munsell_0.5.0     broom_1.0.5      

library(tidyverse)
library(ggpubr)
library(Taxonstand)

rm(list=ls())
output_dir <- "results/02_traits/"

## ------- read in KEW information
kew <- read.delim(file="data/Kew/KEW.txt", header = TRUE, sep = "\t", 
                  stringsAsFactors = FALSE, na.strings = c("", "NA"), 
                  encoding='UTF-8')

nrow(kew)
kew$combinedName <- paste(kew$genus, kew$species, sep = " ")
kew$ploidy[which(kew$ploidy == "-")] <- NA
kew$chromosome_num[which(kew$chromosome_num == "-")] <- NA
names(kew)
# some value has unnoticed space at the ends; 
# we have to remove them, otherwise will cause error when transforming into other types
kew$DNA_1C <- str_trim(kew$DNA_1C, side = c("both", "left", "right"))
kew$ploidy <- str_trim(kew$ploidy, side = c("both", "left", "right"))

## ------- read in BE species and biolFlor information
BE_species <- read.delim("results/01_species_cover/all_sp_mean_cover_each_plot_long.txt")
names(BE_species)
keepBE <- BE_species %>% 
    filter(cover_mean > 0, 
           taxon_status != "Uncertain", 
           taxon_status != "Accepted", 
           taxon_status != "keepTPL") %>% 
    distinct(useful_name) %>% 
    dplyr::pull(useful_name)

tpl_kew <- read.delim("results/02_traits/KEW_Taxonstand.txt",
                      stringsAsFactors = FALSE,
                      sep = "\t", header = TRUE)
kew_original_name <- paste(tpl_kew$Genus, tpl_kew$Species)
kew_tpl_name <- paste(tpl_kew$New.Genus, tpl_kew$New.Species)

a <- data.frame(original_name = kew_original_name, 
                tpl_name = kew_tpl_name, 
                status = tpl_kew$Taxonomic.status) %>% 
    mutate(final = if_else(original_name %in% keepBE, original_name, tpl_name))
kew_corrected <- setNames(a$final, a$original_name)
kew2 <- kew %>% mutate(correctedName = kew_corrected[combinedName])

BE_species_list <- BE_species %>% 
    filter(cover_mean > 0, taxon_status != "Uncertain") %>% 
    dplyr::select(useful_name) %>% 
    mutate(BEname = str_trim(useful_name)) %>% 
    distinct()

biolFlor <- read.delim("results/02_traits/BE_species_BiolFlor_only_median.txt",
                       stringsAsFactors = FALSE) %>%
    dplyr::select(BEname, usedPloidy, usedGS, multiPloidy, multiGS) %>% distinct()
BE_species_traits <- BE_species_list %>% left_join(biolFlor, by= "BEname")
# BE_species_traits$kewGS <- 'N'
# BE_species_traits$kew_usedGS <- NA
# BE_species_traits$kewPloidy <- 'N'
# BE_species_traits$kew_usedPloidy <- NA
BE_species_traits <- BE_species_traits %>% 
    mutate(GS_source = if_else(!is.na(usedGS), "BiolFlor", "NA"), 
           PL_source = if_else(!is.na(usedPloidy), "BiolFlor", "NA"))

# ------- select species lack and information 
biolFlor_no_GS <- BE_species_traits %>% 
    filter(is.na(usedGS)) %>% 
    distinct(BEname) %>% 
    pull(BEname)
length(biolFlor_no_GS) # 206 species lacking genome size information until now 

# ------- select species lack ploidy information 
biolFlor_no_ploidy <- BE_species_traits %>%
    filter(is.na(usedPloidy)) %>% 
    distinct(BEname) %>% 
    pull(BEname)
length(biolFlor_no_ploidy) # 83 species lacking ploidy level information until now

# Fill up the trait information obtained from KEW. 
# For species with multiple records for the trait, we used the median of those records.
# -------------- select Kew species with trait information that is not available in BiolFlor and see if there are species with multiple records
kew_GS <- kew2 %>% filter(correctedName %in% biolFlor_no_GS & !(is.na(DNA_1C))) %>% 
    dplyr::select(correctedName, DNA_1C)
kew_GS_dup <- unique(kew_GS$correctedName[duplicated(kew_GS$correctedName)])
write_csv(tibble(BEname = kew_GS_dup), "results/02_traits/duplicated_GS_KEW.txt")

kew_ploidy <- kew2 %>% filter(correctedName %in% biolFlor_no_ploidy & !(is.na(ploidy))) %>% 
    dplyr::select(correctedName, ploidy)
kew_ploidy_dup <- unique(kew_ploidy$correctedName[duplicated(kew_ploidy$correctedName)])
# write_csv(tibble(BEname = kew_ploidy_dup), "results/02_traits/duplicated_PL_KEW.txt")

length(unique(kew_GS$correctedName)); length(unique(kew_ploidy$correctedName)) # number of records that don't have trait value in BiolFlor but in Kew
# 116 and 42

kew_GS_unique <- setdiff(kew_GS$correctedName, kew_GS_dup) # 102 sp that have unique estimate in KEW
kew_ploidy_unique <- setdiff(kew_ploidy$correctedName, kew_ploidy_dup) # 35 sp that have unique estimate in KEW
map_dbl(list(kew_GS_unique, kew_ploidy_unique), length)

# -------------- fill in KEW info
for (i in seq(nrow(BE_species_traits))){
    sp_name <- BE_species_traits$BEname[i]
    
    if (sp_name %in% kew_GS_unique){
        BE_species_traits$usedGS[i] = as.numeric(kew_GS$DNA_1C[kew_GS$correctedName == sp_name]) * 2
        BE_species_traits$GS_source[i] = "KEW"
        BE_species_traits$multiGS[i] = "N"
    }
    if (sp_name %in% kew_GS_dup){
        kew_GS_values <- kew_GS$DNA_1C[kew_GS$correctedName == sp_name]
        BE_species_traits$usedGS[i] <- median(as.numeric(kew_GS_values) * 2)
        BE_species_traits$GS_source[i] <- "KEW"
        BE_species_traits$multiGS[i] = "Y"
    }
    
    if (sp_name %in% kew_ploidy_unique){
        BE_species_traits$usedPloidy[i] = as.numeric(kew_ploidy$ploidy[kew_ploidy$correctedName == sp_name])
        BE_species_traits$PL_source[i] <- "KEW"
        BE_species_traits$multiPloidy[i] = "N"
    }
    if (sp_name %in% kew_ploidy_dup){
        kew_ploidy_values <- kew_ploidy$ploidy[kew_ploidy$correctedName == sp_name]
        BE_species_traits$usedPloidy[i] = median(as.numeric(kew_ploidy_values))
        BE_species_traits$PL_source[i] <- "KEW"
        BE_species_traits$multiPloidy[i] = "Y"
    }
}


output <- paste0(output_dir, "BE_species_traits_both_median.txt")
BE_species_traits$useful_name <- NULL
write.table(BE_species_traits,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            file = output)


