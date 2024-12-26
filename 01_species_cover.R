#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes ----
# 1. find out what species are recorded by Biodiversity Exploratories (BE)
# 2. standardize the species name
# 3. calculate the mean cover (over 2008~2020) of each species in each plot
# 4. subset resident species in each plot

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


library(Taxonstand)
library(tidyverse)
rm(list=ls())

output_dir <- "results/01_species_cover/"
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)

## 1.what species are recorded by BE ----
# 27386_2_Dataset was downloaded from Biodiversity Exploratories online database
# (www.biodiversity-exploratories.de/en/public-data-bexis/, dataset ID: 27386)

sp_records <- read.delim("data/BE/27386_2_Dataset/27386_2_data.txt", 
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
    mutate(Species = str_trim(Species), # in case there is hidden space at the end
           species_replaced = str_replace_all(Species, "_", " ")) %>% 
    dplyr::select(Useful_EP_PlotID, Year, Species, Cover, species_replaced) %>%
    rename(PlotID = Useful_EP_PlotID)

species_BE <-  sp_records %>% distinct(species_replaced) %>% pull()
length(species_BE) # in total 375 species


## 2. standardize the species name ----
# tpl_species <- TPL(species_BE)
### save the result, so that we don't have to run this step every time;
# output <- paste0(output_dir, "BE_BiolFlor_Taxonstand.txt")
# write.table(tpl_species,
#             col.names = TRUE, row.names = FALSE, sep = "\t",
#             file = output)

tpl_species <- read.delim("results/01_species_cover/BE_BiolFlor_Taxonstand.txt",
                          stringsAsFactors = FALSE,
                          sep = "\t", header = TRUE)
accepted <- tpl_species %>% filter(Taxonomic.status == "Accepted")
nrow(accepted) # 310 species names accepted;

# 375-310=65 unaccepted records
unaccepted <- tpl_species %>% 
    filter(Taxonomic.status != "Accepted") %>% 
    dplyr::select(Taxon, Taxonomic.status, New.Genus, New.Species)
sum(unaccepted$New.Species %in% c("sp", "sp."), na.rm = TRUE) # 34 records are unspecific
sum(is.na(unaccepted$New.Species)) # 1 record species is NA; this record is 'Rhinanthus aggr.'

# 65-34-1=30 records left; manual check will focus on these records
manually_check <- unaccepted %>% filter(New.Species != 'sp', 
                                        New.Species != 'sp.', 
                                        !is.na(New.Species))
# output <- paste0(output_dir, "manually_check.txt")
# write.table(manually_check,
#             col.names = TRUE, row.names = FALSE,
#             sep = "\t", file = output)

# since we draw trait estimates from BiolFlor and KEW,
# I search each of these records in BiolFlor, KEW, and the plant list. 
# If we find this record in one of the databases, keep the BE name; otherwise use TPL name;

# Manually checked results in manually_checked.txt
manually_checked <- read.delim("results/01_species_cover/manually_checked.txt", 
                               sep = "\t", header = TRUE, stringsAsFactors = FALSE)
accepted_part <- accepted %>% 
    mutate(useful_name = Taxon) %>% 
    dplyr::select(Taxon, Taxonomic.status, useful_name)

uncertain_part1 <- unaccepted %>% 
    filter(New.Species == "sp" | New.Species == "sp." | is.na(New.Species)) %>% 
    mutate(useful_name = Taxon, Taxonomic.status = "Uncertain") %>% 
    dplyr::select(Taxon, Taxonomic.status, useful_name)

uncertain_part2 <- manually_checked %>% 
    filter(check == "uncertain") %>% 
    mutate(useful_name = Taxon, Taxonomic.status = "Uncertain") %>% 
    dplyr::select(Taxon, Taxonomic.status, useful_name)

keepBE <- manually_checked %>% 
    filter(str_detect(check, "keepBE")) %>% 
    mutate(useful_name = Taxon, Taxonomic.status = check) %>% 
    dplyr::select(Taxon, Taxonomic.status, useful_name)

keepTPL <- manually_checked %>% 
    filter(str_detect(check, "keepTPL")) %>% 
    mutate(useful_name = paste(New.Genus, New.Species),
           Taxonomic.status = "keepTPL") %>% 
    dplyr::select(Taxon, Taxonomic.status, useful_name)

all_sp_name <- rbind(accepted_part, keepBE, keepTPL, uncertain_part1, uncertain_part2) %>% 
    rename(species_replaced = Taxon,
           taxon_status = Taxonomic.status)
sp_records_checked <- all_sp_name %>% left_join(sp_records, by = "species_replaced")

# sp_records_checked$Species; this column was the very original name of BE species (except removing potential blank at the end )
# sp_records_checked$species_replaced; this column replaced the '_' in sp_records_checked$Species with ' '
# sp_records_checked$useful_name; this column is the species name we used for downstream analysis

# there are two records using "_" at the very beginning; while that should be "-"; 
# explicitly change them here
sp_records_checked$useful_name[sp_records_checked$useful_name == "Capsella bursa pastoris"] <- "Capsella bursa-pastoris"
sp_records_checked$useful_name[sp_records_checked$useful_name == "Silene flos cuculi"] <- "Silene flos-cuculi"

# some records are accepted, but there are "aggr." in the end; 
# here explicitly removed the "aggr."
# "Ononis repens spinosa aggr.", "Rhinanthus.aggr." and "Bromus hordeaceus aggr.incl B commutatus" are not included
index <- str_detect(sp_records_checked$useful_name, "aggr.$") & sp_records_checked$taxon_status %in% c("Accepted", "keepBE_no_aggr.")
sp_records_checked$useful_name[index] <- str_replace(sp_records_checked$useful_name[index], " aggr.$", "")

# some records have 'cf' in between, remove them explicitly
index <- str_detect(sp_records_checked$useful_name, " cf ") & sp_records_checked$taxon_status %in% c("Accepted")
unique(sp_records_checked$useful_name[index])
sp_records_checked$useful_name[sp_records_checked$useful_name == "Galeopsis cf bifida"] <- "Galeopsis bifida"
sp_records_checked$useful_name[sp_records_checked$useful_name == "Allium cf oleraceum"] <- "Allium oleraceuma"
sp_records_checked$useful_name[sp_records_checked$useful_name == "Pinus cf mugo"] <- "Pinus mugo"
# "cf" within "Euphorbia cf helioscopia" and "Trifolium cf montanum" are kept, 
# because there are "Euphorbia helioscopia" and "Trifolium montanum" in the list.

names(sp_records_checked)
# [1] "species_replaced" "taxon_status"     "useful_name"      "PlotID"           "Year"             "Species"         
# [7] "Cover"
# "Species" is the original name of BE records (any empty space at the end of the name was removed)
# "species_replaced" is the '_' in the name replaced by " "
# "useful_name" is name after standardized


## 3. calculate the mean cover (over 2008~2020) of each species in each plot ----
output <- paste0(output_dir, "all_sp_yearly_cover_each_plot_long.txt")
write.table(sp_records_checked,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            file = output)
sp_mean_cover_each_plot0 <- sp_records_checked %>% 
    group_by(useful_name, PlotID) %>% 
    summarise(cover_mean = mean(Cover, na.rm = TRUE))
sp_name_cols <- sp_records_checked %>% 
    dplyr::select(-c(PlotID, Year, Cover)) %>% 
    distinct()
sp_mean_cover_each_plot <- sp_mean_cover_each_plot0 %>% 
    left_join(sp_name_cols, by = 'useful_name')

output <- paste0(output_dir, "all_sp_mean_cover_each_plot_long.txt")
write.table(sp_mean_cover_each_plot, sep = "\t",
            col.names = TRUE, row.names = FALSE, 
            file = output)

# seperate the regions
f_seperate <- function(df, regions, output_dir, sp_type){
    # sp_type <- "all_sp"
    df_regions <- vector(mode = "list")
    for (region in regions){
        # region <- "AEG"
        df_regions[[region]] <- df %>% filter(str_detect(PlotID, region)) 
        region_mean_cover <- paste0(output_dir, region, "_cover_mean_", sp_type, ".txt")
        write.table(df_regions[[region]], sep = "\t", 
                    col.names = TRUE, row.names = FALSE,
                    file = region_mean_cover)
        print(ncol(df_regions[[region]]) - 1) # the first column is plotID
        
        region_not_zero <- paste0(output_dir, region, "_cover_mean_not_zero_", sp_type, ".txt")
        # species present at least in one of the plot in the focal region
        plot_column <- df_regions[[region]][1]
        cover_columns <- df_regions[[region]] %>% 
            dplyr::select(-1) %>% 
            dplyr::select(where(~ sum(., na.rm = TRUE) > 0))
        df_regions[[region]] <- bind_cols(plot_column, cover_columns)
        cat(sprintf("\"%s\" \"%i\"\n", region, ncol(cover_columns)))
        
        write.table(df_regions[[region]], sep = "\t",
                    col.names = TRUE, row.names = FALSE,
                    file = region_not_zero)
    }
    
    return(df_regions)
}

df <- sp_mean_cover_each_plot %>% 
    dplyr::select(useful_name, PlotID, cover_mean) %>% 
    pivot_wider(names_from = useful_name, values_from = cover_mean) %>% 
    arrange(PlotID)
cover_mean <- f_seperate(df, regions, output_dir, "all_sp")
# "AEG" "225"
# "HEG" "251"
# "SEG" "214"


## 4. subset resident species in each plot  ----
# keep only species records that found at least 6 times in a particular plot
a <- sp_records_checked %>% filter(Cover > 0) %>% 
    group_by(useful_name, PlotID) %>% 
    summarise(n = n()) %>% 
    filter(n >= 6)
resident_IDs <- unique(paste0(a$useful_name, "__", a$PlotID))

resident_sp_cover_each_plot0 <- sp_records_checked %>% 
    mutate(ID = paste0(useful_name, "__", PlotID)) %>% 
    filter(ID %in% resident_IDs) %>% 
    # left_join(resident_sp_name_cols, by = "useful_name") %>%
    dplyr::select(ID, Year, Cover) %>% 
    separate(col = ID, into = c("useful_name", "PlotID"), sep = "__") %>% 
    arrange(PlotID) 

output <- paste0(output_dir, "resident_sp_yearly_cover_each_plot_long.txt")
write.table(resident_sp_cover_each_plot0,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            file = output)

resident_sp_mean_cover_each_plot0 <- sp_records_checked %>% 
    mutate(ID = paste0(useful_name, "__", PlotID)) %>% 
    filter(ID %in% resident_IDs) %>% 
    group_by(useful_name, PlotID) %>% 
    summarise(cover_mean = mean(Cover, na.rm = TRUE))

resident_sp_name_cols <- sp_records_checked %>% 
    dplyr::select(-c(PlotID, Year, Cover)) %>% 
    distinct()
resident_sp_mean_cover_each_plot <- resident_sp_mean_cover_each_plot0 %>% 
    left_join(resident_sp_name_cols, by = "useful_name")
output <- paste0(output_dir, "resident_sp_mean_cover_each_plot_long.txt")
write.table(resident_sp_mean_cover_each_plot,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            file = output)
resident_df <- resident_sp_mean_cover_each_plot %>% 
    dplyr::select(useful_name, PlotID, cover_mean) %>% 
    pivot_wider(names_from = useful_name, values_from = cover_mean) %>% 
    arrange(PlotID)
cover_mean <- f_seperate(resident_df, regions, output_dir, "resident_sp")
# "AEG" "156"
# "HEG" "144"
# "SEG" "115"
