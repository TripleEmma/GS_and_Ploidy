#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: deals with trait estimates extracted from BiolFlor ----
# 1. cleaning up the Biolflor dataset 
# 2. assign one estimate for each trait for each species 
# (if multiple records for a single species, use the mean of the records)

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
library(Taxonstand)

rm(list=ls())

output_dir = "results/02_traits/"
### 0. download traits from BiolFlor, 
### selected traits: ploidy level, basic chromosome number, DNA content(pg/2C)
## 1. cleaning up the Biolflor dataset ----
biolFlor_raw <- read.csv2("data/BiolFlor/suchergebnis.csv", 
                          sep = ";", header = FALSE, encoding='UTF-8', 
                          stringsAsFactors = FALSE)

# Usually the first column should be the species name. 
# However, for species have several different trait estimates, 
# it will have several rows and the first column of such rows are empty
# Therefore, we first fill the first column when the species name is missing (with the previous species name).
biolFlor_raw$V1[biolFlor_raw$V1 == ""] <- NA
biolFlor_raw <- biolFlor_raw %>% tidyr::fill(V1)

# The downloaded file is supposed to be separate by ';'. 
# However, the nine column which indicate the quality of the trait also used ';'.
for (i in seq_along(nrow(biolFlor_raw))){ # the 'quality column use ';' to separate the content; this causes confusion. 
    if (biolFlor_raw[i, ncol(biolFlor_raw)] != ""){
        biolFlor_raw[i, 9] <- paste(biolFlor_raw[i, 9], biolFlor_raw[i, 10], sep = " |")
        biolFlor_raw[i, 10] <- biolFlor_raw[i, 11]
        biolFlor_raw[i, 11] <- biolFlor_raw[i, 12]
    }
}

biolFlor_raw$V6[2] <- "reference1"
biolFlor_raw$V11[2] <- "reference2"
biolFlor_raw[, ncol(biolFlor_raw)] <- NULL # delete last column
biolFlor_raw <- biolFlor_raw[-1, ]         # delete first row
names(biolFlor_raw) <- biolFlor_raw[1, ]
biolFlor <- biolFlor_raw[-1, ]

# Some species have several different trait value for genome size, 
# but only one value for ploidy, or vice versa. 
# For such species, the trait with only one value show empty information in the table. 
# We now fill this trait with the last trait value.
for (i in 2:(nrow(biolFlor))){ # some species have several records for one trait but 
    # only one record for the other trait, fill up the trait with value of the last record
    for (j in 2:11){
        if (biolFlor[i, j] == "" & biolFlor[i, 1] == biolFlor[i-1, 1]){
            biolFlor[i, j] <- biolFlor[i-1, j]
        }
    }
}
names(biolFlor) <- c("biolFlorName", "ploidyLevel", "basicChromosomeNumber",
                     "chromosomeNumberGermany", "mostCommonChromosomeRace",
                     "reference1", "DNAcontent_pg_2C", "chromosomeNumber", 
                     "quality", "BENNET_LEITCH_2001", "reference2")

# for traits values that are missing, fill in NA
biolFlor <- biolFlor %>% mutate_all(na_if, "")
# View(biolFlor)
# biolFlorNames <- unique(biolFlor$biolFlorName)
# tpl_biolFlor <- TPL(biolFlorNames)
# # save the result, so that we don't have to run this step every time;
# output <- paste0(output_dir, "BiolFlor_Taxonstand.txt")
# write.table(tpl_biolFlor,
#             col.names = TRUE, row.names = FALSE, sep = "\t",
#             file = output)

################################################################################
# focus on BE species; There should be an easier way to reconcile the names....
BE_species <- read.delim("results/01_species_cover/all_sp_mean_cover_each_plot_long.txt",
                         stringsAsFactors = FALSE)
names(BE_species)
BE_species_list <- BE_species %>% 
    filter(cover_mean > 0, taxon_status != "Uncertain") %>% 
    distinct(useful_name) %>% 
    pull(useful_name)

keepBE <- BE_species %>% 
    filter(cover_mean > 0, 
           taxon_status != "Uncertain", 
           taxon_status != "Accepted", 
           taxon_status != "keepTPL") %>% 
    distinct(useful_name) %>% 
    dplyr::pull(useful_name)

tpl_biolFlor <- read.delim("results/02_traits/BiolFlor_Taxonstand.txt",
                           stringsAsFactors = FALSE,
                           sep = "\t", header = TRUE)

biolFlor_original_name <- paste(tpl_biolFlor$Genus, tpl_biolFlor$Species)
biolFlor_tpl_name <- paste(tpl_biolFlor$New.Genus, tpl_biolFlor$New.Species)
a <- data.frame(original_name = biolFlor_original_name, 
                tpl_name = biolFlor_tpl_name, 
                status = tpl_biolFlor$Taxonomic.status) %>% 
    mutate(final = if_else(original_name %in% keepBE, original_name, tpl_name))
biolFlor_corrected <- setNames(a$final, a$original_name)

# Useful_name could be different from the original species name because of removing things like aggr. and 
# using "standard" name from other database;

f_match <- function(biolFlor, BE_species_list, biolFlor_corrected){
    df_db_used <- biolFlor
    df_db_used$BEname <- NA
    keep_row <- rep(NA, nrow(biolFlor))
    
    for (i in seq(nrow(biolFlor))){
        bi <-  unlist(str_split(df_db_used$biolFlorName[i], " "))
        genus <-  bi[1]
        species <-  bi[2]
        db_name <- paste(genus, species, sep = " ") # since there some names like "Acer monspessulanum L."
        df_db_used$BEname[i] <- biolFlor_corrected[db_name]
        
        if (biolFlor_corrected[db_name] %in% BE_species_list){
            keep_row[i] = i
        }
    }
    
    keep_rows <- keep_row[!is.na(keep_row)]
    df_db_used2 <- df_db_used[keep_rows, ]
    
    return (df_db_used2)
}

biolFlor2 <- f_match(biolFlor, BE_species_list, biolFlor_corrected)
nrow(biolFlor2) # 432  # 475 
length(unique(biolFlor2$biolFlorName)) -length(unique(biolFlor2$BEname))
# 293-263 = 30 original records have been lump together. 
# 475 is much more than 293 or 263 because some species have more than one record for each trait;
# and some species has only one trait available

## 2. decide one value for each trait for each species ----
##----------- Deal with ploidy information
ploidy_df <- biolFlor2 %>%  
    filter(!is.na(ploidyLevel)) %>%   ### only records with ploidy information
    dplyr::select(BEname, ploidyLevel, basicChromosomeNumber, 
                  chromosomeNumberGermany, mostCommonChromosomeRace) %>% 
    distinct() # remove duplicated rows; duplicated because there could be several DNA records for the species
nrow(ploidy_df) # 362
length(unique(ploidy_df$BEname)) 
# 253 BE species that have ploidy level information in BiolFlor 

# b <- read.delim("biolflor.corrected", stringsAsFactors = FALSE) %>% pull(1)
# setdiff(unique(ploidy_df$BEname), b)
# setdiff(b, setdiff(unique(ploidy_df$BEname))
# "Elymus repens" "Silene flos-cuculi" and "Cirsium acaule" were different

unique(ploidy_df$ploidyLevel)
## convert ploidy level string to number
# variable name can't start with number
ploidy_df$ploidyLevel[which(ploidy_df$ploidyLevel == "13-ploid")] <- "ploid_13"
ploidy_df$ploidyLevel[which(ploidy_df$ploidyLevel == "14-ploid")] <- "ploid_14"
ploidy_df$ploidyLevel[which(ploidy_df$ploidyLevel == "16-ploid")] <- "ploid_16"
ploidy_df$ploidyLevel[which(ploidy_df$ploidyLevel == "18-ploid")] <- "ploid_18"
ploidy_df$ploidyLevel[which(ploidy_df$ploidyLevel == "20-ploid")] <- "ploid_20"
ploidy_df$ploidyLevel[which(ploidy_df$ploidyLevel == "24-ploid")] <- "ploid_24"
unique(ploidy_df$ploidyLevel)
ploidy_num <- setNames(c(6, 4, 2, 5, 8, 11, 13, 7, 3, 18, 10, 12, 16), unique(ploidy_df$ploidyLevel))

for (i in seq(nrow(ploidy_df))){
    s <- ploidy_df$ploidyLevel[i]
    ploidy_df$ploidyLevel[i] <- ploidy_num[s]
}
ploidy_df <- transform(ploidy_df, ploidyLevel = as.integer(ploidyLevel))
ploidy_df$usedPloidy <- NA

## indicate if there are multiple records for one species
ploidy_df$multiPloidy <- NA
unique_sp <- ploidy_df %>% 
    dplyr::select(BEname, ploidyLevel) %>% 
    distinct() %>% 
    group_by(BEname) %>% 
    mutate(n = n()) %>% 
    filter(n == 1) %>% 
    pull(BEname)
length(unique_sp) # 204 after correcting BiolFlor name; it was 202 previously (meaning some records are renamed to some other names)

for (i in seq(nrow(ploidy_df))){
    if (ploidy_df$BEname[i] %in% unique_sp){
        ploidy_df$multiPloidy[i] <- 'N'
    }
    else {
        ploidy_df$multiPloidy[i] <- 'Y'
    }
}
duplicate_ploidy_names <- ploidy_df %>% filter(multiPloidy == 'Y') %>% dplyr::select(BEname) %>% distinct()
write_csv(duplicate_ploidy_names, "results/02_traits/duplicated_PL_BiolFlor.txt")

unique_ploidy <- ploidy_df %>% filter(multiPloidy == "N") %>% pull(BEname) %>% unique()
# unique_ploidy duplicateds with unique_sp
all.equal(unique_sp, unique_ploidy)
length(unique(ploidy_df$BEname)) - length(unique_ploidy) # 49 species have multiple ploidy records? previously 50
nrow(duplicate_ploidy_names)

# Some species have several ploidy records and we set up several rules to keep one record for each species.
# (1). if only one record indicate it is found in Germany, use this record to represent this species;
# (2). if several records are found in Germany, the mean of these record to represent this species;
# (3). if no records are found in Germany, check if only one record indicates as most common type; if so, use this record to represent this species;
# (4). if multiple records indicate as most common type, calculate the mean of these records;
# (5). if not (1) or (2) or (3) or (4), calculate the mean of the records

## deal with species with several records
multi_ploidy <- ploidy_df %>% filter(multiPloidy == "Y")
multi_ploidy$dummyGermany <- ifelse(multi_ploidy$chromosomeNumberGermany == "True", 1, 0)
multi_ploidy$dummyMostCommon <- ifelse(multi_ploidy$mostCommonChromosomeRace == "True", 1, 0)
# View(multi_ploidy)

# case(1) and (2)
Germany_ploidy <- (multi_ploidy %>% group_by(BEname) %>%
                       summarise(Germany = sum(dummyGermany)) %>%
                       filter(Germany >= 1))$BEname
length(Germany_ploidy) # 46

# case(3) and (4)
most_common_ploidy <- (multi_ploidy %>% 
                           filter(!BEname %in% Germany_ploidy) %>%
                           group_by(BEname) %>% 
                           summarise(mostCommon = sum(dummyMostCommon)) %>%
                           filter(mostCommon >= 1))$BEname
length(most_common_ploidy) # 0

# case(5)
Germany_most_common <- c(Germany_ploidy, most_common_ploidy)
left_ploidy <- multi_ploidy %>% filter(!BEname %in% Germany_most_common)
length(unique(left_ploidy$BEname)) # 3

## ------------------- 
mean_Germany_ploidy <- multi_ploidy %>% 
    # filter(BEname %in% Germany_ploidy) %>% 
    filter(dummyGermany == 1) %>%  # only records (of each species) occur in Germany are kept
    group_by(BEname) %>% 
    summarise(meanPloidy = round(mean(ploidyLevel), 2))

mean_common_ploidy <- multi_ploidy %>% 
    filter(BEname %in% most_common_ploidy) %>%
    filter(dummyMostCommon == 1) %>%
    group_by(BEname) %>%
    summarise(meanPloidy = round(mean(ploidyLevel), 2))

mean_left_ploidy <- left_ploidy %>% 
    group_by(BEname) %>% 
    summarise(meanPloidy = round(mean(ploidyLevel), 2))

mean_ploidy_df <- rbind(mean_Germany_ploidy, mean_common_ploidy, mean_left_ploidy)

mean_ploidy_vector <- setNames(mean_ploidy_df$meanPloidy, mean_ploidy_df$BEname)

#------ fill in the ploidy_df 
for (i in seq(nrow(ploidy_df))){
    if (ploidy_df$BEname[i] %in% names(mean_ploidy_vector)){
        ploidy_df$usedPloidy[i] <- mean_ploidy_vector[[ploidy_df$BEname[i]]]
    }
    else if(ploidy_df$BEname[i] %in% unique_ploidy){
        ploidy_df$usedPloidy[i] <- round(ploidy_df$ploidyLevel[i], 2)
    }
} # final result is in dataframe ploidy_df


##----------- Deal with genome size information
GS_df <- biolFlor2 %>% 
    filter(!is.na(DNAcontent_pg_2C)) %>%
    dplyr::select(BEname, DNAcontent_pg_2C, chromosomeNumber, 
                  quality, BENNET_LEITCH_2001) %>% distinct()
# remove duplicated rows; duplicated because there could be several ploidy records for the species

GS_df <- transform(GS_df, DNAcontent_pg_2C = as.numeric(DNAcontent_pg_2C))
GS_df$usedGS <- NA
length(unique(GS_df$BEname)) # 130 species have genome information in BiolFlor; previously 128

GS_df$multiGS <- NA
unique_sp <- (GS_df %>% dplyr::select(BEname, DNAcontent_pg_2C) %>% 
                  distinct() %>% 
                  group_by(BEname) %>% 
                  mutate(n=n()) %>% filter(n == 1))$BEname
length(unique_sp) # 82 species have unique GS 

# b <- read.delim("biolflor.corrected", stringsAsFactors = FALSE) %>% pull(1)
# setdiff(unique(GS_df$BEname), b)
# setdiff(b, unique(GS_df$BEname))
# "Elymus repens", "Raphanus raphanistrum"; these two specise names were in BE

for (i in seq(nrow(GS_df))){
    if (GS_df$BEname[i] %in% unique_sp){
        GS_df$multiGS[i] <- 'N'
    }
    else {
        GS_df$multiGS[i] <- 'Y'
    }
}
duplicate_GS_names <- GS_df %>% filter(multiGS == 'Y') %>% dplyr::select(BEname) %>% distinct()
write_csv(duplicate_GS_names, "results/02_traits/duplicated_GS_BiolFlor.txt")

# There are two columns in the dataframe indicate if the record value could be used for downstream analysis.
names(GS_df)
table(GS_df$quality)
GS_df[which(grepl("chromosome number known, but does not occur in German|chromosome number missing in original source | the most common number was added", GS_df$quality)),]
table(GS_df$BENNET_LEITCH_2001)

# The criteria we applied to decide what to do when there are several records for one species are following:
# (1). if several records are found and only one record indicates "occurring in Germany", use this record to represent the species;
# (2). if several records are found "occurring in Germany", use the mean of these records to represent the species;
# (3). if several records are found and none of them indicate "occurring in Germany", see if it is "TRUE" in `BENNET_LEITCH_2001`; if there are several records indicate "TRUE" in this column, use the mean of these records to represent the species;
# (4). if several records are found and none of them indicate "occurring in Germany" nor "TRUE" in `BENNET_LEITCH_2001`, use the mean of all these records to represent the species; 

multi_GS <- GS_df %>% filter(multiGS == "Y")
length(unique(multi_GS$BEname)) # 48 species have multiple GS records 

# case(1) and (2)
Germany_GS <- (multi_GS %>% group_by(BEname) %>% 
                   filter(grepl("chromosome number known and occurring in Germany", quality)))$BEname
length(unique(Germany_GS)) # 45 species; 
# case(3)

BENNET_GS <- (multi_GS %>% 
                  filter(!BEname %in% Germany_GS) %>%
                  group_by(BEname) %>% 
                  filter(grepl("TRUE", BENNET_LEITCH_2001)))$BEname

length(unique(BENNET_GS)) # 0 

Germany_BENNET_GS <- c(Germany_GS, BENNET_GS)
GS_left <- multi_GS %>% filter(!BEname %in% Germany_BENNET_GS) 
length(unique(GS_left$BEname)) # 3

## ------------------- 
mean_Germany_GS <- multi_GS %>% 
    filter(grepl("chromosome number known and occurring in Germany", quality)) %>%
    group_by(BEname) %>% 
    summarise(meanGS = round(mean(DNAcontent_pg_2C, na.rm = TRUE), 2))

mean_BENNET_GS <- multi_GS %>% 
    filter(BEname %in% BENNET_GS) %>%
    filter(grepl("TRUE", BENNET_LEITCH_2001)) %>%
    group_by(BEname) %>% 
    summarise(meanGS = round(mean(DNAcontent_pg_2C, na.rm = TRUE), 2))

mean_left_GS <- GS_left %>%
    group_by(BEname) %>%
    summarise(meanGS = round(mean(DNAcontent_pg_2C, na.rm = TRUE), 2))

mean_GS_df <- rbind(mean_Germany_GS, mean_BENNET_GS, mean_left_GS)
mean_GS_vector <- setNames(mean_GS_df$meanGS, mean_GS_df$BEname)
##------ fill in the GS_df
for (i in seq(nrow(GS_df))){
    if (GS_df$BEname[i] %in% names(mean_GS_vector)){
        GS_df$usedGS[i] <- mean_GS_vector[[GS_df$BEname[i]]]
    }
    else if (GS_df$multiGS[i] == 'N'){
        GS_df$usedGS[i] <- round(GS_df$DNAcontent_pg_2C[i], 2)
    }
} # final result is in dataframe GS_df

##------ merge the ploidy_df and GS_df together
both_df <- ploidy_df %>% full_join(GS_df, by = "BEname")
output <- paste0(output_dir, "BE_species_BiolFlor_only.txt")
write.table(both_df,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            file = output)

num_sp_ploidy <- length(unique(ploidy_df$BEname)) # 253 species
ggplot(ploidy_df, aes(x = usedPloidy)) +
    geom_histogram() + scale_x_continuous(breaks = seq(0, 20, by = 2)) +
    theme_bw() + xlab("Ploidy level") + ylab("# Species") + 
    annotate("text", x = 15, y = 170, size = 5, 
             label = "Ploidy information from BiolFl") +
    annotate("text", x = 15, y = 155, size = 4, 
             label = paste(num_sp_ploidy, "species with ploidy information"))

num_sp_GS <- length(unique(GS_df$BEname)) # 130 species
ggplot(GS_df, aes(x = usedGS)) +
    geom_histogram() + theme_bw() + 
    xlab("Genome size 2C(pg)") + ylab("# Species") + 
    ylim(0, 65) +
    annotate("text", x = 32, y = 62, size = 5, 
             label = "Genome size information from BiolFl") +
    annotate("text", x = 32, y = 57, size = 4, 
             label = paste(num_sp_GS, "species with genome size information"))

