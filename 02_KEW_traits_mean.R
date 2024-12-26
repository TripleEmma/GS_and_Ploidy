#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: retrieves traits value from kew databases if not available from BiolFlor ----
# 1. cleaning up the Kew dataset 
# 2. assign one estimate for each trait for each species 
# 3. combine trait data retrieved from BiolFlor and Kew, and produce plot to display it
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

biolFlor <- read.delim("results/02_traits/BE_species_BiolFlor_only.txt",
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
# For species with multiple records for the trait, we used the mean of those records.
# -------------- select Kew species with trait information that is not available in BiolFlor and see if there are species with multiple records
kew_GS <- kew2 %>% filter(correctedName %in% biolFlor_no_GS & !(is.na(DNA_1C))) %>% 
    dplyr::select(correctedName, DNA_1C)
kew_GS_dup <- unique(kew_GS$correctedName[duplicated(kew_GS$correctedName)])
write_csv(tibble(BEname = kew_GS_dup), "results/02_traits/duplicated_GS_KEW.txt")

kew_ploidy <- kew2 %>% filter(correctedName %in% biolFlor_no_ploidy & !(is.na(ploidy))) %>% 
    dplyr::select(correctedName, ploidy)
kew_ploidy_dup <- unique(kew_ploidy$correctedName[duplicated(kew_ploidy$correctedName)])
write_csv(tibble(BEname = kew_ploidy_dup), "results/02_traits/duplicated_PL_KEW.txt")

length(unique(kew_GS$correctedName)); length(unique(kew_ploidy$correctedName)) 
# number of records that don't have trait value in BiolFlor but in Kew
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
        BE_species_traits$usedGS[i] <- mean(as.numeric(kew_GS_values) * 2)
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
        BE_species_traits$usedPloidy[i] = mean(as.numeric(kew_ploidy_values))
        BE_species_traits$PL_source[i] <- "KEW"
        BE_species_traits$multiPloidy[i] = "Y"
    }
}

# How do these trait values distribute?
# note that the species are those with cover > 0;
ploidys <- BE_species_traits %>% 
    filter(!is.na(usedPloidy)) %>% 
    dplyr::select(BEname, usedPloidy) %>% 
    distinct()
num_sp_ploidy <- nrow(ploidys) # 295
ggplot(ploidys, aes(x = usedPloidy)) +
    geom_histogram() + scale_x_continuous(breaks = seq(0, 20, by = 2)) + 
    theme_bw() + xlab("Ploidy level") + ylab("# Species") +
    annotate("text", x = 14, y = 150, size = 5, label = "Ploidy estimates from BiolFlor & KEW") +
    annotate("text", x = 14, y = 140, size = 4, label = paste(num_sp_ploidy, "species with ploidy estimate"))

GS <- BE_species_traits %>% 
    filter(!is.na(usedGS)) %>% 
    dplyr::select(BEname, usedGS) %>% 
    distinct()
num_sp_GS <- nrow(GS) # 246
ggplot(GS, aes(x = usedGS)) +
    geom_histogram() + theme_bw() + 
    xlab("Genome size 2C(pg)") + ylab("# Species") + 
    annotate("text", x = 32, y = 62, size = 5, label = "Genome size estimates from BiolFlor & KEW") +
    annotate("text", x = 32, y = 57, size = 4, label = paste(num_sp_GS, "species with genome size estimate"))

output <- paste0(output_dir, "BE_species_traits_both.txt")
BE_species_traits$useful_name <- NULL
write.table(BE_species_traits,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            file = output)

plot_theme <- theme(# legend.position = "none",
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    axis.text=element_text(size = 10),
    axis.title=element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.key.size = unit(0.65, 'cm'),
    # legend.position=c(0.55, 0.797),
    legend.background = element_rect(fill = "transparent"),
    legend.box = NULL)

ploidys <- BE_species_traits %>% 
    filter(!is.na(usedPloidy)) %>% 
    dplyr::select(BEname, usedPloidy, PL_source)
ploidys$PL_source <- factor(ploidys$PL_source, levels = c("BiolFlor", "KEW"))
num_sp_ploidy <- nrow(ploidys)
pl_legend_texts <- c(paste0("BiolFlor: ", round(mean(ploidys$PL_source == "BiolFlor")* 100, 2), "%"),
                     paste0("KEW: ", round(mean(ploidys$PL_source == "KEW")* 100, 2), "%"))

g2 <- ggplot(ploidys, aes(x = usedPloidy, fill = PL_source)) +
    geom_histogram(alpha = 0.6, binwidth = 1,) + theme_bw() + 
    scale_x_continuous(breaks = seq(0, max(ploidys$usedPloidy), 2)) +
    # guides(fill=guide_legend(title="Database")) + 
    scale_fill_manual(values = c("#FF9C00", "#008285"),
                      labels = c(pl_legend_texts)) +
    xlab("Ploidy level") + ylab("") +
    annotate("text", x = 13, y = 138.6, size = 4.0, 
             label = paste(num_sp_ploidy, "species with ploidy levle estimate")) +
    plot_theme +
    theme(legend.position=c(0.52, 0.78))

GS <- BE_species_traits %>% 
    filter(!is.na(usedGS)) %>% 
    dplyr::select(BEname, usedGS, GS_source)
GS$GS_source <- factor(GS$GS_source, levels = c("BiolFlor", "KEW"))
num_sp_GS <- nrow(GS)
gs_legend_texts <- c(paste0("BiolFlor: ", round(mean(GS$GS_source == "BiolFlor")* 100, 2), "%"),
                     paste0("KEW: ", round(mean(GS$GS_source == "KEW")* 100, 2), "%"))

g1 <- ggplot(GS, aes(x = usedGS, fill = GS_source)) +
    geom_histogram(alpha = 0.6, bins = 30) + theme_bw() + 
    # guides(fill=guide_legend(title="Database")) + 
    scale_fill_manual(values = c("#FF9C00", "#008285"),
                      labels = c(gs_legend_texts)) +
    xlab("Genome size 2C(pg)") + ylab("# Species") + 
    annotate("text", x = 30, y = 62, size = 4.0, 
             label = paste(num_sp_GS, "species with genome size estimate")) +
    plot_theme +
    theme(legend.position=c(0.452, 0.775))

png("results/06_manuscript/all_sp/Figure1b_all_sp_Trait_histograms.png",
    width = 26, height = 15, units = "cm", res = 600)
ggarrange(g1, g2, nrow = 1)
dev.off()


########## species traits according to regions
rm(list = ls())
traits_df <- read.delim("results/02_traits/BE_species_traits_both.txt", 
                        stringsAsFactors = FALSE, check.names = FALSE)
regions_traits <- function(sp_type){
    cover_dir <- "results/01_species_cover/"
    cover_files <- list.files(path = cover_dir, pattern = paste0("_cover_mean_not_zero_", sp_type))
    regions <- c("AEG", "HEG", "SEG")
    sp <- vector(mode = "list")
    for (i in cover_files){
        region <- unlist(strsplit(i, "_"))[1]
        cover_file <- paste0(cover_dir, i)
        cover_df <- read.delim(cover_file, stringsAsFactors = FALSE, check.names = FALSE)
        sp[[region]] <- names(cover_df)[-1] # the first column is PlotID
    }
    sp_trait <- vector(mode = "list")
    for (region in regions){
        sp_trait[[region]] <- data.frame(species = sp[[region]],
                                         GS = NA, GS_source = NA, 
                                         PL = NA, PL_source = NA)
        GS_vector <- setNames(traits_df$usedGS, traits_df$BEname)
        PL_vector <- setNames(traits_df$usedPloidy, traits_df$BEname)
        GS_source_vector <- setNames(traits_df$GS_source, traits_df$BEname)
        PL_source_vector <- setNames(traits_df$PL_source, traits_df$BEname)
        
        for (j in seq(sp[[region]])){
            sp_record <- sp[[region]][j]
            if (sp_record %in% traits_df$BEname){
                sp_trait[[region]][j, 2] <- GS_vector[[sp_record]]
                sp_trait[[region]][j, 3] <- GS_source_vector[[sp_record]]
                sp_trait[[region]][j, 4] <- PL_vector[[sp_record]]
                sp_trait[[region]][j, 5] <- PL_source_vector[[sp_record]]
            }
        }
    }
    sp_trait_all <- bind_rows(sp_trait, .id = "region")
    output_file <- paste0("results/02_traits/regions_traits_", sp_type, ".txt")
    write.table(sp_trait_all, sep = "\t", 
                row.names = FALSE, output_file)   
}
regions_traits("resident_sp")
regions_traits("all_sp")

##### To draw the same histogram plot for resident species
# "results/02_traits/regions_traits_resident_sp.txt" was created just now (see above)
resident_sp_traits <- read.delim("results/02_traits/regions_traits_resident_sp.txt",
                                 stringsAsFactors = FALSE, sep = "\t")

plot_theme <- theme(# legend.position = "none",
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    axis.text=element_text(size = 10),
    axis.title=element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.key.size = unit(0.65, 'cm'),
    # legend.position=c(0.55, 0.797),
    legend.background = element_rect(fill = "transparent"),
    legend.box = NULL)

GS <- resident_sp_traits %>% filter(!is.na(GS)) %>% 
    dplyr::select(species, GS, GS_source) %>%
    distinct() # get the unique records; some species are present in three regions
GS$GS_source <- factor(GS$GS_source, levels = c("BiolFlor", "KEW"))
num_sp_GS <- nrow(GS) # 172
gs_legend_texts <- c(paste0("BiolFlor: ", round(mean(GS$GS_source == "BiolFlor")* 100, 1), "%"),
                     paste0("KEW: ", round(mean(GS$GS_source == "KEW")* 100, 1), "%"))

g1 <- ggplot(GS, aes(x = GS, fill = GS_source)) +
    geom_histogram(alpha = 0.6, bins = 30) + theme_bw() + 
    # guides(fill=guide_legend(title="Database")) + 
    scale_fill_manual(values = c("#FF9C00", "#008285"),
                      labels = gs_legend_texts) +
    xlab("Genome size 2C(pg)") + ylab("# Species") + 
    annotate("text", x = 23, y = 36.9, size = 4.0, 
             label = paste(num_sp_GS, "species with genome size estimate")) +
    plot_theme +
    theme(legend.position=c(0.4, 0.785))

ploidys <- resident_sp_traits %>% filter(!is.na(PL)) %>% 
    dplyr::select(species, PL, PL_source) %>% 
    distinct()
ploidys$PL_source <- factor(ploidys$PL_source, levels = c("BiolFlor", "KEW"))
num_sp_ploidy <- nrow(ploidys) # 204
pl_legend_texts <- c(paste0("BiolFlor: ", round(mean(ploidys$PL_source == "BiolFlor")* 100, 2), "%"),
                     paste0("KEW: ", round(mean(ploidys$PL_source == "KEW")* 100, 2), "%"))

g2 <- ggplot(ploidys, aes(x = PL, fill = PL_source)) +
    geom_histogram(alpha = 0.6, binwidth = 1,) + theme_bw() + 
    scale_x_continuous(breaks = seq(0, max(ploidys$PL), 2)) +
    # guides(fill=guide_legend(title="Database")) + 
    scale_fill_manual(values = c("#FF9C00", "#008285"),
                      labels = pl_legend_texts) +
    xlab("Ploidy level") + ylab("") +
    annotate("text", x = 11.18, y = 89, size = 4.0, 
             label = paste(num_sp_ploidy, "species with ploidy level estimate")) +
    plot_theme +
    theme(legend.position=c(0.412, 0.791))

png("results/06_manuscript/resident_sp/Figure1b_resident_sp_Trait_histograms.png",
    width = 26, height = 15, units = "cm", res = 600)
ggarrange(g1, g2, nrow = 1)
dev.off()


