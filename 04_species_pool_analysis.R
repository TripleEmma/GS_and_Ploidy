## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: perform species level analysis on all species dataset

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

### species level pool-together analysis
library(tidyverse)
library(phytools)
library(phylolm)
library(caper)
library(ape)
library(geiger)

## species preference for all plots
# species relative cover and proxies; the change is the LUI
rm(list = ls())
phylo_trees <- function(BE_sp){
    # load the whole tree
    tree_daphne <- read.tree(file = "data/phylogeny/DaPhnE_01.tre")
    tree_corrected <- read.delim(file ="data/phylogeny/tree_Taxonstand_corrected.txt", 
                                 stringsAsFactors = FALSE)
    tree_daphne$tip.label <- tree_corrected$final
    
    whole_tree <- tree_daphne
    wanted_tree <- drop.tip(whole_tree, setdiff(whole_tree$tip.label, BE_sp))
    
    ## there are some duplications because of standardizing the names
    duplicated_list <- wanted_tree$tip.label[duplicated(wanted_tree$tip.label)]
    tree_corrected_copy <- tree_corrected %>% 
        mutate(final = if_else(final %in% duplicated_list, original_name, final))
    whole_tree$tip.label <- tree_corrected_copy$final
    wanted_tree <- drop.tip(whole_tree, setdiff(whole_tree$tip.label, BE_sp))
    return(wanted_tree)
}


df_all_plot_cover <- read_delim("results/01_species_cover/resident_sp_mean_cover_each_plot_long.txt") %>% 
    dplyr::select(c(useful_name, PlotID, cover_mean)) %>% 
    pivot_wider(names_from = useful_name, values_from = cover_mean) %>% 
    arrange(PlotID)
df_sp_relative_cover <- as.data.frame(apply(df_all_plot_cover[, -1], 2,
                                            function(x){x/sum(x, na.rm = TRUE)})) %>% 
    mutate(PlotID = df_all_plot_cover$PlotID) %>% 
    dplyr::select(PlotID, everything())

df_proxies <- read_delim("results/03_proxies/global_all_proxies_sp.txt") # replace this

# preference:
proxy_preferences <- vector(mode = "list")
for (i in seq(2, length(df_proxies))){
    # i <- 2
    proxy_preference <- as.data.frame(apply(df_sp_relative_cover[, -1], 2, 
                                            function(x){sum(x * df_proxies[, i], na.rm = TRUE)})) %>% 
        rownames_to_column("species")
    names(proxy_preference)[2] <- names(df_proxies)[i]
    proxy_preferences[[i-1]] <- proxy_preference
}
proxy_preference <- reduce(proxy_preferences, inner_join, by = "species")

to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
proxies <- setdiff(names(df_proxies)[-1], to_remove)
trait_df <- read_delim("results/02_traits/regions_traits_resident_sp.txt") %>% 
    dplyr::select(-region) %>% 
    distinct()

# ------------------------------------- gs
gs_analysis <- trait_df %>% 
    # filter(!is.na(GS), !is.na(PL)) %>% 
    filter(!is.na(GS)) %>% 
    dplyr::select(c(species, GS)) %>% 
    mutate(log_GS = log(GS),
           logGS_scaled = scale(log_GS))
proxy_preference_gs <- proxy_preference %>% 
    filter(species %in% gs_analysis$species) %>% 
    inner_join(gs_analysis, by = "species") %>% 
    column_to_rownames("species")

signal_tree <- phylo_trees(rownames(proxy_preference_gs))
name_Ok <- name.check(signal_tree, proxy_preference_gs)
if (name_Ok != 'OK'){stop("The species name is not matched")}
proxy_preference_gs <- proxy_preference_gs %>% rownames_to_column("species")
proxy_preference_gs <- proxy_preference_gs[match(signal_tree$tip.label, proxy_preference_gs$species), ]
gs_signal <- set_names(proxy_preference_gs$GS, proxy_preference_gs$species)
rownames(proxy_preference_gs) <- NULL
tree_target <- unroot(signal_tree)
signal_lambda <- phylosig(tree_target, gs_signal, 
                          method = "lambda", test = TRUE)
signalL <- signal_lambda$lambda
signalL_pValue <- signal_lambda$P

# not separating regions
GS_regression <- data.frame(env_proxy = rep(NA, length(proxies)),
                            signalL = NA, signalL_p = NA,
                            estimate = NA, pvalue = NA, 
                            upper = NA, lower = NA)

for (n in seq(1, length(proxies))){
    # n <- 10
    proxy <- proxies[n]
    df_tmp <- proxy_preference_gs %>% 
        dplyr::select(species, all_of(proxy), logGS_scaled) %>% 
        rename(proxy = 2) %>% 
        mutate(proxy_scaled = scale(proxy)) %>% 
        column_to_rownames("species")
    
    test_signal <- setNames(df_tmp$proxy, rownames(df_tmp))
    signal_lambda <- phylosig(tree_target, test_signal, 
                              method = "lambda", test = TRUE)
    signalL <- signal_lambda$lambda
    signalL_pValue <- signal_lambda$P
    
    target_mod <- phylolm(proxy_scaled~logGS_scaled, 
                          phy = tree_target, model = "BM", boot = 1000, 
                          data = df_tmp)
    
    mod_summary <- summary(target_mod)
    coefficient <- matrix(mod_summary$coefficients["logGS_scaled", 
                                                   c("Estimate", "p.value", 
                                                     "lowerbootCI", "upperbootCI")], 
                          nrow = 1)
    
    row_num <- n
    GS_regression$env_proxy[row_num] <- proxy
    
    GS_regression$signalL[row_num]<- signalL
    GS_regression$signalL_p[row_num] <- signalL_pValue
    
    GS_regression$estimate[row_num] <- coefficient[1, 1]
    GS_regression$pvalue[row_num] <- coefficient[1, 2]
    GS_regression$upper[row_num] <- coefficient[1, 3]
    GS_regression$lower[row_num] <- coefficient[1, 4]
}

# ------------------------------------- pl
pl_analysis <- trait_df %>% 
    # filter(!is.na(GS), !is.na(PL)) %>% 
    filter(!is.na(PL)) %>% 
    mutate(PL_scaled = scale(PL)) %>% 
    dplyr::select(c(species, PL, PL_scaled))

proxy_preference_pl <- proxy_preference %>% 
    filter(species %in% pl_analysis$species) %>% 
    inner_join(pl_analysis, by = "species") %>% 
    column_to_rownames("species")

signal_tree <- phylo_trees(rownames(proxy_preference_pl))
name_Ok <- name.check(signal_tree, proxy_preference_pl)
if (name_Ok != 'OK'){stop("The species name is not matched")}
proxy_preference_pl <- proxy_preference_pl %>% rownames_to_column("species")
proxy_preference_pl <- proxy_preference_pl[match(signal_tree$tip.label, 
                                                 proxy_preference_pl$species), ]
tree_target <- unroot(signal_tree)
rownames(proxy_preference_pl) <- NULL

# not separating regions
PL_regression <- data.frame(env_proxy = rep(NA, length(proxies)),
                            signalL = NA, signalL_p = NA,
                            estimate = NA, pvalue = NA, 
                            upper = NA, lower = NA)

for (n in seq(1, length(proxies))){
    # n <- 10
    proxy <- proxies[n]
    df_tmp <- proxy_preference_pl %>% 
        dplyr::select(species, all_of(proxy), PL, PL_scaled) %>% 
        rename(proxy = 2) %>% 
        mutate(proxy_scaled = scale(proxy)) %>% 
        column_to_rownames("species")
    
    test_signal <- setNames(df_tmp$proxy, rownames(df_tmp))
    signal_lambda <- phylosig(tree_target, test_signal, 
                              method = "lambda", test = TRUE)
    signalL <- signal_lambda$lambda
    signalL_pValue <- signal_lambda$P
    
    target_mod <- phylolm(proxy_scaled ~ PL_scaled,
                          phy = tree_target, model = "BM", boot = 1000,
                          data = df_tmp)
    
    mod_summary <- summary(target_mod)
    coefficient <- matrix(mod_summary$coefficients["PL_scaled", 
                                                   c("Estimate", "p.value", 
                                                     "lowerbootCI", "upperbootCI")], 
                          nrow = 1)
    
    row_num <- n
    PL_regression$env_proxy[row_num] <- proxy
    
    PL_regression$signalL[row_num]<- signalL
    PL_regression$signalL_p[row_num] <- signalL_pValue
    
    PL_regression$estimate[row_num] <- coefficient[1, 1]
    PL_regression$pvalue[row_num] <- coefficient[1, 2]
    PL_regression$upper[row_num] <- coefficient[1, 3]
    PL_regression$lower[row_num] <- coefficient[1, 4]
}


write_out_table <- function(trait, sp_table){
    analysis_type <- ifelse(trait == "GS", "pooled_GS_phylolm", "pooled_PL_phylolm")
    output <- paste0("results/04_species_level/resident_sp/", analysis_type, ".txt")
    write.table(sp_table, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}
write_out_table("GS", GS_regression)
write_out_table("PL", PL_regression)
