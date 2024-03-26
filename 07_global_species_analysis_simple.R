## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: pull three regions together and perform species level analysis on this dataset.
# (discussion section, simple linear regression)

## Contacts ----
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
library(phytools)
library(phylolm)
library(ape)
library(nortest)
library(ggpubr)
library(geiger)
library(boot)
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

df_proxies <- read_delim("results/03_proxies/global_all_proxies_sp.txt")

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

#------------------------------------- genome size
gs_analysis <- trait_df %>% 
    filter(!is.na(GS), !is.na(PL)) %>% 
    dplyr::select(c(species, GS, PL)) %>% 
    mutate(log_GS = log(GS))
proxy_preference_gs <- proxy_preference %>% 
    filter(species %in% gs_analysis$species) %>% 
    inner_join(gs_analysis, by = "species") %>% 
    column_to_rownames("species")

signal_tree <- phylo_trees(rownames(proxy_preference_gs))
name_Ok <- name.check(signal_tree, proxy_preference_gs)
if (name_Ok != 'OK'){stop("The species name is not matched")}
proxy_preference_gs <- proxy_preference_gs %>% rownames_to_column("species")
proxy_preference_gs <- proxy_preference_gs[match(signal_tree$tip.label, proxy_preference_gs$species), ]
rownames(proxy_preference_gs) <- NULL
target_tree <- unroot(signal_tree)

GS_regression <- data.frame(env_proxy = rep(NA, length(proxies)), 
                            sp_Num = NA, 
                            resid_normality = NA, 
                            estimate = NA, pvalue = NA, 
                            upper = NA, lower = NA)

for (n in seq(1, length(proxies))){
    # n <- 10
    proxy <- proxies[n]
    df_tmp <- proxy_preference_gs %>% 
        dplyr::select(species, all_of(proxy), GS, PL, log_GS)
    names(df_tmp)[2] <- "proxy"
    df_tmp2 <- df_tmp %>% 
        mutate(proxy_scaled = scale(proxy)) %>% 
        column_to_rownames("species")
    
    simple_mod <- phylolm(log_GS ~ proxy_scaled + PL, data = df_tmp2,
                          phy = target_tree, model = "BM", boot = 10000)
    target_mod <- simple_mod
    
    transformed <- chol(solve(vcv(target_tree))) %*% target_mod$residuals
    normality <- lillie.test(transformed)
    mod_summary <- summary(target_mod)
    coefficient <- mod_summary$coefficients["proxy_scaled", 
                                            c("Estimate", "p.value", 
                                              "lowerbootCI", "upperbootCI")]
    
    row_num <- n
    GS_regression$env_proxy[row_num] <- proxy
    GS_regression$sp_Num[row_num] <- nrow(df_tmp)
    GS_regression$resid_normality[row_num] <- normality$p.value
    GS_regression$estimate[row_num] <- coefficient[1]
    GS_regression$pvalue[row_num] <- coefficient[2]
    GS_regression$upper[row_num] <- coefficient[3]
    GS_regression$lower[row_num] <- coefficient[4]
}

# View(GS_regression)  

#------------------------------------- ploidy
pl_analysis <- trait_df %>% 
    filter(!is.na(PL)) %>% 
    dplyr::select(c(species, PL)) %>% 
    mutate(ploidy = if_else(PL > 2, 1, 0))
proxy_preference_pl <- proxy_preference %>% 
    filter(species %in% pl_analysis$species) %>% 
    inner_join(pl_analysis, by = "species") 
fit_logistic_simple <- function(df_tmp, indices, mod_type){
    resampled_data0 <- df_tmp[indices, ]
    resampled_data <- resampled_data0 %>% 
        mutate(proxy_scaled = scale(proxy)) 
    model <- glm(ploidy ~ proxy_scaled, 
                 family="binomial", data = resampled_data)
    return(coef(model))
}

PL_regression <- data.frame(env_proxy = rep(NA, length(proxies)), 
                            sp_Num = NA, 
                            resid_normality = NA, 
                            estimate = NA, pvalue = NA, 
                            upper = NA, lower = NA)

for (n in seq(1, length(proxies))){
    # n <- 10
    proxy <- proxies[n]
    df_tmp <- proxy_preference_pl %>% 
        dplyr::select(species, all_of(proxy), ploidy)
    names(df_tmp)[2] <- "proxy"
    df_tmp2 <- df_tmp %>% 
        mutate(proxy_scaled = scale(proxy)) 
    
    glm_mod_simple <- glm(ploidy ~ proxy_scaled,
                          family="binomial", data = df_tmp2)
    target_mod <- glm_mod_simple
    mod_summary <- summary(target_mod)
    coefficient <- mod_summary$coefficients["proxy_scaled",
                                            c("Estimate", "Pr(>|z|)")]
    transformed <- target_mod$residuals
    normality <- lillie.test(transformed)
    
    fit_logistic_model <- fit_logistic_simple
    num_bootstrap_samples <- 10000
    boot_results <- boot(data = df_tmp2, statistic = fit_logistic_model, R = num_bootstrap_samples)
    conf_intervals <- t(sapply(1:ncol(boot_results$t), function(i) {
        quantile(boot_results$t[, i], c(0.025, 0.975))
    }))
    upper <- conf_intervals[2, 1]
    lower <- conf_intervals[2, 2]
    
    row_num <- n
    PL_regression$env_proxy[row_num] <- proxy
    PL_regression$sp_Num[row_num] <- nrow(df_tmp)
    PL_regression$resid_normality[row_num] <- normality$p.value
    PL_regression$estimate[row_num] <- coefficient[1]
    PL_regression$pvalue[row_num] <- coefficient[2]
    PL_regression$upper[row_num] <- upper
    PL_regression$lower[row_num] <- lower
}
# View(PL_regression)

write_out_table <- function(trait, sp_type, sp_table){
    analysis_type <- ifelse(trait == "GS", "GS_phylolm_", "PL_glm_")
    output <- paste0("justified/global_species/", analysis_type, sp_type, "_simple_global.txt")
    write.table(sp_table, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}
write_out_table("GS", "resident_sp", GS_regression)
write_out_table("PL", "resident_sp", PL_regression)

