## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: pull three regions together and perform community level analysis on this dataset.
# (discussion section; quadratic term included; region is taken as random intercept)

## Contacts ----
# Chongyi Jiang (chongyi.jiang@uni-jena.de)
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
# [1] lmerTest_3.1-3  lme4_1.1-35.1   Matrix_1.6-4    nlme_3.1-164    caper_1.0.3    
# [6] mvtnorm_1.2-4   MASS_7.3-60     boot_1.3-28.1   geiger_2.0.11   ggpubr_0.6.0   
# [11] nortest_1.0-4   phylolm_2.6.2   phytools_2.0-3  maps_3.4.2      ape_5.7-1      
# [16] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
# [21] readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] splines_4.2.0           bit64_4.0.5             vroom_1.6.5            
# [4] foreach_1.5.2           carData_3.0-5           expm_0.999-8           
# [7] cellranger_1.1.0        globals_0.16.2          numDeriv_2016.8-1.1    
# [10] pillar_1.9.0            backports_1.4.1         lattice_0.22-5         
# [13] glue_1.6.2              quadprog_1.5-8          phangorn_2.11.1        
# [16] digest_0.6.33           ggsignif_0.6.4          minqa_1.2.6            
# [19] colorspace_2.1-0        pkgconfig_2.0.3         broom_1.0.5            
# [22] listenv_0.9.0           scales_1.3.0            tzdb_0.4.0             
# [25] optimParallel_1.0-2     timechange_0.2.0        combinat_0.0-8         
# [28] generics_0.1.3          car_3.1-2               withr_2.5.2            
# [31] cli_3.6.2               mnormt_2.1.1            crayon_1.5.2           
# [34] magrittr_2.0.3          readxl_1.4.3            future_1.33.1          
# [37] fansi_1.0.6             parallelly_1.36.0       doParallel_1.0.17      
# [40] rstatix_0.7.2           tools_4.2.0             hms_1.1.3              
# [43] lifecycle_1.0.4         munsell_0.5.0           compiler_4.2.0         
# [46] clusterGeneration_1.3.8 rlang_1.1.2             nloptr_2.0.3           
# [49] grid_4.2.0              iterators_1.0.14        rstudioapi_0.15.0      
# [52] subplex_1.8             igraph_1.6.0            gtable_0.3.4           
# [55] codetools_0.2-19        abind_1.4-5             deSolve_1.40           
# [58] R6_2.5.1                bit_4.0.5               future.apply_1.11.1    
# [61] utf8_1.2.4              fastmatch_1.1-4         stringi_1.8.3          
# [64] parallel_4.2.0          Rcpp_1.0.11             vctrs_0.6.5            
# [67] scatterplot3d_0.3-44    tidyselect_1.2.0        coda_0.19-4    


### community level analysis with plots from all three regions
library(tidyverse)
library(phytools)
library(phylolm)
library(ape)
library(geiger)
library(nlme)
library(nortest)
library(ggpubr)
library(lme4)
library(lmerTest)

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

excluded_plots <- function(all_sp_year, trait_df, trait = "GS", percentage = 0.2){
    # all_sp_year <- df
    cover_df <- all_sp_year %>% 
        dplyr::select(-Year)
    
    # trait = "GS"; percentage = 0.2
    if (trait == "GS"){
        region_trait <- trait_df %>% 
            filter(!is.na(GS), !is.na(PL))
    }else{
        region_trait <- trait_df %>% 
            filter(!is.na(PL))
    }
    
    a <- as.data.frame(t(apply(cover_df[-1], 1, function(x){x/sum(x, na.rm = TRUE)})))
    rownames(a) <- cover_df$PlotID
    missing <- a %>% dplyr::select(-all_of(region_trait$species)) # data frame of species without trait
    missing_cover <- as.data.frame(missing %>% rowSums(na.rm = TRUE)) %>% 
        rownames_to_column("PlotID") 
    names(missing_cover)[2] <- 'missingCover'
    excluded_list <- missing_cover %>% filter(missingCover >  percentage) %>% pull(PlotID)
    return(excluded_list)
}

cover_weighted_yearly <- function(all_sp_year, excluded_trait, 
                                  trait_df, trait = "GS", phylo_tree){
    #  all_sp_year <- df
    "%ni%" <- Negate("%in%")
    excluded <- excluded_trait
    a <- all_sp_year %>% 
        filter(PlotID %ni% excluded) %>% 
        dplyr::select(-Year) %>%  
        column_to_rownames("PlotID") %>% 
        filter(rowSums(is.na(.)) != ncol(.))
    
    region_sp <- names(a)
    
    CWM_gs <- function(col, df, brown){
        df$Col <- col
        pgls <- gls(log_GS ~ 1, data = df, # when only intercept, it is the mean
                    # correlation = brown, weights = varFixed(~Col))
                    correlation = brown, weights = ~Col)
        # predicted_GS <- mean(exp(pgls$fitted))
        predicted_GS <- mean(pgls$fitted)
        return(predicted_GS)
    }
    
    if (trait == "GS"){
        region_trait <- trait_df %>% 
            filter(!is.na(GS), !is.na(PL)) %>% 
            dplyr::select(-contains("PL")) %>% 
            mutate(log_GS = log(GS))
        rownames(region_trait) <- region_trait$species
        target_tree <- phylo_trees(region_trait$species)
        name_Ok <- name.check(target_tree, region_trait)
        if (name_Ok != 'OK'){stop("The species name is not matched")}
        region_trait <- region_trait[match(target_tree$tip.label, region_trait$species), ]
        Species <- region_trait$species
        region_trait <- region_trait %>% rename("Species" = species)
        brown_cor <- corBrownian(value = 1, target_tree, form = ~Species)
        
        b <- a %>% dplyr::select(all_of(Species))
        c <- apply(b, 1, function(x){x/sum(x, na.rm = TRUE)})
        d <- apply(c, 2, function(x){ifelse((is.na(x) | x < 0.000001), 0.000001, x)})
        
        e <- as.data.frame(apply(as.data.frame(d), 2, 
                                 function(x){CWM_gs(col = x, df = region_trait, brown = brown_cor)}))
        e <- e %>% rownames_to_column("PlotID")
        names(e)[2] <- "Plot_log_GS"
        
        # calculate the CWM_pl for these plots using the same species set in those plots
        # the same code as those within else{}
        region_trait <- trait_df %>%
            filter(!is.na(GS), !is.na(PL)) %>%
            dplyr::select(-contains("GS"))
        names(region_trait)[2] <- "trait_value"
        region_trait <- region_trait[match(target_tree$tip.label, region_trait$species), ]
        
        f <- as.data.frame(apply(as.data.frame(d), 2, 
                                 function(x){sum(x * region_trait$trait_value, na.rm = TRUE)})) %>% 
            rownames_to_column("PlotID")
        names(f)[2] <- "Plot_pl"
        
        # merge plot_gs and plot_pl together
        e <- e %>% left_join(f, by = "PlotID")
        
    } else {
        region_trait <- trait_df %>% 
            filter(!is.na(PL)) %>% 
            dplyr::select(-contains("GS"))
        names(region_trait)[2] <- "trait_value"
        b <- a %>% dplyr::select(all_of(region_trait$species))
        c <- as.data.frame(apply(b, 1, function(x){x/sum(x, na.rm = TRUE)}))
        if (all(rownames(c)==region_trait$species)){
            e <- as.data.frame(apply(c, 2, function(x){sum(x * region_trait$trait_value, na.rm = TRUE)})) %>%
                rownames_to_column("PlotID")
            names(e)[2] <- "Plot_pl"
        } else {
            stop("somthing wrong!")
        }
        
    }
    
    return(e)
}

pipeline <- function(year, sp_type, trait){
    # sp_type <- "all_sp"
    (trait_file <- paste0("results/02_traits/regions_traits_", sp_type, ".txt"))
    # (trait_file <- paste0("results/02_traits/regions_traits_", sp_type, "_median.txt"))
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t") %>%
        dplyr::select(-region) %>% 
        distinct()
    
    # year <- 2008; trait <- "GS"
    df <- read_delim("results/01_species_cover/all_sp_yearly_cover_each_plot_long.txt") %>% 
        filter(Year == year) %>% 
        dplyr::select(useful_name, PlotID, Year, Cover) %>% 
        pivot_wider(names_from = useful_name, values_from = Cover)
    # print('I am here')
    excluded_trait <- excluded_plots(df, trait_df, trait = trait, percentage = 0.2)
    # print('I am here 2')
    cover_weighted_trait <- cover_weighted_yearly(df, excluded_trait, trait_df, trait = trait, phylo_tree)
    
    return(cover_weighted_trait)
}

get_last_non_na_column <- function(...) {
    last_non_na_col <- tail(names(c(...))[!is.na(c(...))], 1)
    return(last_non_na_col)
}

CWM <- function(sp_type, years, trait){
    # CWM_trait <- vector(mode = "list")
    # trait <- "GS"; years <- 2008:2020; sp_type <- "all_sp"
    trait_all <- purrr::map(years, pipeline, sp_type, trait)
    names(trait_all) <- paste0("Year", years)
    
    all_trait_merge <- Reduce(function(df1, df2) full_join(df1, df2, by = "PlotID"), trait_all)
    if (trait == "GS"){
        headers <- purrr::map(names(trait_all), 
                              function(x){paste0(x, c("_GS", "_PL"))}) %>% unlist()
        names(all_trait_merge) <- c("PlotID", headers)
        last_non_na_column = pmap_chr(all_trait_merge, get_last_non_na_column) %>% str_sub(end = 8)
        CWM_trait <- all_trait_merge %>% 
            rowwise() %>% 
            mutate(CWM_gs = rowMeans(pick(contains("_GS")), na.rm = TRUE),
                   CWM_pl = rowMeans(pick(contains("_PL")), na.rm = TRUE),
                   n = rowSums(!is.na(pick(-c(PlotID, CWM_gs, CWM_pl)))) / 2
            ) %>% 
            ungroup() %>% 
            mutate(endYEAR = last_non_na_column) %>% 
            filter(n >= 6) %>% 
            dplyr::select(PlotID, CWM_gs, CWM_pl, n, endYEAR) %>% 
            arrange(PlotID)
    }
    
    else{
        names(all_trait_merge) <- c("PlotID", names(trait_all))
        last_non_na_column = pmap_chr(all_trait_merge, get_last_non_na_column)
        CWM_trait <- all_trait_merge %>% 
            rowwise() %>% 
            mutate(CWM_pl = rowMeans(pick(where(is.numeric)), na.rm = TRUE),
                   n = rowSums(!is.na(pick(-c(PlotID, CWM_pl))))
            ) %>% 
            ungroup() %>% 
            mutate(endYEAR = last_non_na_column) %>% 
            filter(n >= 6) %>% 
            # dplyr::select(PlotID, CWM_trait, n) %>% 
            dplyr::select(PlotID, CWM_pl, n, endYEAR) %>% 
            arrange(PlotID)
    }
    return (CWM_trait)
}

community_lm <- function(cwm_df, trait){
    # cwm_df <- CWM_gs
    proxy_df <- read_delim("results/03_proxies/global_all_proxies_community.txt")
    proxies_all <- names(proxy_df)[-c(1, 2)]
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    
    lm_results <- data.frame(env_proxy = rep(NA, length(proxies)), 
                             region_plot_Num = NA, resid_normality = NA, 
                             complexed = NA,
                             estimate = NA, pvalue = NA, upper = NA, lower = NA,
                             sq_estimate = NA, sq_pvalue = NA, 
                             sq_upper = NA, sq_lower = NA)
    n = 0
    
    region_ind <- proxy_df 
    
    df_tmp <- merge(cwm_df, region_ind, by = c("PlotID", "endYEAR"))
    # names(df_tmp)[3] <- "cover_weighted_trait"  ###### check here
    
    for(j in seq(length(proxies))){
        n <- n + 1
        # j <- 1
        # trait <- "GS"
        proxy <- proxies[j]
        
        if (trait == "GS"){
            df_tmp2 <- df_tmp %>% dplyr::select(PlotID, CWM_gs, CWM_pl, all_of(proxy))
            names(df_tmp2)[4] <- "proxy"  ###### check here
            df_tmp3 <- df_tmp2 %>% 
                mutate(proxy_scaled = scale(proxy),
                       proxy_scaled_sq = proxy_scaled^2,
                       region = case_when(grepl('AEG', PlotID) ~ 'AEG', 
                                          grepl('HEG', PlotID) ~ 'HEG', 
                                          grepl('SEG', PlotID) ~ 'SEG'))
            
            simple_mod <- lmer(CWM_gs ~ proxy_scaled + CWM_pl + (1|region), data = df_tmp3)
            complex_mod <- lmer(CWM_gs ~ proxy_scaled + proxy_scaled_sq + CWM_pl + (1|region), data = df_tmp3)
        }
        else{
            df_tmp2 <- df_tmp %>% dplyr::select(PlotID, CWM_pl, all_of(proxy))
            names(df_tmp2)[3] <- "proxy"  ###### check here
            df_tmp3 <- df_tmp2 %>% 
                mutate(proxy_scaled = scale(proxy),
                       proxy_scaled_sq = proxy_scaled^2,
                       region = case_when(grepl('AEG', PlotID) ~ 'AEG', 
                                          grepl('HEG', PlotID) ~ 'HEG', 
                                          grepl('SEG', PlotID) ~ 'SEG'))
            simple_mod <- lmer(CWM_pl ~ proxy_scaled + (1|region), data = df_tmp3)
            complex_mod <- lmer(CWM_pl ~ proxy_scaled + proxy_scaled_sq + (1|region), data = df_tmp3)
        }
        deviance_LRT <- 2 * (logLik(complex_mod) -logLik(simple_mod))
        pValue_LRT <- pchisq(deviance_LRT, df = 1, lower.tail = FALSE)
        complexed <- "N"
        # plot_mod <- simple_mod
        if(pValue_LRT < 0.05){
            complexed <- "Y"
            # plot_mod <- complex_mod
        } 
        
        target_mod <- complex_mod
        normality <- lillie.test(residuals(target_mod))
        mod_summary <- summary(target_mod)
        coefficient <- mod_summary$coefficients["proxy_scaled", c(1, 5)]
        confints <- confint(target_mod)["proxy_scaled", ]
        sq_coefficient <- mod_summary$coefficients["proxy_scaled_sq", c(1, 5)]
        sq_confints <- confint(target_mod)["proxy_scaled_sq", ]
        
        # output table:
        lm_results$region_plot_Num[n] <- nrow(df_tmp3)
        lm_results$env_proxy[n] <- proxy
        lm_results$resid_normality[n] <- normality$p.value
        lm_results$complexed[n] <- complexed
        lm_results$estimate[n] <- coefficient[1]
        lm_results$pvalue[n] <- coefficient[2]
        lm_results$upper[n] <- confints[1]
        lm_results$lower[n] <- confints[2]
        lm_results$sq_estimate[n] <- sq_coefficient[1]
        lm_results$sq_pvalue[n] <- sq_coefficient[2]
        lm_results$sq_upper[n] <- sq_confints[1]
        lm_results$sq_lower[n] <- sq_confints[2]
    }
    return(lm_results)
}

write_out_table <- function(trait, sp_type, sp_analysis){
    analysis_type <- ifelse(trait == "GS", "complex_phylogeny_controled_GS_lmer_", "complex_PL_lmer_")
    output <- paste0("justified/global_community/", analysis_type, sp_type, "_yearly_mean_global.txt")
    write.table(sp_analysis, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}

sp_type <- "all_sp"
years <- 2008:2020
CWM_gs <- CWM(sp_type, years, "GS") 
CWM_pl <- CWM(sp_type, years, "PL")

community_gs_lm <- community_lm(CWM_gs, trait = "GS")
community_pl_lm <- community_lm(CWM_pl, trait = "PL")

write_out_table("GS", sp_type, community_gs_lm)
write_out_table("PL", sp_type, community_pl_lm)
