## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: perform community level analysis on all species dataset (simple linear)
# cover-weighted genomic traits are response and each environmental proxy is predictor
# cover-weighted genome sizes have been corrected for phylogenetic relatness.
# Plots from all regions are pooled together, and region is fitted as random effect

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

### community level analysis with plots from all three regions
library(tidyverse)
library(phytools)
library(phylolm)
library(ape)
library(geiger)
library(nlme)
library(nortest)
library(lme4)
library(lmerTest)

rm(list = ls())
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
proxy_file <- "results/03_proxies/all_proxies_community.txt"
proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
proxies_all <- names(proxy_df)[-c(1, 2)]
to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
proxies <- setdiff(proxies_all, to_remove)
set.seed(0305)

excluded_plots <- function(all_sp_year, trait_df, trait = "GS", percentage = 0.2){
    missing_cover <- vector(mode = "list")
    excluded_list <- vector(mode = "list")
    # all_sp_year <- df
    for (Region in regions){
        # Region = "AEG"
        cover_df <- all_sp_year %>% 
            filter(grepl(Region, PlotID)) %>% 
            dplyr::select(-Year)
        
        # i <- "AEG_cover_mean_not_zero_resident_sp.txt"; trait = "GS"; percentage = 0.2
        if (trait == "GS"){
            region_trait <- trait_df %>% 
                # filter(region == Region, !is.na(GS), !is.na(PL))
                filter(region == Region, !is.na(GS))
        }else{
            region_trait <- trait_df %>% 
                filter(region == Region, !is.na(PL))
        }
        
        a <- as.data.frame(t(apply(cover_df[-1], 1, function(x){x/sum(x, na.rm = TRUE)})))
        rownames(a) <- cover_df$PlotID
        missing <- a %>% dplyr::select(-all_of(region_trait$species)) # data frame of species without trait
        missing_cover[[Region]] <- as.data.frame(missing %>% rowSums(na.rm = TRUE))
    }
    excluded_list <- purrr::map(missing_cover, function(x){rownames(x)[x[, 1] > percentage]})
    return(excluded_list)
}

## module2: calculate cover-weighted GS, PL for each remaining plots
# for GS, we did phylogenetic control; but not for PL
# get phylogentic tree for each region [gs]
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

cover_weighted_yearly <- function(all_sp_year, excluded_trait, 
                                  trait_df, trait = "GS", phylo_tree){
    plot_trait <- vector(mode = "list")
    "%ni%" <- Negate("%in%")
    for (Region in regions){   
        # Region = "AEG"; all_sp_year <- df
        excluded <- excluded_trait[[Region]]
        cover_df_to_weight <- all_sp_year %>% 
            filter(grepl(Region, PlotID)) %>% 
            filter(PlotID %ni% excluded) %>% 
            dplyr::select(-Year) %>%  
            column_to_rownames("PlotID") %>% 
            filter(rowSums(is.na(.)) != ncol(.))
        region_sp <- names(cover_df_to_weight)
        
        CWM_gs <- function(col, df, brown){
            # df <- region_trait; col <- d$AEG01; brown <- brown_cor
            df$Col <- col
            corrected <- gls(log_GS ~ 1, data = df, # when only intercept, it is the mean
                             correlation = brown)
            
            # predicted_GS <- mean(exp(pgls$fitted))
            GS_corrected <- residuals(corrected) + coef(corrected)
            predicted_GS <- sum(GS_corrected * df$Col)
            return(predicted_GS)
        }
        
        if (trait == "GS"){
            region_trait <- trait_df %>% 
                # filter(region == Region, !is.na(GS), !is.na(PL)) %>% 
                # dplyr::select(-contains("PL")) %>% 
                filter(region == Region, !is.na(GS)) %>% 
                mutate(log_GS = log(GS))
            rownames(region_trait) <- region_trait$species
            target_tree <- phylo_trees(region_trait$species)
            name_Ok <- name.check(target_tree, region_trait)
            if (name_Ok != 'OK'){stop("The species name is not matched")}
            region_trait <- region_trait[match(target_tree$tip.label, region_trait$species), ]
            Species <- region_trait$species
            region_trait <- region_trait %>% rename("Species" = species)
            brown_cor <- corBrownian(value = 1, target_tree, form = ~Species)
            
            b <- cover_df_to_weight %>% dplyr::select(all_of(Species))
            c <- apply(b, 1, function(x){x/sum(x, na.rm = TRUE)})
            d <- apply(c, 2, function(x){ifelse((is.na(x) | x < 0.000001), 0.000001, x)})
            d <- as.data.frame(d)
            if (all(rownames(d) == region_trait$Species)){
                e <- as.data.frame(apply(d, 2, 
                                         function(x){CWM_gs(col = x, df = region_trait, brown = brown_cor)}))
                e <- e %>% rownames_to_column("PlotID")
                names(e)[2] <- "Plot_log_GS"
            }
            else{
                stop('something wrong!!!')
            }
            
        } else {
            region_trait <- trait_df %>% 
                filter(region == Region, !is.na(PL)) %>% 
                dplyr::select(-contains("GS"))
            names(region_trait)[3] <- "trait_value"
            b <- cover_df_to_weight %>% dplyr::select(all_of(region_trait$species))
            c <- as.data.frame(apply(b, 1, function(x){x/sum(x, na.rm = TRUE)}))
            if (all(rownames(c)==region_trait$species)){
                e <- as.data.frame(apply(c, 2, function(x){sum(x * region_trait$trait_value, na.rm = TRUE)})) %>%
                    rownames_to_column("PlotID")
                names(e)[2] <- "Plot_pl"
            } else {
                stop("somthing wrong!")
            }
            
        }
        plot_trait[[Region]] <- e
    }
    return(plot_trait)
}

pipeline <- function(year, sp_type, trait){
    (trait_file <- paste0("results/02_traits/regions_traits_", sp_type, ".txt"))
    # (trait_file <- paste0("results/02_traits/regions_traits_", sp_type, "_median.txt"))
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    # year <- 2008; trait <- "GS"
    df <- read_delim("results/01_species_cover/all_sp_yearly_cover_each_plot_long.txt") %>% 
        filter(Year == year) %>% 
        dplyr::select(useful_name, PlotID, Year, Cover) %>% 
        pivot_wider(names_from = useful_name, values_from = Cover)
    
    excluded_trait <- excluded_plots(df, trait_df, trait = trait, percentage = 0.2)
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
    
    list_of_big_tables <- lapply(trait_all, function(tables) {
        bind_rows(tables)  # Merge and add an optional Source column
    })
    all_trait_merge <- Reduce(function(df1, df2) full_join(df1, df2, by = "PlotID"),
                              list_of_big_tables)
    
    names(all_trait_merge) <- c("PlotID", names(trait_all))
    last_non_na_column = pmap_chr(all_trait_merge, get_last_non_na_column)
    CWM_trait <- all_trait_merge %>% 
        rowwise() %>% 
        mutate(CWM_trait = rowMeans(pick(where(is.numeric)), na.rm = TRUE),
               n = rowSums(!is.na(pick(-c(PlotID, CWM_trait))))
        ) %>% 
        ungroup() %>% 
        mutate(endYEAR = last_non_na_column) %>% 
        filter(n >= 6) %>% 
        dplyr::select(PlotID, CWM_trait, n, endYEAR) %>% 
        arrange(PlotID)
    
    return (CWM_trait)
}

## module3: regression analysis (lm)
# we use bootstrap to get confidence interval; construct boot function
bootstrap_lmer <- function(data, formula, group_var, num_bootstrap, complex_formula = NULL) {
    # Store bootstrap results
    bootstrap_results <- matrix(NA, nrow = num_bootstrap, ncol = 0)
    groups <- unique(data[[group_var]])
    
    # Perform bootstrap
    for (i in 1:num_bootstrap) {
        cat("Bootstrap iteration:", i, "\n")
        
        repeat {
            # Resample data with replacement
            # resampled_data <- data[sample(nrow(data), replace = TRUE), ]
            
            resampled_data_list <- lapply(groups, function(group) {
                subset_data <- data[data[[group_var]] == group, ]
                subset_data[sample(nrow(subset_data), replace = TRUE), ]
            })
            
            # Combine resampled data from all groups
            resampled_data <- do.call(rbind, resampled_data_list)
            
            # Try to fit the model
            fit <- tryCatch({
                if (!is.null(complex_formula)) {
                    lmer(complex_formula, data = resampled_data)
                } else {
                    lmer(formula, data = resampled_data)
                }
            }, error = function(e) {
                NULL
            })
            
            # Check if the model is valid (not singular)
            if (!is.null(fit) && !isSingular(fit)) {
                # Store fixed effect coefficients
                coeffs <- fixef(fit)
                
                # Initialize the matrix with the correct number of columns
                if (ncol(bootstrap_results) == 0) {
                    bootstrap_results <- matrix(NA, nrow = num_bootstrap, ncol = length(coeffs))
                    colnames(bootstrap_results) <- names(coeffs)
                }
                
                bootstrap_results[i, ] <- coeffs
                break
            } else {
                cat("Singular fit or error - resampling...\n")
            }
        }
    }
    
    # Return bootstrap results
    return(as.data.frame(bootstrap_results))
}

calculate_ci <- function(bootstrap_results) {
    apply(bootstrap_results, 2, function(x) {
        quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
    })
}

community_lmer <- function(cwm_df, trait){
    # cwm_df <- CWM_gs
    proxy_df <- read_delim("results/03_proxies/global_all_proxies_community.txt")
    proxies_all <- names(proxy_df)[-c(1, 2)]
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    
    lm_results_simple <- data.frame(env_proxy = rep(NA, length(proxies)), 
                                    estimate = NA, pvalue = NA, 
                                    upper = NA, lower = NA)
    
    lm_results_complex <- data.frame(env_proxy = rep(NA, length(proxies)), 
                                     bic_diff = NA, complexed = NA,
                                     estimate = NA, pvalue = NA, 
                                     upper = NA, lower = NA,
                                     sq_estimate = NA, sq_pvalue = NA, 
                                     sq_upper = NA, sq_lower = NA)
    
    n = 0
    
    region_ind <- proxy_df 
    
    df_tmp <- merge(cwm_df, region_ind, by = c("PlotID", "endYEAR"))
    # names(df_tmp)[3] <- "cover_weighted_trait"  ###### check here
    
    for(j in seq(length(proxies))){
        n <- n + 1
        # j <- 1; trait <- "GS"
        proxy <- proxies[j]
        
        df_tmp2 <- df_tmp %>% dplyr::select(PlotID, CWM_trait, all_of(proxy))
        names(df_tmp2)[3] <- "proxy"  ###### check here
        df_tmp3 <- df_tmp2 %>% 
            mutate(proxy_scaled = scale(proxy),
                   CWM_trait0 = CWM_trait,
                   CWM_trait = scale(CWM_trait0),
                   proxy_scaled_sq = proxy_scaled^2,
                   region = case_when(grepl('AEG', PlotID) ~ 'AEG', 
                                      grepl('HEG', PlotID) ~ 'HEG', 
                                      grepl('SEG', PlotID) ~ 'SEG'))
        
        simple_mod <- lmer(CWM_trait ~ proxy_scaled + (1|region), data = df_tmp3)
        complex_mod <- lmer(CWM_trait ~ proxy_scaled + proxy_scaled_sq + (1|region), data = df_tmp3)
        
        complexed <- "N"
        target_mod <- simple_mod
        
        mod_summary <- summary(target_mod)
        coefficient <- mod_summary$coefficients["proxy_scaled", c(1, 5)]
        
        num_bootstrap_samples <- 1000
        bootstrap_results_simple <- bootstrap_lmer(
            data = df_tmp3,
            formula = CWM_trait ~ proxy_scaled + (1 | region),
            group_var = "region",
            num_bootstrap = num_bootstrap_samples
        )
        conf_intervals <- calculate_ci(bootstrap_results_simple)

        # output simple model table:
        lm_results_simple$env_proxy[n] <- proxy
        
        lm_results_simple$estimate[n] <- coefficient[1]
        lm_results_simple$pvalue[n] <- coefficient[2]
        lm_results_simple$upper[n] <- conf_intervals[1, 2]
        lm_results_simple$lower[n] <- conf_intervals[2, 2]
        
        bic_simple <- BIC(simple_mod)
        bic_complex <- BIC(complex_mod)
        bic_diff <- bic_complex - bic_simple
        # if(pValue_LRT < 0.05){
        if(bic_diff < 0){
            complexed <- "Y"
            target_mod <- complex_mod
            normality <- lillie.test(residuals(target_mod))
            mod_summary <- summary(target_mod)
            coefficient <- mod_summary$coefficients["proxy_scaled", c(1, 5)]
            sq_coefficient <- mod_summary$coefficients["proxy_scaled_sq", c(1, 5)]
            bootstrap_results_complex <- bootstrap_lmer(
                data = df_tmp3,
                formula = CWM_trait ~ proxy_scaled + proxy_scaled_sq + (1|region),
                group_var = "region",
                num_bootstrap = num_bootstrap_samples
            )
            conf_intervals <- calculate_ci(bootstrap_results_complex)
            
            lm_results_complex$sq_estimate[n] <- sq_coefficient[1]
            lm_results_complex$sq_pvalue[n] <- sq_coefficient[2]
            lm_results_complex$sq_upper[n] <- conf_intervals[1, 3]
            lm_results_complex$sq_lower[n] <- conf_intervals[2, 3]
        } 
        
        # output complex model table:
        lm_results_complex$env_proxy[n] <- proxy
        lm_results_complex$bic_diff[n] <- bic_diff
        lm_results_complex$complexed[n] <- complexed
        
        lm_results_complex$estimate[n] <- coefficient[1]
        lm_results_complex$pvalue[n] <- coefficient[2]
        lm_results_complex$upper[n] <- conf_intervals[1, 2]
        lm_results_complex$lower[n] <- conf_intervals[2, 2]
        
    }
    
    return(list(lm_table_s = lm_results_simple,
                lm_table_c = lm_results_complex))
}


write_out_table <- function(trait, sp_analysis, type){
    
    trait_type <- ifelse(trait == "GS", "pool_GS_lmer_", "pool_PL_lmer_")
    analysis_type <- ifelse(type == "simple", "simple", "complex")
    output <- paste0("results/05_community_level/all_sp/", trait_type, analysis_type, ".txt")
    write.table(sp_analysis, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}

#####################################
sp_type <- "all_sp"
years <- 2008:2020
CWM_gs <- CWM(sp_type, years, "GS") 
CWM_pl <- CWM(sp_type, years, "PL")

community_gs_lmer <- community_lmer(CWM_gs, trait = "GS")
community_pl_lmer <- community_lmer(CWM_pl, trait = "PL")

write_out_table("GS", community_gs_lmer[[1]], 'simple')
write_out_table("GS", community_gs_lmer[[2]], 'complex')
write_out_table("PL", community_pl_lmer[[1]], 'simple')
write_out_table("PL", community_pl_lmer[[2]], 'complex')
