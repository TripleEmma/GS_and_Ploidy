## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: perform species level analysis on resident species dataset (incorporate quadratic terms)
# When the dependent variable is ploidy, a target environmental proxy (scaled) and its quadratic terms are independent variables . 
# When the dependent variable is genome size, ploidy level is included as a third independent variable. 
# When the dependent variable is genome size, we control for phylogenetic relationships. 

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
#     [1] geiger_2.0.11   boot_1.3-28.1   ggpubr_0.6.0    nortest_1.0-4   nlme_3.1-164   
# [6] caper_1.0.3     mvtnorm_1.2-4   MASS_7.3-60     phylolm_2.6.2   phytools_2.0-3 
# [11] maps_3.4.2      ape_5.7-1       lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1  
# [16] dplyr_1.1.4     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1   
# [21] ggplot2_3.4.4   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] subplex_1.8             Rcpp_1.0.11             lattice_0.22-5         
# [4] listenv_0.9.0           digest_0.6.33           foreach_1.5.2          
# [7] utf8_1.2.4              parallelly_1.36.0       R6_2.5.1               
# [10] backports_1.4.1         coda_0.19-4             pillar_1.9.0           
# [13] rlang_1.1.2             rstudioapi_0.15.0       car_3.1-2              
# [16] phangorn_2.11.1         Matrix_1.6-4            combinat_0.0-8         
# [19] igraph_1.6.0            munsell_0.5.0           broom_1.0.5            
# [22] compiler_4.2.0          numDeriv_2016.8-1.1     pkgconfig_2.0.3        
# [25] mnormt_2.1.1            optimParallel_1.0-2     globals_0.16.2         
# [28] tidyselect_1.2.0        expm_0.999-8            quadprog_1.5-8         
# [31] codetools_0.2-19        fansi_1.0.6             future_1.33.1          
# [34] tzdb_0.4.0              withr_2.5.2             grid_4.2.0             
# [37] gtable_0.3.4            lifecycle_1.0.4         magrittr_2.0.3         
# [40] scales_1.3.0            carData_3.0-5           future.apply_1.11.1    
# [43] cli_3.6.2               stringi_1.8.3           ggsignif_0.6.4         
# [46] scatterplot3d_0.3-44    doParallel_1.0.17       generics_0.1.3         
# [49] vctrs_0.6.5             fastmatch_1.1-4         deSolve_1.40           
# [52] iterators_1.0.14        tools_4.2.0             glue_1.6.2             
# [55] hms_1.1.3               abind_1.4-5             parallel_4.2.0         
# [58] timechange_0.2.0        colorspace_2.1-0        rstatix_0.7.2          
# [61] clusterGeneration_1.3.8


### species level analysis (quadratic regression)
library(tidyverse)
library(phytools)
library(phylolm)
library(caper)
library(ape)
library(nlme)
library(nortest)
library(ggpubr)
library(geiger)
library(lmtest)
library(boot)

rm(list = ls())

## phylogentic tree function
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

## assign each species an environmental preference (cover weighted)
# stress proxy; cover data
cover_data <- function(cover_type, proxy_file){
    # proxy_file <- "results/03_proxies/all_proxies_sp.txt"
    # sp_type <- "resident_sp"
    # cover_type <- paste0("cover_mean_not_zero_", sp_type, ".txt")
    sp_cover_proxies <- vector(mode = "list")
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies <- names(proxy_df)
    cover_dir <- "results/01_species_cover/"
    (cover_file <- list.files(path = cover_dir, pattern = cover_type))
    for (i in cover_file){
        # i <- "AEG_cover_mean_not_zero_resident_sp.txt"
        region <- unlist(strsplit(i, split = "_"))[1]
        region_proxy <- proxy_df %>% filter(str_detect(PlotID, region))
        (cover_region <- paste0(cover_dir, i))
        a <- read.delim(cover_region, stringsAsFactors = FALSE, check.names = FALSE)
        a <- a[match(region_proxy$PlotID, a$PlotID), ]
        b <- as.data.frame(apply(a[-1], 2, function(x){x/sum(x, na.rm = TRUE)})) # this gives species relative cover across the region
        # a[-1] removed the first column
        rownames(b) <- a$PlotID
        c <- vector(mode = "list")
        
        if(all(rownames(b) == region_proxy$PlotID)){
            for (j in 2:ncol(proxy_df)){
                # j = 2
                ind <- proxies[j]
                d <- colSums(apply(b, 2, function(x){x * region_proxy[, j]}), na.rm = TRUE)
                c[[ind]] <- as.data.frame(d) %>% rownames_to_column("species")
                names(c[[ind]]) <- c("species", ind)
            }
        } else{
            stop("The order of the plots is not correct!!!")
        }
        sp_cover_proxies[[region]] <- c %>% reduce(inner_join, by = "species")
    }
    return(sp_cover_proxies)
}

GS_analysis <- function(regions, region_color, proxy_file, trait_file, phylo_trees, sp_cover_proxies){
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    GS_regression <- data.frame(region = rep(NA, length(regions) * (length(proxies)-1)), 
                                region_sp_Num = NA, GS_sp_Num = NA, 
                                GS_PL_sp_Num = NA, env_proxy = NA,
                                resid_normality = NA, 
                                deviance_LRT = NA, pValue_LRT = NA,
                                complexed = NA,
                                estimate = NA, pvalue = NA, 
                                upper = NA, lower = NA,
                                sq_estimate = NA, sq_pvalue = NA, 
                                sq_upper = NA, sq_lower = NA)
    n <- 0
    norm_plots <- vector(mode = "list")
    Homoscedasticity_plots <- vector(mode = "list")
    for (i_ in seq_along(regions)){
        # i_ <- 2
        i <- regions[i_]
        norm_plots[[i]] <- vector(mode = "list") # each proxy would have normality plots
        Homoscedasticity_plots[[i]] <- vector(mode = "list") # same as normality plots
        
        region_df <- trait_df %>% filter(region == i) %>% 
            dplyr::select(region, species, GS, PL) %>% 
            filter(!is.na(GS)) %>% 
            mutate(log_GS = log(GS))
        region_df0 <- region_df
        rownames(region_df) <- region_df$species
        signal_tree <- phylo_trees(region_df$species)
        name_Ok <- name.check(signal_tree, region_df)
        if (name_Ok != 'OK'){stop("The species name is not matched")}
        region_df <- region_df[match(signal_tree$tip.label, region_df$species), ]
        signal_trait <- setNames(region_df$GS, region_df$species)
        
        region_df <- region_df %>% filter(!is.na(PL))
        cat("GS range: ", range(region_df$GS), "\n")
        # sp_cover_proxies_target <- sp_cover_proxies[[i]] %>% filter(species %in% region_df$species) # sp with trait available
        df_tmp <- region_df %>% inner_join(sp_cover_proxies[[i]], by = "species") 
        target_tree <- phylo_trees(df_tmp$species)
        target_tree <- unroot(target_tree)
        df_tmp <- df_tmp[match(target_tree$tip.label, df_tmp$species), ]
        
        for (j in seq(2, length(proxies))){
            n <- n + 1
            # j <- 2
            proxy <- proxies[j]
            df_tmp2 <- df_tmp %>% dplyr::select(species, region, all_of(proxy), GS, PL, log_GS)
            names(df_tmp2)[3] <- "proxy"
            df_tmp3 <- df_tmp2 %>% 
                mutate(proxy_scaled = scale(proxy),
                       proxy_sq_scaled = proxy_scaled^2) %>% 
                column_to_rownames("species")
            
            simple_mod <- phylolm(log_GS ~ proxy_scaled + PL, data = df_tmp3,
                                  phy = target_tree, model = "BM", boot = 10000)
            complex_mod <- phylolm(log_GS ~ proxy_scaled + proxy_sq_scaled + PL, data = df_tmp3,
                                   phy = target_tree, model = "BM", boot = 10000)
            
            deviance_LRT <- 2 * (complex_mod$logLik - simple_mod$logLik)
            pValue_LRT <- pchisq(deviance_LRT, df = 1, lower.tail = FALSE)
            complexed <- "N"
            plot_mod <- simple_mod
            target_mod <- complex_mod
            if(pValue_LRT < 0.05){
                complexed <- "Y"
                plot_mod <- complex_mod}
            
            transformed <- chol(solve(vcv(target_tree))) %*% target_mod$residuals
            normality <- lillie.test(transformed)
            mod_summary <- summary(target_mod)
            
            coefficient <- mod_summary$coefficients["proxy_scaled", 
                                                    c("Estimate", "p.value", 
                                                      "lowerbootCI", "upperbootCI")]
            sq_coefficient <- mod_summary$coefficients["proxy_sq_scaled", 
                                                       c("Estimate", "p.value", 
                                                         "lowerbootCI", "upperbootCI")]
            
            # output table:
            row_num <- n
            GS_regression$region[row_num] <- i
            GS_regression$region_sp_Num[row_num] <- nrow(trait_df %>% filter(region == i))
            GS_regression$GS_sp_Num[row_num] <- nrow(region_df0)
            GS_regression$GS_PL_sp_Num[row_num] <- nrow(df_tmp)
            GS_regression$env_proxy[row_num] <- proxy
            GS_regression$deviance_LRT[row_num] <- deviance_LRT
            GS_regression$pValue_LRT[row_num] <- pValue_LRT
            GS_regression$resid_normality[row_num] <- normality$p.value
            GS_regression$complexed[row_num] <- complexed
            GS_regression$estimate[row_num] <- coefficient[1]
            GS_regression$pvalue[row_num] <- coefficient[2]
            GS_regression$upper[row_num] <- coefficient[3]
            GS_regression$lower[row_num] <- coefficient[4]
            GS_regression$sq_estimate[row_num] <- sq_coefficient[1]
            GS_regression$sq_pvalue[row_num] <- sq_coefficient[2]
            GS_regression$sq_upper[row_num] <- sq_coefficient[3]
            GS_regression$sq_lower[row_num] <- sq_coefficient[4]
            
            # output plot:
            y_predicted = plot_mod$fitted.values
            transformed <- chol(solve(vcv(target_tree))) %*% plot_mod$residuals
            residual_df <- data.frame(Residual = transformed, predicted = y_predicted)
            
            plot_theme <- theme_classic() +
                theme(axis.text.x = element_text(size = 16),
                      axis.text.y = element_text(size = 16))
            
            norm_plot <- ggplot(residual_df, aes(sample = Residual)) +
                stat_qq(colour = region_color[i]) +
                stat_qq_line(colour = region_color[i]) +
                labs(x = NULL, y = NULL) +
                plot_theme
            norm_plots[[i]][[j-1]] <- norm_plot
            
            Homoscedasticity_plot <- ggplot(residual_df, aes(predicted, Residual)) +
                geom_point(colour = region_color[i], size = 2) +
                geom_hline(yintercept = 0, linetype = "dashed") +
                # ggtitle(paste0("Residual Homoscedasticity (", regions2[i], ')')) +
                labs(x = NULL, y = NULL) +
                plot_theme
            Homoscedasticity_plots[[i]][[j-1]] <- Homoscedasticity_plot
        }
    }
    return (list(GS_table = GS_regression,
                 GS_norm_plot = norm_plots,
                 GS_homocedasticity_plot = Homoscedasticity_plots))
}

PL_analysis <- function(regions, region_color, proxy_file, trait_file, phylo_trees, sp_cover_proxies){

    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    PL_regression <- data.frame(region = rep(NA, length(regions) * (length(proxies)-1)), 
                                region_sp_Num = NA, PL_sp_Num = NA, 
                                diploid_num = NA, polyploid_num = NA,
                                env_proxy = NA,
                                resid_normality = NA, 
                                deviance_LRT = NA, pValue_LRT = NA, 
                                complexed = NA,
                                estimate = NA, pvalue = NA, 
                                upper = NA, lower = NA,
                                sq_estimate = NA, sq_pvalue = NA, 
                                sq_upper = NA, sq_lower = NA)
    
    n <- 0
    norm_plots <- vector(mode = "list")
    Homoscedasticity_plots <- vector(mode = "list")
    
    fit_logistic_complex <- function(df_tmp3, indices, mod_type){
        resampled_data0 <- df_tmp3[indices, ]
        resampled_data <- resampled_data0 %>% 
            mutate(proxy_scaled = scale(proxy),
                   proxy_sq_scaled = proxy_scaled^2)
        model <- glm(ploidy ~ proxy_scaled + proxy_sq_scaled, 
                     family="binomial", data = resampled_data)
        return(coef(model))
    }
    
    for (i_ in seq_along(regions)){
        # i_ <- 3
        i <- regions[i_]
        norm_plots[[i]] <- vector(mode = "list")
        Homoscedasticity_plots[[i]] <- vector(mode = "list")
        region_df <- trait_df %>% filter(region == i) %>% 
            dplyr::select(region, species, PL, PL_source) %>% 
            filter(!is.na(PL)) %>% mutate(ploidy = if_else(PL > 2, 1, 0))
        cat("ploidy level range: ", range(region_df$PL), "\n")
        
        # sp_cover_proxies_target <- sp_cover_proxies[[i]] %>% filter(species %in% region_df$species)
        df_tmp <- region_df %>% inner_join(sp_cover_proxies[[i]], by = "species")

        for (j in seq(2, length(proxies))){
            n <- n + 1
            # j <- 10
            proxy <- proxies[j]
            df_tmp2 <- df_tmp %>% 
                dplyr::select(species, region, all_of(proxy), ploidy, PL)
            names(df_tmp2)[3] <- "proxy"
            row.names(df_tmp2) <- NULL
            df_tmp3 <- df_tmp2 %>% 
                mutate(proxy_scaled = scale(proxy),
                       proxy_sq_scaled = proxy_scaled^2) %>% 
                column_to_rownames("species")
            
            glm_mod_simple <- glm(ploidy ~ proxy_scaled,
                                  family="binomial", data = df_tmp3)
            glm_mod_complex <- glm(ploidy ~ proxy_scaled + proxy_sq_scaled,
                                   family="binomial", data = df_tmp3)

            # lr_test_result <- lrtest(glm_mod_complex, glm_mod_simple)
            # deviance_LRT <- lr_test_result[["Chisq"]][2]
            # pValue_LRT <- lr_test_result[["Pr(>Chisq)"]][2]
            deviance_LRT <- 2 * (logLik(glm_mod_complex) - logLik(glm_mod_simple))
            pValue_LRT <- pchisq(deviance_LRT, df = 1, lower.tail = FALSE)
            
            complexed <- 'N'
            plot_mod <- glm_mod_simple
            target_mod <- glm_mod_complex
            
            if(pValue_LRT < 0.05){
                complexed <- 'Y'
                plot_mod <- glm_mod_complex}
            
            mod_summary <- summary(target_mod)
            coefficient <- mod_summary$coefficients["proxy_scaled",
                                                    c("Estimate", "Pr(>|z|)")]
            sq_coefficient <- mod_summary$coefficients["proxy_sq_scaled",
                                                       c("Estimate", "Pr(>|z|)")]
            
            fit_logistic_model <- fit_logistic_complex
            num_bootstrap_samples <- 10000
            boot_results <- boot(data = df_tmp3, 
                                 statistic = fit_logistic_model, 
                                 R = num_bootstrap_samples)
            conf_intervals <- t(sapply(1:ncol(boot_results$t), function(i) {
                quantile(boot_results$t[, i], c(0.025, 0.975))
            }))
            upper <- conf_intervals[2, 1]
            lower <- conf_intervals[2, 2]
            sq_upper <- conf_intervals[3, 1]
            sq_lower <- conf_intervals[3, 2]
            
            transformed <- target_mod$residuals
            normality <- lillie.test(transformed)
            
            # output table:
            row_num <- n
            PL_regression$region[row_num] <- i
            PL_regression$region_sp_Num[row_num] <- nrow(trait_df %>% filter(region == i))
            PL_regression$PL_sp_Num[row_num] <- nrow(region_df)
            PL_regression$diploid_num[row_num] <- sum(region_df$ploidy == 0)
            PL_regression$polyploid_num[row_num] <- sum(region_df$ploidy == 1)
            PL_regression$env_proxy[row_num] <- proxy
            PL_regression$resid_normality[row_num] <- normality$p.value
            PL_regression$deviance_LRT[row_num] <- deviance_LRT
            PL_regression$pValue_LRT[row_num] <- pValue_LRT
            PL_regression$complexed[row_num] <- complexed
            PL_regression$estimate[row_num] <- coefficient[1]
            PL_regression$pvalue[row_num] <- coefficient[2]
            PL_regression$upper[row_num] <- upper
            PL_regression$lower[row_num] <- lower
            PL_regression$sq_estimate[row_num] <- sq_coefficient[1]
            PL_regression$sq_pvalue[row_num] <- sq_coefficient[2]
            PL_regression$sq_upper[row_num] <- sq_upper
            PL_regression$sq_lower[row_num] <- sq_lower
            
            # output plot:
            y_predicted <- plot_mod$fitted.values
            transformed <- plot_mod$residuals
            residual_df <- data.frame(Residual = plot_mod$residuals, 
                                      predicted = y_predicted,
                                      trans_resid = transformed)
            
            plot_theme <- theme_classic() +
                theme(axis.text.x = element_text(size = 16),
                      axis.text.y = element_text(size = 16))
            
            norm_plot <- ggplot(residual_df, aes(sample = trans_resid)) +
                stat_qq(colour = region_color[i]) +
                stat_qq_line(colour = region_color[i]) +
                # ggtitle(paste0("Residual QQ-plot (", regions2[i], ')')) +
                labs(x = NULL, y = NULL) +
                plot_theme
            norm_plots[[i]][[j-1]] <- norm_plot
            
            Homoscedasticity_plot <- ggplot(residual_df, aes(predicted, trans_resid)) +
                geom_point(colour = region_color[i], size = 2) +
                geom_hline(yintercept = 0, linetype = "dashed") +
                # ggtitle(paste0("Residual Homoscedasticity (", regions2[i], ')')) +
                theme_classic() +
                labs(x = NULL, y = NULL) +
                plot_theme
            Homoscedasticity_plots[[i]][[j-1]] <- Homoscedasticity_plot
        }
    }
    return(list(PL_table = PL_regression,
                PL_norm_plot = norm_plots,
                PL_homocedasticity_plot = Homoscedasticity_plots))
}

write_out_table <- function(sp_type, trait, sp_result){
    sp_table <- sp_result[[1]]
    analysis_type <- ifelse(trait == "GS", "GS_phylolm_", "PL_glm_")
    output <- paste0("results/04_species_level/", sp_type, "/", analysis_type, sp_type, "_quadratic.txt")
    # output <- paste0("results/04_species_level/", sp_type, "/", analysis_type, sp_type, "_quadratic_reverse.txt")
    write.table(sp_table, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}

save_plots <- function(sp_result, trait, sp_type, proxy_file){
    # sp_result <- sp_gs_result
    # sp_result <- sp_pl_result
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    
    norm_plots <- sp_result[[2]]
    Homoscedasticity_plots <- sp_result[[3]]
    
    for (i in seq(2, length(proxies))){
        j <- i -1
        norm_plots_all <- ggarrange(norm_plots$AEG[[j]],
                                    norm_plots$HEG[[j]],
                                    norm_plots$SEG[[j]],
                                    nrow = 3, ncol = 1)
        combinded_norm <- annotate_figure(norm_plots_all,
                                          bottom = text_grob("X", size = 15, x = 0.55, y = 0.9),
                                          left = text_grob("Y", size = 15, x = 1, y = 0.55, rot = 90),
                                          top = text_grob(paste0("[", trait, "] Residual QQ-plot (", proxies[j], ")"),
                                                          size = 16, x = 0.32, y = 0.6))
        
        
        Homoscedasticity_plots_all <- ggarrange(Homoscedasticity_plots$AEG[[j]],
                                                Homoscedasticity_plots$HEG[[j]],
                                                Homoscedasticity_plots$SEG[[j]],
                                                nrow = 3, ncol = 1)
        
        combinded_homoscedasticity <- annotate_figure(Homoscedasticity_plots_all,
                                                      bottom = text_grob("Predicted possiblity", size = 15, x = 0.55, y = 0.9),
                                                      left = text_grob("Residual", size = 15, x = 1, y = 0.55, rot = 90),
                                                      top = text_grob(paste0("[", trait, "] Residual Homoscedasticity (", proxies[j], ")"),
                                                                      size = 16, x = 0.39, y = 0.6))
        
        a <- ggarrange(combinded_norm, combinded_homoscedasticity)
        residual_output <- paste0("results/06_manuscript/", sp_type, "/residuals/",
                                  sp_type, "_", trait, "_", proxies[j],".png")
        png(residual_output,
            width = 32, height = 36, units = "cm", res = 600)
        print(a)
        dev.off()
    }
}

main_fun <- function(sp_type){
    set.seed(1909)
    # sp_type <- "resident_sp"
    # sp_type <- "all_sp"
    regions <- c("AEG", "HEG", "SEG")
    regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
    region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)
    proxy_file <- "results/03_proxies/all_proxies_sp.txt"
    cover_type <- paste0("cover_mean_not_zero_", sp_type, ".txt")
    trait_file <- paste0("results/02_traits/regions_traits_", sp_type, ".txt")
    sp_cover_proxies <- cover_data(cover_type, proxy_file)
    sp_gs_result <- GS_analysis(regions, region_color, proxy_file, trait_file, phylo_trees, sp_cover_proxies)
    sp_pl_result <- PL_analysis(regions, region_color, proxy_file, trait_file, phylo_trees, sp_cover_proxies)
    save_plots(sp_gs_result, "GS", sp_type, proxy_file)
    save_plots(sp_pl_result, "PL", sp_type, proxy_file)
    write_out_table(sp_type, "GS", sp_gs_result)
    write_out_table(sp_type, "PL", sp_pl_result)
    
}

# walk(c("resident_sp", "all_sp"), main_fun)
walk(c("resident_sp"), main_fun)
