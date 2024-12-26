## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: perform species level analysis on resident species dataset (phylogenetic linear regression models)
# different niche scores are responses and each genomic trait is predictor
# each region will have 20 different models (10 niche score * 2 genomic traits)

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

### species level analysis
library(tidyverse)
library(phylolm)
library(ape)
library(nlme)
library(geiger)
library(patchwork) 
library(ggpubr)

set.seed(1909)
rm(list = ls())

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
    # proxy_file <- "results/03_proxies/all_proxies_sp.txt" # (sp --> species level)
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
sp_type <- "resident_sp"
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)
proxy_file <- "results/03_proxies/all_proxies_sp.txt"
cover_type <- paste0("cover_mean_not_zero_", sp_type, ".txt")
trait_file <- paste0("results/02_traits/regions_traits_", sp_type, ".txt")
sp_cover_proxies <- cover_data(cover_type, proxy_file)

ploting <- function(fitted_values, transformed_residuals, region){
    region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), 
                             c("AEG", "HEG", "SEG"))
    
    results_df <- data.frame(
        FittedValues = fitted_values,
        TransformedResiduals = transformed_residuals
    )
    
    # Plot Fitted Values vs. Transformed Residuals
    p1 <- ggplot(results_df, aes(x = FittedValues, y = TransformedResiduals)) +
        geom_point(color = region_color[region]) +
        geom_hline(yintercept = 0, color = region_color[region], linetype = "dashed") +
        labs(x = "Fitted Values", y = "Transformed Residuals", 
             title = "Fitted Values vs. Transformed Residuals") +
        theme_minimal()
    
    # Create a QQ plot for the transformed residuals
    p2 <- ggplot(data.frame(Residuals = transformed_residuals), aes(sample = Residuals)) +
        stat_qq(color = region_color[region]) +
        stat_qq_line(color = region_color[region]) +
        labs(title = "QQ Plot of Transformed Residuals") +
        theme_minimal()
    p <- p1|p2
    return (p)
}
write_out_table <- function(sp_type, trait, sp_result){
    sp_table <- sp_result[[1]]
    analysis_type <- paste0(trait, "_phylolm_")
    output <- paste0("results/04_species_level/", sp_type, "/", analysis_type, 
                     sp_type, ".txt")
    write.table(sp_table, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}
save_plots <- function(sp_type, trait, sp_result, proxy_file){
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)[-1]
    
    plots <- sp_result[[2]]
    
    for (j in seq(1, length(proxies))){
        # j <- i - 1
        plots_all <- ggarrange(plots$AEG[[j]],
                               plots$HEG[[j]],
                               plots$SEG[[j]],
                               nrow = 3, ncol = 1)
        
        residual_output <- paste0("results/04_species_level/", sp_type, "/",
                                  sp_type, "_", trait, "_", proxies[j],".png")
        png(residual_output,
            width = 32, height = 36, units = "cm", res = 600)
        print(plots_all)
        dev.off()
    }
}

sp_GS_analysis <- function(regions, region_color, proxy_file, trait_file, 
                           phylo_trees, sp_cover_proxies){
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd",
                   "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    GS_regression <- data.frame(region = rep(NA, length(regions) * (length(proxies)-1)), 
                                region_sp_Num = NA, GS_sp_Num = NA, env_proxy = NA,
                                signalL = NA, signalL_p = NA,
                                signalK = NA, signalK_p = NA,
                                resi_signal = NA, resi_P = NA,
                                estimate = NA, pvalue = NA, 
                                lower = NA,upper = NA
    )
    
    n <- 0
    plots <- vector(mode = "list", length = 3)
    
    for (i_ in seq_along(regions)){
        # i_ <- 1
        i <- regions[i_]
        plots[[i]] <- vector(mode = "list") # each proxy would have normality plots
        
        region_df <- trait_df %>% filter(region == i) %>% 
            dplyr::select(region, species, GS) %>% 
            filter(!is.na(GS)) %>% 
            mutate(log_GS = log(GS)) %>% 
            mutate(logGS_scaled = scale(log_GS)) %>% 
            inner_join(sp_cover_proxies[[i]], by = "species") 
        
        tree_target <- phylo_trees(region_df$species)
        df_tmp <- region_df[match(tree_target$tip.label, region_df$species), ]
        rownames(df_tmp) <- df_tmp$species
        name_Ok <- name.check(tree_target, df_tmp)
        
        if (length(name_Ok) != 1){stop("The species name is not matched")}
        
        tree_target <- unroot(tree_target)
        cat("GS range: ", range(df_tmp$GS), "\n")
        
        for (j in seq(2, length(proxies))){ # because the first is plotID
            n <- n + 1
            # j <- 3
            proxy <- proxies[j]
            df_tmp2 <- df_tmp %>% dplyr::select(species, region, all_of(proxy), 
                                                GS, log_GS, logGS_scaled)
            names(df_tmp2)[3] <- "proxy"
            df_tmp2$proxy_scaled <- scale(df_tmp2$proxy)
            
            
            test_signal <- setNames(df_tmp2$proxy, df_tmp2$species)
            signal_lambda <- phylosig(tree_target, test_signal, 
                                      method = "lambda", test = TRUE)
            signalL <- signal_lambda$lambda
            signalL_pValue <- signal_lambda$P
            
            signal_K <- phylosig(tree_target, test_signal, 
                                 test = TRUE, nsim = 1000)
            signalK <- signal_K$K
            signalK_pValue <- signal_K$P
            
            fit_simple <- phylolm(proxy_scaled ~ logGS_scaled,
                                  phy = tree_target, model = "BM", boot = 1000,
                                  data = df_tmp2)
            fit_target <- fit_simple
            
            residuals_signalL <- phylosig(tree_target, 
                                          residuals(fit_target), 
                                          method = "lambda", test = TRUE)
            # residuals_signalK <- phylosig(tree_target, 
            #                              residuals(fit_target), 
            #                              method = "K", test = TRUE)
            
            transformed_residuals <- chol(solve(vcv(tree_target))) %*% fit_target$residuals
            fitted_values <- fitted(fit_target)
            p <- ploting(fitted_values, transformed_residuals, i)
            
            plots[[i]][[j-1]] <- p
            
            # get statistics 
            mod_summary <- summary(fit_target)
            
            coefficient <- matrix(mod_summary$coefficients["logGS_scaled", 
                                                           c("Estimate", "p.value", 
                                                             "lowerbootCI", "upperbootCI")], 
                                  nrow = 1)
            
            
            # output table:
            row_num <- n
            GS_regression$region[row_num] <- i
            GS_regression$region_sp_Num[row_num] <- nrow(trait_df %>% filter(region == i))
            GS_regression$GS_sp_Num[row_num] <- nrow(region_df)
            GS_regression$env_proxy[row_num] <- proxy
            
            GS_regression$signalL[row_num]<- signalL
            GS_regression$signalL_p[row_num] <- signalL_pValue
            GS_regression$signalK[row_num]<- signalK
            GS_regression$signalK_p[row_num] <- signalK_pValue
            
            GS_regression$resi_signal[row_num] <- residuals_signalL$lambda
            GS_regression$resi_P[row_num] <- residuals_signalL$P
            
            GS_regression$estimate[row_num] <- coefficient[1, 1]
            GS_regression$pvalue[row_num] <- coefficient[1, 2]
            GS_regression$lower[row_num] <- coefficient[1, 3]
            GS_regression$upper[row_num] <- coefficient[1, 4]
            
        }
        
    }
    return (list(GS_table = GS_regression,
                 GS_residual_plot = plots))
}

sp_PL_analysis <- function(regions, region_color, proxy_file, trait_file, 
                           phylo_trees, sp_cover_proxies){
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd",
                   "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    PL_regression <- data.frame(region = rep(NA, length(regions) * (length(proxies)-1)), 
                                region_sp_Num = NA, PL_sp_Num = NA, 
                                env_proxy = NA, 
                                signalL = NA, signalL_p = NA,
                                signalK = NA, signalK_p = NA,
                                resi_signal = NA, resi_P = NA,
                                estimate = NA, pvalue = NA, 
                                lower = NA, upper = NA)
    
    n <- 0
    plots <- vector(mode = "list", length = 3)
    
    for (i_ in seq_along(regions)){
        # i_ <- 2
        i <- regions[i_]
        plots[[i]] <- vector(mode = "list") # each proxy would have normality plots
        
        region_df <- trait_df %>% filter(region == i) %>% 
            dplyr::select(region, species, PL) %>% 
            filter(!is.na(PL)) %>% 
            inner_join(sp_cover_proxies[[i]], by = "species") 
        
        tree_target <- phylo_trees(region_df$species)
        df_tmp <- region_df[match(tree_target$tip.label, region_df$species), ]
        rownames(df_tmp) <- df_tmp$species
        name_Ok <- name.check(tree_target, df_tmp)
        
        if (length(name_Ok) != 1){stop("The species name is not matched")}
        
        tree_target <- unroot(tree_target)
        cat("PL range: ", range(df_tmp$PL), "\n")
        
        for (j in seq(2, length(proxies))){ # because the first is plotID
            n <- n + 1
            # j <- 8
            proxy <- proxies[j]
            df_tmp2 <- df_tmp %>% dplyr::select(species, region, 
                                                all_of(proxy), PL)
            names(df_tmp2)[3] <- "proxy"
            df_tmp2$proxy_scaled <- scale(df_tmp2$proxy)
            df_tmp2$PL_scaled <- scale(df_tmp2$PL)
            
            test_signal <- setNames(df_tmp2$proxy, df_tmp2$species)
            signal_lambda <- phylosig(tree_target, test_signal, 
                                      method = "lambda", test = TRUE)
            signalL <- signal_lambda$lambda
            signalL_pValue <- signal_lambda$P
            
            signal_K <- phylosig(tree_target, test_signal, 
                                 test = TRUE, nsim = 1000)
            signalK <- signal_K$K
            signalK_pValue <- signal_K$P
            
            fit_simple <- phylolm(proxy_scaled ~ PL_scaled,
                                  phy = tree_target, model = "BM", boot = 1000,
                                  data = df_tmp2)
            fit_target <- fit_simple
            
            # transform the residuals and plot them to see the assumptions
            residuals_signal <- phylosig(tree_target, 
                                         residuals(fit_target), 
                                         method = "lambda", test = TRUE)
            
            transformed_residuals <- chol(solve(vcv(tree_target))) %*% fit_target$residuals
            fitted_values <- fitted(fit_target)
            p <- ploting(fitted_values, transformed_residuals, i)
            plots[[i]][[j-1]] <- p
            
            # get statistics 
            mod_summary <- summary(fit_target)
            coefficient <- matrix(mod_summary$coefficients["PL_scaled", 
                                                           # coefficient <- matrix(mod_summary$coefficients["PL", 
                                                           c("Estimate", "p.value", 
                                                             "lowerbootCI", "upperbootCI")], 
                                  nrow = 1)
            
            # output table:
            row_num <- n
            PL_regression$region[row_num] <- i
            PL_regression$region_sp_Num[row_num] <- nrow(trait_df %>% filter(region == i))
            PL_regression$PL_sp_Num[row_num] <- nrow(region_df)
            PL_regression$env_proxy[row_num] <- proxy
            
            PL_regression$signalL[row_num]<- signalL
            PL_regression$signalL_p[row_num] <- signalL_pValue
            PL_regression$signalK[row_num]<- signalK
            PL_regression$signalK_p[row_num] <- signalK_pValue
            PL_regression$resi_signal[row_num] <- residuals_signal$lambda
            PL_regression$resi_P[row_num] <- residuals_signal$P
            
            PL_regression$estimate[row_num] <- coefficient[1, 1]
            PL_regression$pvalue[row_num] <- coefficient[1, 2]
            PL_regression$lower[row_num] <- coefficient[1, 3]
            PL_regression$upper[row_num] <- coefficient[1, 4]
            
        }
        
    }
    return (list(PL_table = PL_regression,
                 PL_residual_plot = plots))
}

sp_both_analysis <- function(regions, region_color, proxy_file, trait_file, 
                             phylo_trees, sp_cover_proxies){
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd",
                   "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    both_regression <- data.frame(region = rep(NA, length(regions) * (length(proxies)-1)), 
                                  region_sp_Num = NA, trait_sp_Num = NA, 
                                  env_proxy = NA,
                                  signalL = NA, signalL_p = NA,
                                  signalK = NA, signalK_p = NA,
                                  resi_signal = NA, resi_P = NA,
                                  
                                  PL_estimate = NA, PL_pvalue = NA, 
                                  PL_lower1 = NA, PL_upper1 = NA,
                                  PL_lower2 = NA, PL_upper2 = NA, 
                                  
                                  GS_estimate = NA, GS_pvalue = NA, 
                                  GS_lower1 = NA, GS_upper1 = NA,
                                  GS_lower2 = NA, GS_upper2 = NA)
    
    n <- 0
    
    for (i_ in seq_along(regions)){
        # i_ <- 2
        i <- regions[i_]
        
        region_df <- trait_df %>% filter(region == i) %>% 
            dplyr::select(region, species, PL, GS) %>% 
            filter(!is.na(PL), !is.na(GS)) %>% 
            inner_join(sp_cover_proxies[[i]], by = "species") 
        
        tree_target <- phylo_trees(region_df$species)
        df_tmp <- region_df[match(tree_target$tip.label, region_df$species), ]
        rownames(df_tmp) <- df_tmp$species
        name_Ok <- name.check(tree_target, df_tmp)
        
        if (length(name_Ok) != 1){stop("The species name is not matched")}
        
        tree_target <- unroot(tree_target)
        cat("PL range: ", range(df_tmp$PL), "\n")
        cat("GS range: ", range(df_tmp$GS), "\n")
        
        for (j in seq(2, length(proxies))){ # because the first is plotID
            n <- n + 1
            # j <- 8
            proxy <- proxies[j]
            df_tmp2 <- df_tmp %>% dplyr::select(species, region, all_of(proxy), 
                                                PL, GS)
            names(df_tmp2)[3] <- "proxy"
            df_tmp2$proxy_scaled <- scale(df_tmp2$proxy)
            df_tmp2$PL_scaled <- scale(df_tmp2$PL)
            df_tmp2$logGS_scaled <- scale(log(df_tmp2$GS))
            
            test_signal <- setNames(df_tmp2$proxy, df_tmp2$species)
            signal_lambda <- phylosig(tree_target, test_signal, 
                                      method = "lambda", test = TRUE)
            signalL <- signal_lambda$lambda
            signalL_pValue <- signal_lambda$P
            
            signal_K <- phylosig(tree_target, test_signal, 
                                 test = TRUE, nsim = 1000)
            signalK <- signal_K$K
            signalK_pValue <- signal_K$P
            
            fit_simple <- phylolm(proxy_scaled ~ PL_scaled + logGS_scaled,
                                  phy = tree_target, model = "BM", boot = 1000,
                                  data = df_tmp2)
            fit_target <- fit_simple
            
            # transform the residuals and plot them to see the assumptions
            residuals_signal <- phylosig(tree_target, 
                                         residuals(fit_target), 
                                         method = "lambda", test = TRUE)
            
            
            # get statistics 
            mod_summary <- summary(fit_target)
            coefficient <- matrix(mod_summary$coefficients[c("PL_scaled", "logGS_scaled"), 
                                                           # coefficient <- matrix(mod_summary$coefficients["PL", 
                                                           c("Estimate", "p.value", 
                                                             "lowerbootCI", "upperbootCI")], 
                                  nrow = 2)
            
            # output table:
            row_num <- n
            both_regression$region[row_num] <- i
            both_regression$region_sp_Num[row_num] <- nrow(trait_df %>% filter(region == i))
            both_regression$trait_sp_Num[row_num] <- nrow(region_df)
            both_regression$env_proxy[row_num] <- proxy
            
            both_regression$signalL[row_num]<- signalL
            both_regression$signalL_p[row_num] <- signalL_pValue
            both_regression$signalK[row_num]<- signalK
            both_regression$signalK_p[row_num] <- signalK_pValue
            both_regression$resi_signal[row_num] <- residuals_signal$lambda
            both_regression$resi_P[row_num] <- residuals_signal$P
            
            both_regression$PL_estimate[row_num] <- coefficient[1, 1]
            both_regression$PL_pvalue[row_num] <- coefficient[1, 2]
            both_regression$PL_lower1[row_num] <- coefficient[1, 3]
            both_regression$PL_upper1[row_num] <- coefficient[1, 4]
            
            
            both_regression$GS_estimate[row_num] <- coefficient[2, 1]
            both_regression$GS_pvalue[row_num] <- coefficient[2, 2]
            both_regression$GS_lower1[row_num] <- coefficient[2, 3]
            both_regression$GS_upper1[row_num] <- coefficient[2, 4]
            
            both_regression$PL_lower2[row_num] <- quantile(fit_target$bootstrap[,'PL_scaled'], 0.05)
            both_regression$PL_upper2[row_num] <- quantile(fit_target$bootstrap[,'PL_scaled'], 0.95)
            
            both_regression$GS_lower2[row_num] <- quantile(fit_target$bootstrap[,'logGS_scaled'], 0.05)
            both_regression$GS_upper2[row_num] <- quantile(fit_target$bootstrap[,'logGS_scaled'], 0.95)
        }
        
    }
    return (list(both_regression))
}


### final results output
sp_gs_result <- sp_GS_analysis(regions, region_color, proxy_file, trait_file, 
                               phylo_trees, sp_cover_proxies)
write_out_table(sp_type, 'GS', sp_gs_result)
save_plots(sp_type, 'GS', sp_gs_result, proxy_file)


sp_pl_result <- sp_PL_analysis(regions, region_color, proxy_file, trait_file, 
                               phylo_trees, sp_cover_proxies)
write_out_table(sp_type, 'PL', sp_pl_result )
save_plots(sp_type, 'PL', sp_pl_result, proxy_file)    


sp_both_result <- sp_both_analysis(regions, region_color, proxy_file, trait_file, 
                                   phylo_trees, sp_cover_proxies)

write_out_table(sp_type, 'both', sp_both_result)





