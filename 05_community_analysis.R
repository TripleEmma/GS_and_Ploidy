## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: perform community level analysis on all species dataset (simple linear)
# cover-weighted genomic traits are response and each environmental proxy is predictor
# cover-weighted genome sizes have been corrected for phylogenetic relatness.
# A plot is excluded from a specific year if 20% of its cover that year was attributed to species with missing trait value

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

### community level analysis
library(tidyverse)
library(geiger)
library(nlme)
library(nortest)
library(patchwork) 
library(ggpubr)
library(boot)
library(caper)

rm(list = ls())
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
proxy_file <- "results/03_proxies/all_proxies_community.txt"
proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
proxies_all <- names(proxy_df)[-c(1, 2)]
to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
proxies <- setdiff(proxies_all, to_remove)
set.seed(0305)

# module1: get a list of excluded plots (because of not enough cover) in each region;
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
            GS_corrected <- residuals(corrected) + coef(corrected)
            predicted_GS <- sum(GS_corrected * df$Col)
            return(predicted_GS)
        }
        
        if (trait == "GS"){
            region_trait <- trait_df %>% 
                filter(region == Region, !is.na(GS)) %>% 
                mutate(log_GS = log(GS))
            rownames(region_trait) <- region_trait$species
            target_tree <- phylo_trees(region_trait$species)
            # target_tree <- unroot(target_tree)
            
            name_Ok <- name.check(target_tree, region_trait)
            if (name_Ok != 'OK'){stop("The species name is not matched")}
            region_trait <- region_trait[match(target_tree$tip.label, region_trait$species), ]
            Species <- region_trait$species
            region_trait <- region_trait %>% rename("Species" = species)
            brown_cor <- corBrownian(value = 1, phy = target_tree, form = ~Species)
            
            b <- cover_df_to_weight %>% dplyr::select(all_of(Species))
            c <- apply(b, 1, function(x){x/sum(x, na.rm = TRUE)})
            d <- apply(c, 2, function(x){ifelse((is.na(x) | x < 0.000001), 0.000001, x)})
            d <- as.data.frame(d)
            if (all(rownames(d) == region_trait$Species)){
                e <- as.data.frame(apply(d, 2, 
                                         function(x){CWM_gs(col = x, df = region_trait, brown = brown_cor)}))
                # function(x){CWM_gs(col = x, df = region_trait, brown = target_tree)}))
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
    CWM_trait <- vector(mode = "list")
    # trait <- "GS"; years <- 2008:2020; sp_type <- "all_sp"
    trait_all <- purrr::map(years, pipeline, sp_type, trait)
    names(trait_all) <- paste0("Year", years)
    
    for (i in regions){
        # i <- "AEG"
        all_trait_region_list <- lapply(trait_all, function(x) x[[i]])
        purrr::map_dbl(all_trait_region_list, nrow) # each year some plots were excluded because more than 20% cover were from species with missing trait
        # There were three plots missing from each year
        all_trait_region <- Reduce(function(df1, df2) full_join(df1, df2, by = "PlotID"), all_trait_region_list)
        names(all_trait_region) <- c("PlotID", names(trait_all))
        last_non_na_column = pmap_chr(all_trait_region, get_last_non_na_column)
        
        CWM_trait[[i]] <- all_trait_region %>% 
            rowwise() %>% 
            mutate(CWM_trait = rowMeans(pick(where(is.numeric)), na.rm = TRUE),
                   n = rowSums(!is.na(pick(-c(PlotID, CWM_trait))))
            ) %>% 
            ungroup() %>% 
            mutate(endYEAR = last_non_na_column) %>% 
            filter(n >= 6) %>% 
            dplyr::select(PlotID, CWM_trait, n, endYEAR) %>% 
            arrange(PlotID)
    }
    return (CWM_trait)
}

# Figure1d
community_trait_plots <- function(sp_type, cover_weighted_GS, cover_weighted_PL, region_color){
    CWM_gs <- bind_rows(cover_weighted_GS, .id = "Region")
    CWM_pl <- bind_rows(cover_weighted_PL, .id = "Region")
    print(paste("GS:", range(CWM_gs$CWM_trait)))
    print(paste("GS mean:", mean(CWM_gs$CWM_trait)))
    print(paste("PL:", range(CWM_pl$CWM_trait)))
    print(paste("PL mean:", mean(CWM_pl$CWM_trait)))
    set.seed(1909)
    plots <- function(CWM_trait, trait){
        ylabel <- "Cover weighted ploidy level"
        text_pos <- 6
        if (trait == "GS"){
            ylabel <- "Cover weighted genome size (log-transformed)"
            text_pos <- max(CWM_trait$CWM_trait) + 0.1
        }
        
        plots_num <- CWM_trait %>% count(Region)
        ggplot(CWM_trait, aes(Region, CWM_trait, color = Region)) +
            geom_boxplot(size = 0.8) + geom_jitter(size = 1.4, alpha = 0.5) +
            geom_text(data = plots_num, aes(y = text_pos, label = n), size = 5, color = 'black') + 
            scale_color_manual(values = region_color) + 
            ylab(ylabel) + 
            xlab(NULL) + theme_bw() +
            scale_x_discrete(labels = c("Alb", "Hai", "Sch")) +
            theme(legend.position = "none",
                  axis.text = element_text(size = 11),
                  axis.title = element_text(size = 12),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(linetype = "dashed"),
                  plot.margin = unit(margin(10, 10, 10, 12), "cm"))
    }
    p1 <- plots(CWM_gs, "GS")
    p2 <- plots(CWM_pl, "PL")
    
    output <- paste0("results/06_manuscript/", sp_type, "/Figure1d_",
                     sp_type, "_community_trait_box.png")
    png(output, width = 20, height = 12, units = "cm", res = 600)
    print(ggarrange(p1, p2))
    dev.off()
}

ploting <- function(target_mod, region, region_color){
    # region <- i
    # region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), 
    #                          c("AEG", "HEG", "SEG"))
    fitted_values <- fitted(target_mod)
    residuals <- target_mod$residuals
    results_df <- data.frame(
        FittedValues = fitted_values,
        res = residuals
    )
    
    # Plot Fitted Values vs. Transformed Residuals
    p1 <- ggplot(results_df, aes(x = FittedValues, y = res)) +
        geom_point(color = region_color[region]) +
        geom_hline(yintercept = 0, color = region_color[region], linetype = "dashed") +
        labs(x = "Fitted Values", y = "Residuals", 
             title = "Fitted Values vs. Residuals") +
        theme_minimal()
    
    # Create a QQ plot for the transformed residuals
    p2 <- ggplot(data.frame(Residuals = residuals), aes(sample = Residuals)) +
        stat_qq(color = region_color[region]) +
        stat_qq_line(color = region_color[region]) +
        labs(title = "QQ Plot of Residuals") +
        theme_minimal()
    p <- p1|p2
    return (p)
}
save_plots <- function(sp_type, trait, community_plots, proxy_file){
    # community_result <- community_gs_lm
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)[-c(1, 2)]
    
    plots <- community_plots
    
    for (i in seq(1, length(proxies))){
        # j <- i - 1
        plots_all <- ggarrange(plots$AEG[[i]],
                               plots$HEG[[i]],
                               plots$SEG[[i]],
                               nrow = 3, ncol = 1)
        
        residual_output <- paste0("results/05_community_level/", sp_type, "/",
                                  trait, "_", proxies[i],".png")
        png(residual_output,
            width = 32, height = 36, units = "cm", res = 600)
        print(plots_all)
        dev.off()
    }
}

write_out_table <- function(sp_type, trait, community_table, type){
    # community_table <- community_result[[1]]
    trait_type <- ifelse(trait == "GS", "phylogeny_controled_GS_lm_", "PL_lm_")
    analysis_type <- ifelse(type == "simple", "simple", "complex")
    output <- paste0("results/05_community_level/", sp_type, "/", analysis_type, 
                     "_", trait_type, "yearly_mean.txt")
    write.table(community_table, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}
## module3: regression analysis (lm)
# we use bootstrap to get confidence interval; construct boot function

boot_choice <- function(complexed = 'N'){
    
    fit_logistic_model <- function(df_tmp3, indices){
        resampled_data0 <- df_tmp3[indices, ]
        model <- lm(CWM_trait ~ proxy_scaled, data = resampled_data0)
        return(coef(model))
    }
    
    if (complexed == "Y"){
        fit_logistic_model <- function(df_tmp3, indices){
            resampled_data0 <- df_tmp3[indices, ]
            model <- lm(CWM_trait ~ proxy_scaled + proxy_scaled_sq, data = resampled_data0)
            return(coef(model))
        }
    }
    
    return(fit_logistic_model)
}

community_lm <- function(regions, proxy_df, cwm_df, trait, region_color){
    # cwm_df <- CWM_gs
    all_proxies_boxplot <- vector(mode = "list")
    proxies_all <- names(proxy_df)[-c(1, 2)]
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    # proxies <- c("rawFertilization")
    lm_results_simple <- data.frame(region = rep(NA, length(regions) * length(proxies)), 
                                    region_plot_Num = NA, env_proxy = NA, 
                                    resid_normality = NA, 
                                    estimate = NA, pvalue = NA, 
                                    upper = NA, lower = NA)
    
    lm_results_complex <- data.frame(region = rep(NA, length(regions) * length(proxies)), 
                                     region_plot_Num = NA, env_proxy = NA, 
                                     resid_normality = NA, 
                                     bic_diff = NA, complexed = NA,
                                     estimate = NA, pvalue = NA, 
                                     upper = NA, lower = NA,
                                     sq_estimate = NA, sq_pvalue = NA, 
                                     sq_upper = NA, sq_lower = NA)
    n = 0
    plots <- vector(mode = "list")
    
    for (i_ in seq(length(regions))){
        # i_ <- 1
        i <- regions[i_]
        region_ind <- proxy_df %>% filter(grepl(i, PlotID))
        
        df_tmp <- merge(cwm_df[[i]], region_ind, by = c("PlotID", "endYEAR"))
        # names(df_tmp)[3] <- "cover_weighted_trait"  ###### check here
        all_proxies_boxplot[[i]] <- df_tmp
        print(paste(trait, " range:", range(df_tmp$CWM_trait)))
        
        for(j in seq(length(proxies))){
            n <- n + 1
            # j <- 5; trait <- "gs"
            proxy <- proxies[j]
            cat(i, proxy, "\n")
            
            df_tmp2 <- df_tmp %>% dplyr::select(PlotID, CWM_trait, all_of(proxy))
            names(df_tmp2)[3] <- "proxy"  ###### check here
            df_tmp3 <- df_tmp2 %>% 
                mutate(proxy_scaled = scale(proxy),
                       CWM_trait0 = CWM_trait,
                       CWM_trait = scale(CWM_trait0),
                       proxy_scaled_sq = proxy_scaled^2)
            simple_mod <- lm(CWM_trait ~ proxy_scaled, data = df_tmp3)
            complex_mod <- lm(CWM_trait ~ proxy_scaled + proxy_scaled_sq, data = df_tmp3)
            
            complexed <- "N"
            target_mod <- simple_mod
            
            normality <- lillie.test(residuals(target_mod))
            mod_summary <- summary(target_mod)
            coefficient <- mod_summary$coefficients["proxy_scaled", c(1, 4)]
            
            fit_logistic_model <-  boot_choice(complexed)
            num_bootstrap_samples <- 1000
            boot_results <- boot(data = df_tmp3, 
                                 statistic = fit_logistic_model, 
                                 R = num_bootstrap_samples)
            conf_intervals <- t(sapply(1:ncol(boot_results$t), function(i) {
                quantile(boot_results$t[, i], c(0.025, 0.975))
            }))
            
            # output simple model table:
            lm_results_simple$region[n] <- i
            lm_results_simple$region_plot_Num[n] <- nrow(df_tmp3)
            lm_results_simple$env_proxy[n] <- proxy
            lm_results_simple$resid_normality[n] <- normality$p.value
            
            lm_results_simple$estimate[n] <- coefficient[1]
            lm_results_simple$pvalue[n] <- coefficient[2]
            lm_results_simple$upper[n] <- conf_intervals[2, 1]
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
                coefficient <- mod_summary$coefficients["proxy_scaled", c(1, 4)]
                sq_coefficient <- mod_summary$coefficients["proxy_scaled_sq", c(1, 4)]
                fit_logistic_model <-  boot_choice(complexed)
                num_bootstrap_samples <- 1000
                boot_results <- boot(data = df_tmp3, 
                                     statistic = fit_logistic_model, 
                                     R = num_bootstrap_samples)
                conf_intervals <- t(sapply(1:ncol(boot_results$t), function(i) {
                    quantile(boot_results$t[, i], c(0.025, 0.975))
                }))
                
                lm_results_complex$sq_estimate[n] <- sq_coefficient[1]
                lm_results_complex$sq_pvalue[n] <- sq_coefficient[2]
                lm_results_complex$sq_upper[n] <- conf_intervals[3, 1]
                lm_results_complex$sq_lower[n] <- conf_intervals[3, 2]
            } 
            
            # output complex model table:
            lm_results_complex$region[n] <- i
            lm_results_complex$region_plot_Num[n] <- nrow(df_tmp3)
            lm_results_complex$env_proxy[n] <- proxy
            lm_results_complex$resid_normality[n] <- normality$p.value
            lm_results_complex$bic_diff[n] <- bic_diff
            lm_results_complex$complexed[n] <- complexed
            
            lm_results_complex$estimate[n] <- coefficient[1]
            lm_results_complex$pvalue[n] <- coefficient[2]
            lm_results_complex$upper[n] <- conf_intervals[2, 1]
            lm_results_complex$lower[n] <- conf_intervals[2, 2]
            
            p <- ploting(target_mod, i, region_color)
            plots[[i]][[j]] <- p
            
        }
    }
    
    return (list(lm_table_s = lm_results_simple,
                 lm_table_c = lm_results_complex,
                 lm_residual_plot = plots,
                 proxy_boxplot = all_proxies_boxplot))
}


community_proxy_plots <- function(proxy_boxplot, region_color, trait){
    # proxy_boxplot <- community_gs_lm$proxy_boxplot
    # region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)
    names(region_color) <- NULL
    new_label <- c("Grazing", "Mowing", "Fertilization", "LUI", 
                   "Patch_size", "Area_weighted_patch_size", "Patch_number", 
                   "Edge_density", "NN_distance", "Hanski_connectivity")
    
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", 
                   "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd", 
                   "PlotID", "endYEAR", "CWM_trait", "n")
    
    CWM_proxy <- bind_rows(proxy_boxplot, .id = "region") %>% 
        mutate(Region = case_when(str_detect(PlotID, "AEG") ~ "Alb", 
                                  str_detect(PlotID, "HEG") ~ "Hai",
                                  str_detect(PlotID, "SEG") ~ "Sch")) %>% 
        dplyr::select(-all_of(to_remove), -region) %>% 
        pivot_longer(-Region, names_to = 'env_proxies', values_to = 'intensity') 
    
    CWM_proxy$env_proxies <- as.character(CWM_proxy$env_proxies)
    a <- setNames(new_label, unique(CWM_proxy$env_proxies))
    
    b <- CWM_proxy %>% mutate(new_label = a[CWM_proxy$env_proxies])
    b$new_label <- factor(b$new_label, levels = new_label)
    b$region <- factor(b$Region, levels = unique(b$Region))
    
    p <- ggplot(b, aes(x=Region, y=intensity, color = Region)) +
        geom_boxplot(size = 0.8) + 
        geom_jitter(size = 1.4, alpha = 0.5) +
        theme_bw() + xlab("") + ylab("") +
        facet_wrap(~new_label, scales = 'free_y') +
        scale_color_manual(values = region_color) + 
        theme(legend.position = "none",
              axis.text.x = element_text(size = 18),
              axis.text.y = element_text(size = 15),
              strip.text = element_text(size = 16),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(linetype = "dashed"))
    output <- paste0("results/06_manuscript/FigureS10/community_proxies_", trait, ".png")
    png(output, width = 36, height = 32, units = "cm", res = 600)
    print(p)
    dev.off()
}

#####################################
sp_type <- "all_sp"
years <- 2008:2020
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), 
                         c("AEG", "HEG", "SEG"))
CWM_gs <- CWM(sp_type, years, "GS") 
community_gs_lm <- community_lm(regions, proxy_df, CWM_gs, trait = "GS", region_color)
write_out_table(sp_type, "GS", community_gs_lm[[1]], "simple")
write_out_table(sp_type, "GS", community_gs_lm[[2]], "complex")
save_plots(sp_type, "GS", community_gs_lm[[3]], proxy_file)
community_proxy_plots(community_gs_lm[[4]], region_color, "GS")


CWM_pl <- CWM(sp_type, years, "PL")
community_pl_lm <- community_lm(regions, proxy_df, CWM_pl, trait = "PL", region_color)
write_out_table(sp_type, "PL", community_pl_lm[[1]], "simple")
write_out_table(sp_type, "PL", community_pl_lm[[2]], "complex")
save_plots(sp_type, "PL", community_pl_lm[[3]], proxy_file)
community_proxy_plots(community_pl_lm[[4]], region_color, "PL")

community_trait_plots(sp_type, CWM_gs, CWM_pl, region_color)

# write.table(bind_rows(CWM_gs), 'communityGS.txt')
# a <- bind_rows(CWM_gs) %>% 
#     mutate(region = case_when(str_detect(PlotID, 'AEG')~'AEG',
#                               str_detect(PlotID, 'HEG')~'HEG',
#                               str_detect(PlotID, 'SEG')~'SEG'))
# 
# ggplot(a, aes(x = CWM_pl, y = CWM_gs, color = region)) + geom_line()
