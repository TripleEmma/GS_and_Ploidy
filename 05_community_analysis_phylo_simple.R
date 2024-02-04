## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: perform community level analysis on all species dataset (simple linear)
# When the dependent variable is cover-weighted-ploidy, a target environmental proxy is the only independent variable (scaled). 
# When the dependent variable is cover-weighted-ploidy, ploidy level is included as a second independent variable. 
# When the dependent variable is cover-weighted-ploidy, we control for phylogenetic relationships by using gls(). 

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
# [1] ggpubr_0.6.0    nortest_1.0-4   nlme_3.1-164    geiger_2.0.11   phytools_2.0-3 
# [6] maps_3.4.2      ape_5.7-1       lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1  
# [11] dplyr_1.1.4     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1   
# [16] ggplot2_3.4.4   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.11             subplex_1.8             mvtnorm_1.2-4          
# [4] lattice_0.22-5          digest_0.6.33           foreach_1.5.2          
# [7] utf8_1.2.4              R6_2.5.1                backports_1.4.1        
# [10] coda_0.19-4             pillar_1.9.0            rlang_1.1.2            
# [13] rstudioapi_0.15.0       car_3.1-2               phangorn_2.11.1        
# [16] Matrix_1.6-4            combinat_0.0-8          labeling_0.4.3         
# [19] bit_4.0.5               igraph_1.6.0            munsell_0.5.0          
# [22] broom_1.0.5             compiler_4.2.0          numDeriv_2016.8-1.1    
# [25] pkgconfig_2.0.3         mnormt_2.1.1            optimParallel_1.0-2    
# [28] tidyselect_1.2.0        expm_0.999-8            quadprog_1.5-8         
# [31] codetools_0.2-19        fansi_1.0.6             crayon_1.5.2           
# [34] tzdb_0.4.0              withr_2.5.2             MASS_7.3-60            
# [37] grid_4.2.0              gtable_0.3.4            lifecycle_1.0.4        
# [40] magrittr_2.0.3          scales_1.3.0            vroom_1.6.5            
# [43] carData_3.0-5           cli_3.6.2               stringi_1.8.3          
# [46] farver_2.1.1            ggsignif_0.6.4          scatterplot3d_0.3-44   
# [49] doParallel_1.0.17       generics_0.1.3          vctrs_0.6.5            
# [52] cowplot_1.1.2           fastmatch_1.1-4         deSolve_1.40           
# [55] iterators_1.0.14        tools_4.2.0             bit64_4.0.5            
# [58] glue_1.6.2              hms_1.1.3               abind_1.4-5            
# [61] parallel_4.2.0          timechange_0.2.0        colorspace_2.1-0       
# [64] rstatix_0.7.2           clusterGeneration_1.3.8

library(tidyverse)
library(geiger)
library(nlme)
library(nortest)
library(ggpubr)

rm(list = ls())
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)
proxy_file <- "results/03_proxies/all_proxies_community.txt"
proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
proxies_all <- names(proxy_df)[-c(1, 2)]
to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
proxies <- setdiff(proxies_all, to_remove)

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
                filter(region == Region, !is.na(GS), !is.na(PL))
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
        a <- all_sp_year %>% 
            filter(grepl(Region, PlotID)) %>% 
            filter(PlotID %ni% excluded) %>% 
            dplyr::select(-Year) %>%  
            column_to_rownames("PlotID") %>% 
            filter(rowSums(is.na(.)) != ncol(.))
        region_sp <- names(a)
        
        CWM_gs <- function(col, df, brown){
            df$Col <- col
            pgls <- gls(log_GS ~ 1, data = df, # when only intercept, it is the mean
                        correlation = brown, weights = ~Col)
            # predicted_GS <- mean(exp(pgls$fitted))
            predicted_GS <- mean(pgls$fitted)
            return(predicted_GS)
        }
        
        if (trait == "GS"){
            region_trait <- trait_df %>% 
                filter(region == Region, !is.na(GS), !is.na(PL)) %>% 
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
            d <- as.data.frame(d)
            if (all(rownames(d) == region_trait$Species)){
                e <- as.data.frame(apply(d, 2, 
                                         function(x){CWM_gs(col = x, df = region_trait, brown = brown_cor)}))
                e <- e %>% rownames_to_column("PlotID")
                names(e)[2] <- "Plot_log_GS"}
            else{stop("something wrong!!!")}
            
            # calculate the CWM_pl for these plots using the same species set in those plots
            # the same code as those within else{}
            region_trait <- trait_df %>%
                filter(region == Region, !is.na(GS), !is.na(PL)) %>%
                dplyr::select(-contains("GS"))
            names(region_trait)[3] <- "trait_value"
            region_trait <- region_trait[match(target_tree$tip.label, region_trait$species), ]
            
            f <- as.data.frame(apply(d, 2, 
                                     function(x){sum(x * region_trait$trait_value, na.rm = TRUE)})) %>% 
                rownames_to_column("PlotID")
            names(f)[2] <- "Plot_pl"
            
            # merge plot_gs and plot_pl together
            e <- e %>% left_join(f, by = "PlotID")
            
        } else {
            region_trait <- trait_df %>% 
                filter(region == Region, !is.na(PL)) %>% 
                dplyr::select(-contains("GS"))
            names(region_trait)[3] <- "trait_value"
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
    CWM_trait <- vector(mode = "list")
    # trait <- "GS"; years <- 2008:2020; sp_type <- "all_sp"
    trait_all <- purrr::map(years, pipeline, sp_type, trait)
    names(trait_all) <- paste0("Year", years)
    
    for (i in regions){
        # i <- "AEG"
        all_trait_region_list <- lapply(trait_all, function(x) x[[i]])
        map_dbl(all_trait_region_list, nrow) # each year some plots were excluded because more than 20% cover were from species with missing trait
        # There were three plots missing from each year
        all_trait_region <- Reduce(function(df1, df2) full_join(df1, df2, by = "PlotID"), all_trait_region_list)
        if (trait == "GS"){
            headers <- purrr::map(names(trait_all), 
                                  function(x){paste0(x, c("_GS", "_PL"))}) %>% unlist()
            names(all_trait_region) <- c("PlotID", headers)
            last_non_na_column = pmap_chr(all_trait_region, get_last_non_na_column) %>% str_sub(end = 8)
            CWM_trait[[i]] <- all_trait_region %>% 
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
            names(all_trait_region) <- c("PlotID", names(trait_all))
            last_non_na_column = pmap_chr(all_trait_region, get_last_non_na_column)
            CWM_trait[[i]] <- all_trait_region %>% 
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
        
        
    }
    return (CWM_trait)
}

# plot the CWM-gs and CWM-pl value distributions in each region
community_trait_plots <- function(sp_type, cover_weighted_GS, cover_weighted_PL, region_color){
    CWM_gs <- bind_rows(cover_weighted_GS, .id = "Region")
    CWM_pl <- bind_rows(cover_weighted_PL, .id = "Region")
    names(CWM_gs)[3] <- "plot_trait"
    names(CWM_pl)[3] <- "plot_trait"
    print(paste("GS:", range(CWM_gs$plot_trait)))
    print(paste("GS mean:", mean(CWM_gs$plot_trait)))
    print(paste("PL:", range(CWM_pl$plot_trait)))
    print(paste("PL mean:", mean(CWM_pl$plot_trait)))
    set.seed(1909)
    plots <- function(CWM_trait, trait){
        ylabel <- "Cover weighted ploidy level"
        text_pos <- 6
        if (trait == "GS"){
            ylabel <- "Cover weighted genome size (log-transformed)"
            text_pos <- max(CWM_trait$plot_trait) + 0.1
        }
        
        plots_num <- CWM_trait %>% count(Region)
        ggplot(CWM_trait, aes(Region, plot_trait, color = Region)) +
            geom_violin(size = 0.8) + geom_jitter(size = 1.4, alpha = 0.5) +
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
                     sp_type, "_community_trait_violin.png")
    png(output, width = 20, height = 12, units = "cm", res = 600)
    print(ggarrange(p1, p2))
    dev.off()
}

# module3: regression analysis (lm)
community_lm <- function(regions, proxy_df, cwm_df, trait, region_color){
    # cwm_df <- CWM_gs
    # norm_plots <- vector(mode = "list")
    # Homoscedasticity_plots <- vector(mode = "list")
    proxies_all <- names(proxy_df)[-c(1, 2)]
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    
    lm_results <- data.frame(region = rep(NA, length(regions) * length(proxies)), 
                             region_plot_Num = NA, env_proxy = NA, 
                             resid_normality = NA, 
                             estimate = NA, pvalue = NA, 
                             upper = NA, lower = NA)
    n = 0
    for (i_ in seq(length(regions))){
        # i_ <- 1
        i <- regions[i_]
        region_ind <- proxy_df %>% filter(grepl(i, PlotID))
        
        df_tmp <- merge(cwm_df[[i]], region_ind, by = c("PlotID", "endYEAR"))
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
                           proxy_scaled_sq = proxy_scaled^2)
                simple_mod <- lm(CWM_gs ~ proxy_scaled + CWM_pl, data = df_tmp3)
            }
            else{
                df_tmp2 <- df_tmp %>% dplyr::select(PlotID, CWM_pl, all_of(proxy))
                names(df_tmp2)[3] <- "proxy"  ###### check here
                df_tmp3 <- df_tmp2 %>% 
                    mutate(proxy_scaled = scale(proxy),
                           proxy_scaled_sq = proxy_scaled^2)
                simple_mod <- lm(CWM_pl ~ proxy_scaled, data = df_tmp3)
            }
            target_mod <- simple_mod
            mod_summary <- summary(target_mod)
            coefficient <- mod_summary$coefficients["proxy_scaled", c(1, 4)]
            confints <- confint(target_mod)["proxy_scaled", ]
            normality <- lillie.test(residuals(target_mod))
            
            # output table:
            lm_results$region[n] <- i
            lm_results$region_plot_Num[n] <- nrow(df_tmp3)
            lm_results$env_proxy[n] <- proxy
            lm_results$resid_normality[n] <- normality$p.value
            lm_results$estimate[n] <- coefficient[1]
            lm_results$pvalue[n] <- coefficient[2]
            lm_results$upper[n] <- confints[1]
            lm_results$lower[n] <- confints[2]
        }
    }
    return(lm_results)
}

write_out_table <- function(trait, sp_type, sp_analysis){
    analysis_type <- ifelse(trait == "GS", "simple_phylogeny_controled_GS_lm_", "simple_PL_lm_")
    output <- paste0("results/05_community_level/", sp_type, "/", analysis_type, sp_type, "_yearly_mean.txt")
    write.table(sp_analysis, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}

#####################################
sp_type <- "all_sp"
years <- 2008:2020
CWM_gs <- CWM(sp_type, years, "GS") 
CWM_pl <- CWM(sp_type, years, "PL")
community_trait_plots(sp_type, CWM_gs, CWM_pl, region_color)

community_gs_lm <- community_lm(regions, proxy_df, CWM_gs, trait = "GS", region_color)
community_pl_lm <- community_lm(regions, proxy_df, CWM_pl, trait = "PL", region_color)

write_out_table("GS", sp_type, community_gs_lm)
write_out_table("PL", sp_type, community_pl_lm)
