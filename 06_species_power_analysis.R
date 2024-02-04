## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: perform species level power analysis 
# For genome size, we incorporate phylogenetic signals

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
library(phytools)
library(phylolm)
library(ape)
library(nlme)
library(nortest)
library(ggpubr)
library(geiger)

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

GS_analysis <- function(regions, region_color, trait_file, phylo_trees){
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    powers <- vector(mode = "list", length = length(regions))
    
    for (i_ in seq_along(regions)){
        # i_ <- 1
        i <- regions[i_]
        region_df <- trait_df %>% filter(region == i) %>% 
            dplyr::select(region, species, GS, PL) %>% 
            filter(!is.na(GS)) %>% 
            mutate(log_GS = log(GS))
        rownames(region_df) <- region_df$species
        signal_tree <- phylo_trees(region_df$species)
        name_Ok <- name.check(signal_tree, region_df)
        if (name_Ok != 'OK'){stop("The species name is not matched")}
        region_df <- region_df[match(signal_tree$tip.label, region_df$species), ]
        signal_trait <- setNames(region_df$GS, region_df$species)
        
        region_df <- region_df %>% filter(!is.na(PL))
        cat("GS range: ", range(region_df$GS), "\n")
        df_tmp <- region_df
        target_tree <- phylo_trees(df_tmp$species)
        target_tree <- unroot(target_tree)
        
        
        df_tmp2 <- df_tmp %>% dplyr::select(species, region, GS, PL, log_GS)
        df_tmp3 <- df_tmp2[match(target_tree$tip.label, df_tmp2$species), ]

        simple_mod <- phylolm(log_GS ~ PL, data = df_tmp3,
                              phy = target_tree, model = "BM")
        target_mod <- simple_mod
        
        transformed <- chol(solve(vcv(target_tree))) %*% target_mod$residuals
        # transformed <- target_mod$residuals
        y_predicted = target_mod$fitted.values
        
        var_res = var(transformed)
        
        # var_res = var(target_mod$residuals)
        power_test <- function(r, df_tmp3, y_predicted, var_res){
            sigma = matrix(c(1 * var_res, r * sqrt(var_res), r * sqrt(var_res), 1), ncol = 2)
            ps <- replicate(10000, {
                new <- MASS::mvrnorm(nrow(df_tmp3), mu = c(0, 0), Sigma = sigma)
                new_df <- data.frame(region = df_tmp3$region, proxy_scaled = new[, 2], 
                                     log_GS = y_predicted + new[, 1], PL = df_tmp3$PL)
                rownames(new_df) <- rownames(df_tmp3)
                new_mod <- phylolm(log_GS ~ proxy_scaled + PL , data = new_df, 
                                   phy = target_tree, model = "BM")
                sig <- summary(new_mod)$coefficient["proxy_scaled", "p.value"]
            })
            return (mean(ps<0.05))
        }
        r <- seq(0.1, 1, 0.05)
        powers[[i]] <- purrr::map(r, ~power_test(.x, df_tmp3, y_predicted, var_res))
    }
    return(powers)
}

PL_analysis <- function(regions, region_color, trait_file, phylo_trees){
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    powers <- vector(mode = "list", length = length(regions))
    
    for (i_ in seq_along(regions)){
        # i_ <- 1
        i <- regions[i_]
        cat(i, "\n")
        region_df <- trait_df %>% filter(region == i) %>% 
            dplyr::select(region, species, PL, PL_source) %>% 
            filter(!is.na(PL)) %>% mutate(ploidy = if_else(PL > 2, 1, 0))
        
        df_tmp2 <- region_df %>% dplyr::select(species, region, ploidy, PL)
        row.names(df_tmp2) <- NULL
        df_tmp3 <- df_tmp2 %>% 
            column_to_rownames("species")
        
        simple_mod <- glm(ploidy ~ 1, family="binomial", data = df_tmp3)
        target_mod <- simple_mod
        
        transformed <- target_mod$residuals
        y_predicted <- target_mod$fitted.values
        
        var_res = var(transformed)
        
        # var_res = var(target_mod$residuals)
        # cat("start power analysis\n")
        power_test <- function(r, df_tmp3, y_predicted, var_res){
            sigma = matrix(c(1 * var_res, r * sqrt(var_res), r * sqrt(var_res), 1), ncol = 2)
            ps <- replicate(10000, {
                new <- MASS::mvrnorm(nrow(df_tmp3), mu = c(0, 0), Sigma = sigma)
                log_scale <- log(y_predicted / (1-y_predicted)) + new[, 1]
                inv_log_scale <- arm::invlogit(log_scale)
                new_PL <- rbinom(nrow(df_tmp3), 1, inv_log_scale)
                new_df <- data.frame(region = df_tmp3$region,  proxy_scaled = new[, 2], 
                                     ploidy = new_PL)
                rownames(new_df) <- rownames(df_tmp3)

                new_mod <- glm(ploidy ~ proxy_scaled, family="binomial", data = new_df)
                sig <- summary(new_mod)$coefficient["proxy_scaled", "Pr(>|z|)"]
            })
            return (mean(ps<0.05))
        }
        r <- seq(0.1, 1, 0.05)
        powers[[i]] <- purrr::map(r, ~power_test(.x, df_tmp3, y_predicted, var_res))
        
    }
    return(powers)
}

write_out_table <- function(trait, sp_type, power_result){
    analysis_type <- ifelse(trait == "GS", "GS_phylolm_", "PL_")
    r = seq(0.1, 1, 0.05)
    df <- data.frame(power = unlist(power_result),
                     Region = rep(regions, each = length(r)),
                     r = r)
    rownames(df) <- NULL
    
    output <- paste0("results/04_species_level/", analysis_type, sp_type, "_power_test.txt")
    write.table(df, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}

plot_powers <- function(trait, powers, regions, level = "sp"){
    
    r = seq(0.1, 1, 0.05)
    df <- data.frame(power = unlist(powers),
                     Region = rep(regions, each = length(r)),
                     r = r)
    
    plot_theme <- theme(# legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 16, face = "bold"),

        # legend.background = element_rect(fill = "transparent"),
        # legend.text = element_text(size = 10),
        # legend.title = element_blank(),
        legend.key.size = unit(0.65, 'cm'),
        legend.spacing = unit(0.5, "cm")
    )
    Title <- "Ploidy Level"
    if (trait == "GS"){
        Title <- "Genome Size"
    }
    Level <- "(Species level)"
    if (level == "community"){
        Level <- "(Cover-weighted-Mean)"
    }
    Title <- paste(Title, Level)
    
    p <- ggplot(df, aes(x = r, y = power, color = Region)) + 
        geom_line(linewidth = 1) +
        scale_color_manual(name = "Region", 
                           values = c("#f8766d", "#7cae00", "#00bfc4"),
                           labels = c("Alb", "Hai", "Sch"),
                           breaks = regions) +
        ggtitle(Title) +
        xlab("Correlation") + ylab("Power") + 
        theme_bw() +
        plot_theme
    
    return(p)
}

set.seed(1909)
sp_type <- "resident_sp"
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)
trait_file <- paste0("results/02_traits/regions_traits_", sp_type, ".txt")
sp_gs_result <- GS_analysis(regions, region_color, trait_file, phylo_trees)
sp_pl_result <- PL_analysis(regions, region_color, trait_file, phylo_trees)
write_out_table("GS", sp_type, sp_gs_result)
write_out_table("PL", sp_type, sp_pl_result)
gs_plot <- plot_powers("GS", sp_gs_result, regions)
pl_plot <- plot_powers("PL", sp_pl_result, regions)
output <- paste0("results/06_manuscript/Figure4/Figure4_species_level_power_test.png")
png(output, width = 24, height = 14, units = "cm", res = 600)
ggarrange(gs_plot, pl_plot, nrow = 1,
          common.legend = TRUE, legend="bottom")
dev.off()

# walk(c("all_sp", "resident_sp"), main_fun)
