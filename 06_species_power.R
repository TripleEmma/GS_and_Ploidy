## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: species level power analysis

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


### species level power analysis
library(tidyverse)
library(phylolm)
library(ape)
library(nlme)
library(geiger)

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
sp_type <- "resident_sp"
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)
proxy_file <- "results/03_proxies/all_proxies_sp.txt"
cover_type <- paste0("cover_mean_not_zero_", sp_type, ".txt")
trait_file <- paste0("results/02_traits/regions_traits_", sp_type, ".txt")
sp_cover_proxies <- cover_data(cover_type, proxy_file)

sp_analysis <- function(regions, region_color, proxy_file, trait_file, 
                        phylo_trees, sp_cover_proxies){
    proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
    proxies_all <- names(proxy_df)
    to_remove <- c("Grazing_sd", "Mowing_sd", "Fertilization_sd", "LUI_sd", 
                   "C.PX.100.Gd", "C.shape.ind.Gd")
    proxies <- setdiff(proxies_all, to_remove)
    trait_df <- read.delim(trait_file, stringsAsFactors = FALSE, sep = "\t")
    
    GS_power <- data.frame(region = rep(NA, 
                                        length(regions) * (length(proxies)-1) * length(seq(0, 1, 0.05))), 
                           env_proxy = NA, r = NA, power = NA)
    
    PL_power <- data.frame(region = rep(NA, 
                                        length(regions) * (length(proxies)-1) * length(seq(0, 1, 0.05))), 
                           env_proxy = NA, r = NA, power = NA)
    
    n = 0
    for (i_ in seq_along(regions)){
        # i_ <- 2
        i <- regions[i_]
        
        ### GS
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
        
        ### PL 
        region_df_pl <- trait_df %>% filter(region == i) %>% 
            dplyr::select(region, species, PL) %>% 
            filter(!is.na(PL)) %>% 
            mutate(PL_scaled = scale(PL)) %>% 
            inner_join(sp_cover_proxies[[i]], by = "species") 
        
        tree_target_pl <- phylo_trees(region_df_pl$species)
        df_tmp_pl <- region_df_pl[match(tree_target_pl$tip.label, region_df_pl$species), ]
        rownames(df_tmp_pl) <- df_tmp_pl$species
        name_Ok_pl <- name.check(tree_target_pl, df_tmp_pl)
        
        if (length(name_Ok_pl) != 1){stop("The species name is not matched")}
        
        tree_target_pl <- unroot(tree_target_pl)
        cat("PL range: ", range(df_tmp_pl$PL), "\n")
        
        
        for (j in seq(2, length(proxies))){ # because the first is plotID
            n <- n + 1
            # j <- 2
            proxy <- proxies[j]
            
            ###GS
            df_tmp2 <- df_tmp %>% dplyr::select(species, region, all_of(proxy),
                                                GS, log_GS, logGS_scaled)
            names(df_tmp2)[3] <- "proxy"
            df_tmp2$proxy_scaled <- scale(df_tmp2$proxy)
            
            target_mod <- phylolm(proxy_scaled ~ 1, data = df_tmp2,
                                  phy = tree_target, model = "BM")
            
            transformed <- chol(solve(vcv(tree_target))) %*% target_mod$residuals
            y_predicted = target_mod$fitted.values
            var_res = var(transformed)
            
            ### PL
            df_tmp2_pl <- df_tmp_pl %>% dplyr::select(species, region, all_of(proxy), 
                                                      PL, PL_scaled)
            names(df_tmp2_pl)[3] <- "proxy"
            df_tmp2_pl$proxy_scaled <- scale(df_tmp2_pl$proxy)
            
            target_mod_pl <- phylolm(proxy_scaled ~ 1, data = df_tmp2_pl,
                                     phy = tree_target_pl, model = "BM")
            
            transformed_pl <- chol(solve(vcv(tree_target_pl))) %*% target_mod_pl$residuals
            y_predicted_pl = target_mod_pl$fitted.values
            var_res_pl = var(transformed_pl)
            
            
            power_test <- function(r, df_tmp2, y_predicted, var_res, tree_target){
                # r <- 0.5
                sigma = matrix(c(1 * var_res, r * sqrt(var_res), r * sqrt(var_res), 1), ncol = 2)
                # ps <- replicate(10000, {
                ps <- replicate(5000, {
                    new <- MASS::mvrnorm(nrow(df_tmp2), mu = c(0, 0), Sigma = sigma)
                    
                    new_df <- data.frame(region = df_tmp2$region, feature = new[, 2], 
                                         proxy = y_predicted + new[, 1])
                    rownames(new_df) <- rownames(df_tmp2)
                    
                    new_mod <- phylolm(proxy~feature, data = new_df,
                                       phy = tree_target, model = "BM")
                    sig <- summary(new_mod)$coefficient["feature", "p.value"]
                })
                return (mean(ps<0.05))
            }
            
            r <- seq(0, 1, 0.05)
            row_num_start <- n * length(seq(0, 1, 0.05)) - (length(seq(0, 1, 0.05)) - 1)
            row_num <- n * length(seq(0, 1, 0.05))
            
            ### GS
            powers <- purrr::map(r, ~power_test(.x, df_tmp2, y_predicted,
                                                var_res, tree_target))
            GS_power$region[row_num_start: row_num] <- i
            GS_power$env_proxy[row_num_start: row_num] <- proxy
            GS_power$r[row_num_start: row_num] <- r
            GS_power$power[row_num_start: row_num] <- unlist(powers)
            
            ### PL
            powers_pl <- purrr::map(r, ~power_test(.x, df_tmp2_pl, y_predicted_pl, 
                                                   var_res_pl, tree_target_pl))
            PL_power$region[row_num_start: row_num] <- i
            PL_power$env_proxy[row_num_start: row_num] <- proxy
            PL_power$r[row_num_start: row_num] <- r
            PL_power$power[row_num_start: row_num] <- unlist(powers_pl)
        }
    }
    
    return (list(GS_power, PL_power))
}

sp_result <- sp_analysis(regions, region_color, proxy_file, trait_file, 
                         phylo_trees, sp_cover_proxies)

write_out_table <- function(sp_type, trait, sp_result){
    sp_table <- sp_result
    analysis_type <- ifelse(trait == "GS", "GS_power", "PL_power") # GS_power.txt
    output <- paste0("results/04_species_level/", sp_type, "/", analysis_type, ".txt")
    write.table(sp_table, sep = "\t",
                col.names = TRUE, row.names = FALSE,
                file = output)
}


### final results output
write_out_table(sp_type, 'GS', sp_result[[1]])
write_out_table(sp_type, 'PL', sp_result[[2]])

