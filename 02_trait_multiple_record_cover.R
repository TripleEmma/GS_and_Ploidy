#---------------------------------------------------------------------#
## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: address whether the species with multiple entries mainly species with higher cover values 
## (that could have large impact on the outcome of the community-level analyses) ?
# this will produce the Figure S2

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


library(tidyverse)

rm(list = ls())
df <- read_delim("results/02_traits/BE_species_traits_both.txt") %>% 
    dplyr::select(BEname, multiPloidy, multiGS) %>% 
    rename(species = BEname)
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
plots <- vector(mode = 'list', length=6)
n = 0
for (i in regions){
    # i <- 'AEG'
    region_file <- paste0("results/01_species_cover/", i, "_cover_mean_not_zero_all_sp.txt")
    df_region <- read_delim(region_file)
    df_region_relative <- as.data.frame(t(apply(df_region[-1], 1, function(x){x/sum(x, na.rm = TRUE)}))) %>% 
        pivot_longer(everything(), names_to = "species", values_to = "cover") %>% 
        filter(cover > 0) %>% 
        left_join(df) %>% 
        mutate(PL=case_when(multiPloidy == "Y" ~ 2,
                            multiPloidy == 'N' ~ 1, 
                            is.na(multiPloidy) ~0)) %>% 
        mutate(GS=case_when(multiGS == "Y" ~ 2,
                            multiGS == 'N' ~ 1, 
                            is.na(multiGS) ~0)) %>% 
        dplyr::select(species, cover, GS, PL)
    bins <- 80
    df_region_relative$bin <- cut(df_region_relative$cover, breaks = bins, labels = FALSE)
    
    n = n + 1
    counts_gs <- df_region_relative %>%
        group_by(bin, GS) %>%
        summarise(count = n())
    counts_gs$log_count <- log10(counts_gs$count)
    
    gs_title <- paste0(regions2[[i]], " -- Genome size")
    region_gs <- ggplot(counts_gs, aes(x = bin, y = log_count, fill = factor(GS))) +
        geom_col(position = "stack", color = "black", alpha = 0.5) +
        xlab('Cover proportion') +
        ylab('log10(# Plots)') +
        ggtitle(gs_title) +
        theme_minimal() +
        theme(
            legend.position = 'bottom',
            legend.title = element_blank(), 
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 17, face = "bold")) +
        scale_fill_manual(values = c("0" = "#F8766D", "1" = "#00BA38", "2" = "#619CFF"),
                          labels = c('No trait available', 
                                     'Unique trait record', 
                                     'Multiple trait records')) 
    plots[[n]] <- region_gs
    
    n = n + 1
    counts_pl <- df_region_relative %>%
        group_by(bin, PL) %>%
        summarise(count = n())
    counts_pl$log_count <- log10(counts_pl$count)
    
    pl_title <- paste0(regions2[[i]], " -- Ploidy level")
    region_pl <- ggplot(counts_pl, aes(x = bin, y = log_count, fill = factor(PL))) +
        geom_col(position = "stack", color = "black", alpha = 0.5) +
        xlab('Cover proportion') +
        ylab('log10(# Plots)') +
        ggtitle(pl_title) +
        theme_minimal() +
        theme(
            legend.position = 'bottom',
            legend.title = element_blank(), 
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 17, face = "bold")) +
        scale_fill_manual(values = c("0" = "#F8766D", "1" = "#00BA38", "2" = "#619CFF"),
                          labels = c('No trait available', 
                                     'Unique trait record', 
                                     'Multiple trait records')) 
    plots[[n]] <- region_pl
    
}
plots_together <- ggarrange(plots[[1]], plots[[2]],
                            plots[[3]], plots[[4]],
                            plots[[5]], plots[[6]],
                            nrow = 3, ncol = 2,
                            common.legend = TRUE,
                            legend = "bottom")
output <- paste0("results/06_manuscript/FigureS2/FigureS2", "_relative_cover_trait.png")
png(output,
    width = 26, height = 30, units = "cm", res = 600)
print(plots_together)
dev.off()

