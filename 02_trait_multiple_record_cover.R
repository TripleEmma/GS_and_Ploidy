#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: address whether the species with multiple entries mainly species with higher cover values 
## (that could have large impact on the outcome of the community-level analyses) ?
# this will produce the Figure S2

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
# [1] ggpubr_0.6.0    Taxonstand_2.4  pbapply_1.7-2   lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1  
# [7] dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.4.4  
# [13] tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] pillar_1.9.0      compiler_4.2.0    tools_4.2.0       bit_4.0.5         lifecycle_1.0.4  
# [6] gtable_0.3.4      timechange_0.3.0  pkgconfig_2.0.3   rlang_1.1.3       cli_3.6.2        
# [11] rstudioapi_0.15.0 parallel_4.2.0    gridExtra_2.3     withr_3.0.0       generics_0.1.3   
# [16] vctrs_0.6.5       hms_1.1.3         cowplot_1.1.3     bit64_4.0.5       grid_4.2.0       
# [21] tidyselect_1.2.0  glue_1.7.0        R6_2.5.1          rstatix_0.7.2     fansi_1.0.6      
# [26] vroom_1.6.5       carData_3.0-5     farver_2.1.1      car_3.1-2         tzdb_0.4.0       
# [31] magrittr_2.0.3    backports_1.4.1   scales_1.3.0      abind_1.4-5       colorspace_2.1-0 
# [36] ggsignif_0.6.4    labeling_0.4.3    utf8_1.2.4        stringi_1.8.3     munsell_0.5.0    
# [41] broom_1.0.5       crayon_1.5.2   


library(tidyverse)
library(ggpubr)

rm(list = ls())
df <- read_delim("results/02_traits/BE_species_traits_both.txt") %>% 
    dplyr::select(BEname, multiPloidy, multiGS) %>% 
    rename(species = BEname)
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hainich", "Schorfheide"), regions)
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
        ylab('log10(# Species Contributing to Cover Proportion)') +
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
        ylab('log10(# Species Contributing to Cover Proportion)') +
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

