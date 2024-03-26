## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: produced Figure S7

## Contacts ----
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
#     [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] corrplot_0.92       gridExtra_2.3       gtable_0.3.4        data.table_1.14.10 
# [5] ggpubr_0.6.0        ggrepel_0.9.4       ggVennDiagram_1.4.9 scales_1.3.0       
# [9] lubridate_1.9.3     forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4        
# [13] purrr_1.0.2         readr_2.1.4         tidyr_1.3.0         tibble_3.2.1       
# [17] ggplot2_3.4.4       tidyverse_2.0.0    
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.11       pillar_1.9.0      compiler_4.2.0    tools_4.2.0      
# [5] lifecycle_1.0.4   timechange_0.2.0  pkgconfig_2.0.3   rlang_1.1.2      
# [9] cli_3.6.2         rstudioapi_0.15.0 withr_2.5.2       generics_0.1.3   
# [13] vctrs_0.6.5       hms_1.1.3         tidyselect_1.2.0  glue_1.6.2       
# [17] R6_2.5.1          rstatix_0.7.2     fansi_1.0.6       carData_3.0-5    
# [21] car_3.1-2         tzdb_0.4.0        magrittr_2.0.3    backports_1.4.1  
# [25] abind_1.4-5       colorspace_2.1-0  ggsignif_0.6.4    utf8_1.2.4       
# [29] stringi_1.8.3     munsell_0.5.0     broom_1.0.5      


library(tidyverse)
library(ggrepel)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)

########################### Figure 2 and 3 ###########################
# Forest plot for species level (Figure 2) and community level (Figure 3)
rm(list = ls())
forest_plot <- function(input_simple, input_quadratic, 
                        global_simple, global_complex,
                        trait = "PL", level = "sp"){
    # trait <- "GS"; level = "sp"
    # input_simple <- "results/04_species_level/resident_sp/GS_phylolm_resident_sp_simple.txt"
    # input_quadratic <- "results/04_species_level/resident_sp/GS_phylolm_resident_sp_quadratic.txt"
    # global_simple <- "justified/global_species/GS_phylolm_resident_sp_simple_global.txt"
    # global_complex <- "justified/global_species/GS_phylolm_resident_sp_complex_global.txt"
    
    # trait <- "PL"; level = "sp"
    # input_simple <- "results/04_species_level/resident_sp/PL_glm_resident_sp_simple.txt"
    # input_quadratic <- "results/04_species_level/resident_sp/PL_glm_resident_sp_quadratic.txt"
    # global_simple <- "justified/global_species/PL_glm_resident_sp_simple_global.txt"
    # global_complex <- "justified/global_species/PL_glm_resident_sp_complex_global.txt"
    
    
    # trait <- "GS"; level = "community"
    # input_simple <- "results/05_community_level/all_sp/simple_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt"
    # input_quadratic <- "results/05_community_level/all_sp/complex_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt"
    
    # trait <- "PL"; level = "community"
    # input_simple <- "results/05_community_level/all_sp/simple_PL_lm_all_sp_yearly_mean.txt"
    # input_quadratic <- "results/05_community_level/all_sp/complex_PL_lm_all_sp_yearly_mean.txt"
    df_global <- read_delim(global_simple)
    df_quadratic_global <- read_delim(global_complex)
    df_global$Quadratic <- df_quadratic_global$complexed
    df_global$region <- 'All'
    df_global <- df_global %>% dplyr::select(c(region, env_proxy, estimate, pvalue, upper, lower, Quadratic))
    
    df <- read_delim(input_simple)
    df_quadratic <- read_delim(input_quadratic)
    df$Quadratic <- df_quadratic$complexed
    df <- df %>% dplyr::select(c(region, env_proxy, estimate, pvalue, upper, lower, Quadratic))
    
    df <- rbind(df, df_global)
    df$region <- factor(df$region, levels = c('AEG', 'HEG', 'SEG', 'All'))
    previous_names <- c("Grazing", "Mowing", "Fertilization", "LUI", "Patch.Area", 
                        "C.patch.area.w.Gd", "C.patch.no.Gd", "C.edge.dens.Gd", 
                        "C.NN.dist.Gd", "C.Hanski.Gd")
    new_names <- c("Grazing", "Mowing", "Fertilization", "LUI", "Patch_size", 
                   "  Area_weighted\npatch_size  ", "Patch_number", "Edge_density",
                   "NN_distance", "Hanski\nconnectivity")
    env_names <- setNames(new_names, previous_names)
    y_order <- new_names
    
    df2 <- df %>% 
        filter(env_proxy %in% all_of(previous_names)) %>% 
        mutate(env_proxy = env_names[env_proxy]) %>%
        mutate(env_proxy = factor(env_proxy, levels = y_order)) %>% 
        arrange(rev(env_proxy), rev(region)) %>%
        mutate(index = c(1:nrow(.))) 
    
    insert_rows <- function(df, everyother, directions = "down"){
        insert_row <- data.frame(matrix(NA, nrow = nrow(df)/everyother, ncol = length(names(df))))
        names(insert_row) <- names(df)
        insert_row$index <- c(1:nrow(insert_row)) * everyother + 0.5
        df_new <- rbind(df, insert_row) %>% 
            arrange(index) %>% 
            tidyr::fill(env_proxy, .direction = directions) %>% 
            mutate(index = c(1:nrow(.)))
        return (df_new)
    }
    
    
    df2_1 <- insert_rows(df2, 4)
    df2_2 <- insert_rows(df2_1, 5)
    df2_3 <- insert_rows(df2_2, 6, directions = "up")
    df2_4 <- insert_rows(df2_3, 7, directions = "up")
    df2_5 <- df2_4 %>% 
        slice(1:(n() - 2)) %>% 
        add_row(.before = 1) %>% 
        add_row(.before = 1) %>% 
        tidyr::fill(env_proxy, .direction = "up") %>% 
        mutate(index = c(1:nrow(.)))
    df3 <- df2_5 %>% 
        mutate(sig = case_when(pvalue >= 0.05 ~ '',
                               pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                               pvalue < 0.01 & pvalue >= 0.001 ~ '**', 
                               pvalue < 0.001 ~ '***')) %>% 
        pivot_longer(c(upper, lower), values_to = "bounds", names_to = "xmin_xmax") %>% 
        dplyr::select(-xmin_xmax)
    
    add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.045
    add_range <-  0.2
    if (trait == "GS"){
        add_range <- 0.15
    }
    
    x_range <- max(abs(range(df3$bounds, na.rm = TRUE))) + add_range 
    
    strip_data <- df2_5 %>% dplyr::select(env_proxy) %>% 
        mutate(xmin = -x_range, xmax = x_range,
               y_position = rev(1:nrow(.)), 
               ymin = y_position - 0.5,
               ymax = y_position + 0.5,
               fill = rep(c("a", "b"), each = 8, length = nrow(.))) %>% 
        pivot_longer(cols = c(xmin, xmax), values_to = "x", names_to = "xmin_xmax") %>% 
        dplyr::select(-xmin_xmax)
    
    # trait <- "GS"
    # level <- "community"
    y_lab <- NULL
    Title <- "Ploidy Level"
    if (trait == "GS"){
        Title <- "Genome Size"
        y_lab <- "Environmental Proxies"
    }
    Level <- "(Species level)"
    if (level == "community"){
        Level <- "(Cover-weighted-Mean)"
    }
    Title <- paste(Title, Level)
    # y_text <- str_wrap(df2$env_proxy[seq(2, nrow(df2), 3)], width = 10)
    y_text <- str_wrap(df2$env_proxy[seq(2, nrow(df2), 4)], width = 1)
    
    if (trait == "GS" & level == "sp"){df3$sig_pos <- 0.04}
    if (trait == "GS" & level == "community"){df3$sig_pos <- 0.02}
    if (trait == "PL" & level == "sp"){df3$sig_pos <- 0.08}
    if (trait == "PL" & level == "community"){df3$sig_pos <- 0.08}
    
    plot_theme <- theme(axis.text.y = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        axis.title = element_text(size = 15),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.position = "bottom")
    
    if (trait == "PL") {
        plot_theme <- theme(axis.text.y = element_text(size = 12, margin = margin(r = 8, unit = "pt")),
                            axis.text.x = element_text(size = 12),
                            axis.title = element_text(size = 15),
                            plot.title = element_text(size = 16, face = "bold"),
                            axis.ticks.y = element_blank(),
                            legend.position = "bottom")
    }
    df4 <- df3 %>% filter(complete.cases(.))
    p <- ggplot(data = df3, aes(x = bounds, y = index, color = region, group = index)) +
        geom_ribbon(data = strip_data, 
                    aes(x = x, ymin = ymin, ymax = ymax, 
                        group = y_position, fill = fill),
                    inherit.aes = FALSE) + 
        scale_fill_manual(name = NULL, breaks = c("a", "b"),
                          values = c("#EFF0F0", "white")) +
        scale_x_continuous(#limits=c(-x_range, x_range),
            expand=c(0,0)) +
        scale_y_continuous(limits=c(0.5, nrow(df2_5)+0.5),
                           breaks=seq(4.5, nrow(df2_5), 8),
                           labels = y_text,
                           expand=c(0,0)) +
        geom_text(aes(x = add_text, y = index, label = round(estimate, 3)), 
                  size = 2.6, show.legend = FALSE) +
        # geom_text(aes(x = add_text + 0.01, 
        geom_text(aes(x = add_text + sig_pos, y = index - 0.5, label = sig),
                  hjust = 0, size = 6, , show.legend = FALSE) +
        geom_line(data = df4, mapping = aes(linetype =  Quadratic), linewidth = 1) +
        geom_point(size = 0.0001, show.legend = FALSE) +
        geom_point(aes(x = estimate)) + 
        labs(x = expression(paste("Standardized coefficient (", beta, ")")),
             y = y_lab) + 
        guides(fill = "none") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_color_manual(name = "Region",
                           # "#F8766D" "#7CAE00" "#00BFC4" "#C77CFF"
                           values = c("#f8766d", "#7cae00", "#00bfc4", "#C77CFF"),
                           labels = c("Alb", "Hai", "Sch", "All"),
                           breaks = c("AEG", "HEG", "SEG", "All")) +
        plot_theme +
        ggtitle(Title) 
    return(p)
}
##### sp
sp_inputfiles <- c("results/04_species_level/resident_sp/GS_phylolm_resident_sp_simple.txt",
                   "results/04_species_level/resident_sp/GS_phylolm_resident_sp_quadratic.txt",
                   "justified/global_species/GS_phylolm_resident_sp_simple_global.txt",
                   "justified/global_species/GS_phylolm_resident_sp_complex_global.txt",
                   "results/04_species_level/resident_sp/PL_glm_resident_sp_simple.txt",
                   "results/04_species_level/resident_sp/PL_glm_resident_sp_quadratic.txt",
                   "justified/global_species/PL_glm_resident_sp_simple_global.txt",
                   "justified/global_species/PL_glm_resident_sp_complex_global.txt")

inputfiles <- sp_inputfiles
level <- "sp"
sp_type <- "residenet_sp"
GS_plot <- forest_plot(inputfiles[1], inputfiles[2], 
                       inputfiles[3], inputfiles[4], 
                       trait = "GS", level)
output_dir <- "justified/global_species/"
output_gs <- paste0(output_dir, sp_type, "_GS_forest_plot.png")
ggsave(output_gs, width = 17, height = 22, units = "cm")
PL_plot <- forest_plot(inputfiles[5], inputfiles[6], 
                       inputfiles[7], inputfiles[8], 
                       trait = "PL", level)
output_pl <- paste0(output_dir, sp_type, "_PL_forest_plot.png")
ggsave(output_pl, width = 17, height = 22, units = "cm")


##### community_lmer 
community_inputfiles <- c("results/05_community_level/all_sp/simple_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt",
                          "results/05_community_level/all_sp/complex_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt",
                          "justified/global_community/simple_phylogeny_controled_GS_lmer_all_sp_yearly_mean_global.txt",
                          "justified/global_community/complex_phylogeny_controled_GS_lmer_all_sp_yearly_mean_global.txt",
                          "results/05_community_level/all_sp/simple_PL_lm_all_sp_yearly_mean.txt",
                          "results/05_community_level/all_sp/complex_PL_lm_all_sp_yearly_mean.txt",
                          "justified/global_community/simple_PL_lmer_all_sp_yearly_mean_global.txt",
                          "justified/global_community/complex_PL_lmer_all_sp_yearly_mean_global.txt")
inputfiles <- community_inputfiles
level <- "community"
sp_type <- "all_sp"
output_dir <- "justified/global_community/"
GS_plot <- forest_plot(inputfiles[1], inputfiles[2], 
                       inputfiles[3], inputfiles[4],
                       trait = "GS", level)
output_gs <- paste0(output_dir, sp_type, "_GS_forest_plot_lmer.png")
ggsave(output_gs, width = 17, height = 22, units = "cm")
PL_plot <- forest_plot(inputfiles[5], inputfiles[6], 
                       inputfiles[7], inputfiles[8],
                       trait = "PL", level)
output_pl <- paste0(output_dir, sp_type, "_PL_forest_plot_lmer.png")
ggsave(output_pl, width = 17, height = 22, units = "cm")




