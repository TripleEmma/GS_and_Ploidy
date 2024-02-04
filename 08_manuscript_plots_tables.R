## Title ----
# Genome size and ploidy do not predict plant responses to land-use intensity and habitat fragmentation in temperate grassland

## Purposes: produced plots and tables

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
library("scales")
library("ggVennDiagram")
library(ggrepel)
library(ggpubr)
library(data.table)
library(gtable)
library(grid)
library(gridExtra)

regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)

########################### Figure S1 ###########################
# this figure is produced by script "02_traits_mean_median_diff.R"

########################### Figure S2 ###########################
# this figure is produced by script "02_multiple_trait_cover.R"

########################### Figure S3 ###########################
# FigureS3: correlation plot indicating landscape metrics are highly correlated with others
library(corrplot)
######### EP
rm(list = ls())
output_dir <- "results/06_manuscript/FigureS3/"
EP_path <- "data/BE/Landscape/EP/"
(EP_pattern <- list.files(path = EP_path, pattern = "_2008_LUagg2_EP\\.txt"))
list_EP <- vector(mode = "list")
for (i in EP_pattern){
    # i <- "Alb_2008_LUagg2_EP.txt"
    print(i)
    EP_file <- paste0(EP_path, i)
    region <- unlist(strsplit(i, '_'))[1]
    year <- unlist(strsplit(i, '_'))[2]
    resolution <- unlist(strsplit(i, '_'))[3]
    title <- paste0(region, "-", year)
    list_name <- paste0(region, "_", year, "_", resolution)
    tmp_df <- read.delim(EP_file, sep = " ", 
                         na.strings = c("NA", "", "Inf"),
                         stringsAsFactors = FALSE)
    print(unique(tmp_df$LUagg2))
    
    output_file <- paste0(output_dir, "EP_", list_name, ".png")
    png(output_file, width = 12, height = 8, units = "cm", res = 600)
    M <- cor(tmp_df[, 4:length(tmp_df)], use = "complete.obs")
    M.sig = cor.mtest(M)
    corrplot(M, method = 'number', type = "lower",
             tl.cex = 0.3, tl.col = "black",
             number.cex = 0.3, cl.cex = 0.3,
             mar = c(0, 0, 0, 0)
    )
    mtext(title, at=5, line = 2, cex = 1)
    dev.off()
}
######### buffer
output_dir <- "results/06_manuscript/FigureS3/"
buffer_path <- "data/BE/Landscape/buffer/"
(buffer_pattern <- list.files(path = buffer_path, pattern = "_2008_LUagg2_1000\\.txt"))
list_buffer <- vector(mode = "list")
for (i in buffer_pattern){
    # print(i)
    # i <- "Sch_1960_LUagg2_2000.txt"
    buffer_file <- paste0(buffer_path, i)
    preffix <- unlist(strsplit(i, '\\.'))[1]
    region <- unlist(strsplit(preffix, '_'))[1]
    year <- unlist(strsplit(preffix, '_'))[2]
    resolution <- unlist(strsplit(preffix, '_'))[3]
    buffer <- unlist(strsplit(preffix, '_'))[4]
    title <- paste0(region, "_", year, "_", resolution, "_", buffer, "m")
    tmp_df <- read.delim(buffer_file, na.strings = c("NA", "", "Inf"),
                         sep = " ", stringsAsFactors = FALSE) %>% 
        dplyr::select(starts_with("L."), starts_with("Shannon."), ends_with("Gd"))
    
    # if(resolution=="LUagg1"){
    #   a <- tmp_df %>% dplyr::select(starts_with("L."), starts_with("Shannon."), ends_with("Gr"))
    # } else {
    #   a <- tmp_df %>% dplyr::select(starts_with("L."), starts_with("Shannon."), ends_with("Gd"))
    # }
    
    output_file <- paste0(output_dir, "buffer_", title, ".png")
    png(output_file, width = 12, height = 8, units = "cm", res = 600)
    
    # output_file <- paste0(output_dir, "buffer_", resolution, "/",title, ".tiff")
    # tiff(output_file, width = 1200, height = 962)
    M <- cor(tmp_df, use = "complete.obs")
    
    corrplot(M, method = 'number', type = "lower",
             tl.cex = 0.3, tl.col = "black",
             number.cex = 0.3, cl.cex = 0.3,
             mar = c(0, 0, 0, 0)
             # title = title,
             # tl.srt=0,
             # p.mat = M.sig$p, sig.level = .05
    )
    title <- paste0(region, "_", year)
    mtext(title, at=10, line = 1, cex = 1)
    dev.off()
}


########################### Table S1 ###########################
## this table is produced by 03_proxies_selected.R

########################### Table S2 ###########################
rm(list = ls())
resident_df <- read.delim("results/01_species_cover/resident_sp_mean_cover_each_plot_long.txt") %>% 
    mutate(region = case_when(str_detect(PlotID, "AEG") ~ "Alb_resident_sp",
                              str_detect(PlotID, "HEG") ~ "Hai_resident_sp",
                              str_detect(PlotID, "SEG") ~ "Sch_resident_sp")) %>% 
    dplyr::select(c(useful_name, region)) %>% 
    unique() %>% mutate(useful_name = str_trim(useful_name)) %>% 
    pivot_wider(names_from = region, values_from = region, values_fn = ~1, values_fill = 0)

BE_sp <- read.delim("results/01_species_cover/all_sp_yearly_cover_each_plot_long.txt") %>% 
    filter(Cover > 0) %>% 
    mutate(region = case_when(str_detect(PlotID, "AEG") ~ "Alb",
                              str_detect(PlotID, "HEG") ~ "Hai",
                              str_detect(PlotID, "SEG") ~ "Sch")) %>% 
    dplyr::select(c(Species, taxon_status, useful_name, region)) %>% 
    unique() %>% mutate(useful_name = str_trim(useful_name)) %>% 
    pivot_wider(names_from = region, values_from = region, values_fn = ~1, values_fill = 0) %>% 
    left_join(resident_df, by = "useful_name") %>% 
    mutate_all(~replace_na(., 0))

traits <- read.delim("results/02_traits/BE_species_traits_both.txt")
names(traits)[1] <- "useful_name"
all_info <- merge(BE_sp, traits, by = "useful_name", all.x = TRUE) %>% 
    dplyr::select(c(Species, useful_name,
                    usedPloidy, usedGS, GS_source, PL_source,
                    Alb, Hai, Sch, 
                    Alb_resident_sp, Hai_resident_sp, Sch_resident_sp)) %>% 
    rename(BE_record = Species, SpeciesName = useful_name, GenomeSize = usedGS, PloidyLevel=usedPloidy)

output_file <- "results/06_manuscript/TableS2.txt"
write.table(all_info, file = output_file, sep = "\t", 
            row.names = FALSE, col.names = TRUE)



########################### Figure S4 ###########################
# this plot is produced by script "04_stress_selected.R"

########################### Figure 1 ###########################
##### Figure 1a: A Venn diagram indicating the overlap of species found in different regions
rm(list = ls())
venn_sp_overlap <- function(sp_type){
    serach_pattern <- paste0("cover_mean_not_zero_", sp_type, "\\.txt")
    save_pattern <- "results/06_manuscript/resident_sp/Figure1a_resident_sp_region_venn.png"
    # save_pattern <- paste0("results/06_manuscript/", sp_type, "/Figure1a_region_venn.png")
    
    if(sp_type == "all_sp"){
        save_pattern <- "results/06_manuscript/all_sp/Figure1a_all_sp_region_venn.png"}
    
    input_dir <- "results/01_species_cover"
    
    (sp_occurred <- list.files(path = input_dir, 
                               pattern = serach_pattern))
    sp_occurrence <- vector(mode = "list")
    for (i in sp_occurred){
        region <- unlist(strsplit(i, "_"))[1]
        f <- paste0(input_dir, "/", i)
        a <- read.delim(f, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, check.names = FALSE)
        sp_occurrence[[region]] <- names(a)[-1]
    }
    venn_name <- paste(c("Alb:", "Hai:", "Sch:"),
                       c(length(sp_occurrence$AEG), 
                         length(sp_occurrence$HEG), 
                         length(sp_occurrence$SEG)))
    venn_sp <- ggVennDiagram(sp_occurrence, label_alpha = 0, 
                             set_size = 5, label_size = 5,
                             category.names = venn_name) + 
        scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
        scale_color_brewer(palette = 5) + 
        theme(legend.position = "none")
    
    venn_sp
    ggsave(
        save_pattern,
        venn_sp, width = 15, height = 9)                  
}
walk(c("all_sp", "resident_sp"), venn_sp_overlap)
##### Figure 1a: A Venn diagram indicating the overlap of species found in different regions


##### Figure 1b: histogram of species trait
# this plot is produced by script "02_KEW_traits.R"
##### Figure 1b: histogram of species trait

##### Figure 1c: violin plots for trait distribution (species level)
rm(list = ls())
trait_violin <- function(sp_type){
    # sp_type <- "resident_sp"
    regions <- c("AEG", "HEG", "SEG")
    regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
    region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)
    (input_file <- paste0("results/02_traits/regions_traits_", sp_type,".txt"))
    # "region"	"species"	"GS"	"GS_source"	"PL"	"PL_source"
    sp_type_trait <- read.delim(input_file, sep = "\t", 
                                stringsAsFactors = FALSE)
    
    plots <- function(sp_type_trait, trait){
        ylabel <- "Ploidy level"
        region_sp_num <- sp_type_trait %>% 
            filter(!is.na(PL)) %>% 
            count(region)
        trait_df <- sp_type_trait %>% 
            filter(!is.na(PL)) %>% 
            rename(sp_trait = "PL")
        text_pos <- max(trait_df$sp_trait) + 0.5
        
        if (trait == "GS"){
            ylabel <- "Genome size (2C/pg)"
            region_sp_num <- sp_type_trait %>% 
                filter(!is.na(GS)) %>% 
                count(region)
            trait_df <- sp_type_trait %>% 
                filter(!is.na(GS)) %>% 
                rename(sp_trait = "GS")
            text_pos <- max(trait_df$sp_trait) + 1
        }
        
        ggplot(trait_df, aes(region, sp_trait, color = region)) +
            geom_violin(size = 0.8) + geom_jitter(size = 1.4, alpha = 0.5) +
            geom_text(data = region_sp_num, aes(y = text_pos, label = n), size = 5, color = 'black') + 
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
    gs_plot <- plots(sp_type_trait, "GS")
    pl_plot <- plots(sp_type_trait, "PL")
    output <- paste0("results/06_manuscript/", sp_type, "/Figure1c_", sp_type, "_SP_level_trait_violin.png")
    png(output, width = 20, height = 12, units = "cm", res = 600)
    print(ggarrange(gs_plot, pl_plot))
    dev.off()
}
walk(c("all_sp", "resident_sp"), trait_violin)
##### Figure 1c: violin plots for trait distribution (species level)

##### Figure 1d: violin plots for trait distribution (community level) 
# this plot is produced by script "05_community_analysis_phylo.R"
##### Figure 1d: violin plots for trait distribution (community level)


########################### Figure 2 and 3 ###########################
# Forest plot for species level (Figure 2) and community level (Figure 3)
rm(list = ls())
forest_plot <- function(input_simple, input_quadratic, trait = "PL", level = "sp"){
    # trait <- "GS"; level = "sp"
    # input_simple <- "results/04_species_level/resident_sp/GS_phylolm_resident_sp_simple.txt"
    # input_quadratic <- "results/04_species_level/resident_sp/GS_phylolm_resident_sp_quadratic.txt"
    
    # trait <- "PL"; level = "sp"
    # input_simple <- "results/04_species_level/resident_sp/PL_phyloglm_resident_sp_simple.txt"
    # input_quadratic <- "results/04_species_level/resident_sp/PL_phyloglm_resident_sp_quadratic.txt"
    
    # trait <- "GS"; level = "community"
    # input_simple <- "results/05_community_level/all_sp/simple_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt"
    # input_quadratic <- "results/05_community_level/all_sp/complex_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt"
    
    # trait <- "PL"; level = "community"
    # input_simple <- "results/05_community_level/all_sp/simple_PL_lm_all_sp_yearly_mean.txt"
    # input_quadratic <- "results/05_community_level/all_sp/complex_PL_lm_all_sp_yearly_mean.txt"
    
    
    df <- read_delim(input_simple)
    df_quadratic <- read_delim(input_quadratic)
    df$Quadratic <- df_quadratic$complexed
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
    df2_1 <- insert_rows(df2, 3)
    df2_2 <- insert_rows(df2_1, 4)
    df2_3 <- insert_rows(df2_2, 5, directions = "up")
    df2_4 <- insert_rows(df2_3, 6, directions = "up")
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
               fill = rep(c("a", "b"), each = 7, length = nrow(.))) %>% 
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
    y_text <- str_wrap(df2$env_proxy[seq(2, nrow(df2), 3)], width = 1)
    
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
                           breaks=seq(4, nrow(df2_5), 7),
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
                           values = c("#f8766d", "#7cae00", "#00bfc4"),
                           labels = c("Alb", "Hai", "Sch"),
                           breaks = c("AEG", "HEG", "SEG")) +
        plot_theme +
        ggtitle(Title) 
    return(p)
}
##### sp
sp_inputfiles <- c("results/04_species_level/resident_sp/GS_phylolm_resident_sp_simple.txt",
                   "results/04_species_level/resident_sp/GS_phylolm_resident_sp_quadratic.txt",
                   "results/04_species_level/resident_sp/PL_glm_resident_sp_simple.txt",
                   "results/04_species_level/resident_sp/PL_glm_resident_sp_quadratic.txt")
inputfiles <- sp_inputfiles
level <- "sp"
sp_type <- "residenet_sp"
GS_plot <- forest_plot(inputfiles[1], inputfiles[2], trait = "GS", level)
output_dir <- "results/06_manuscript/"
output_gs <- paste0(output_dir, "/Figure2/Figure2_", level, "_", sp_type, "_GS_forest_plot.png")
ggsave(output_gs, width = 17, height = 20, units = "cm")
PL_plot <- forest_plot(inputfiles[3], inputfiles[4], trait = "PL", level)
output_pl <- paste0(output_dir, "/Figure2/Figure2_", level, "_", sp_type, "_PL_forest_plot.png")
ggsave(output_pl, width = 17, height = 20, units = "cm")

##### community 
community_inputfiles <- c("results/05_community_level/all_sp/simple_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt",
                          "results/05_community_level/all_sp/complex_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt",
                          "results/05_community_level/all_sp/simple_PL_lm_all_sp_yearly_mean.txt",
                          "results/05_community_level/all_sp/complex_PL_lm_all_sp_yearly_mean.txt")
inputfiles <- community_inputfiles
level <- "community"
sp_type <- "all_sp"
output_dir <- "results/06_manuscript/"
GS_plot <- forest_plot(inputfiles[1], inputfiles[2], trait = "GS", level)
output_gs <- paste0(output_dir, "/Figure3/Figure3_", level, "_", sp_type, "_GS_forest_plot.png")
ggsave(output_gs, width = 17, height = 20, units = "cm")
PL_plot <- forest_plot(inputfiles[3], inputfiles[4], trait = "PL", level)
output_pl <- paste0(output_dir, "/Figure3/Figure3_", level, "_", sp_type, "_PL_forest_plot.png")
ggsave(output_pl, width = 17, height = 20, units = "cm")


########################### Figure S5 (discussion overlap species) ###########################
##### overlap_sp
sp_inputfiles <- c("results/04_species_level/overlap_sp/GS_phylolm_resident_sp_overlap_sp_simple.txt",
                   "results/04_species_level/overlap_sp/GS_phylolm_resident_sp_overlap_sp_complex.txt",
                   "results/04_species_level/overlap_sp/PL_glm_resident_sp_overlap_sp_simple.txt",
                   "results/04_species_level/overlap_sp/PL_glm_resident_sp_overlap_sp_complex.txt")
inputfiles <- sp_inputfiles
level <- "sp"
sp_type <- "residenet_sp"
GS_plot <- forest_plot(inputfiles[1], inputfiles[2], trait = "GS", level)
output_dir <- "results/06_manuscript/"
output_gs <- paste0(output_dir, "/FigureS5/FigureS5_", level, "_", sp_type, "_GS_forest_plot.png")
ggsave(output_gs, width = 17, height = 20, units = "cm")
PL_plot <- forest_plot(inputfiles[3], inputfiles[4], trait = "PL", level)
output_pl <- paste0(output_dir, "/FigureS5/FigureS5_", level, "_", sp_type, "_PL_forest_plot.png")
ggsave(output_pl, width = 17, height = 20, units = "cm")


########################### Figure 4 ###########################
# this figure is produced by power analysis scripts (06_species_power_analysis.R 
# and 07_community_power_analysis.R)



