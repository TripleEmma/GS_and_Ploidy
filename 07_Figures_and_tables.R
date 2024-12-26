## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: Produce figures and tables used in the manuscript

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

library(tidyverse)
library("scales")
library("ggVennDiagram")
library(ggrepel)
library(ggpubr)
library(data.table)
library(gtable)
library(grid)
library(gridExtra)
rm(list = ls())

regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)
regions_full <- setNames(c("Alb", "Hainich", "Schorfheide"), regions2)
########################### Figure S1 ###########################
# this figure is produced by script "02_traits_mean_median_diff.R"

########################### Figure S2 ###########################
# this figure is produced by script "02_multiple_trait_cover.R"

########################### Figure S3 ###########################
# FigureS3: correlation plot indicating landscape metrics are highly correlated with others
library(corrplot)
######### EP
output_dir <- "results/06_manuscript/FigureS3/"
EP_path <- "data/BE/Landscape/EP/"
(EP_pattern <- list.files(path = EP_path, pattern = "_2008_LUagg2_EP\\.txt"))
list_EP <- vector(mode = "list")
for (i in EP_pattern){
    # i <- "Sch_2008_LUagg2_EP.txt"
    print(i)
    EP_file <- paste0(EP_path, i)
    region <- unlist(strsplit(i, '_'))[1]
    year <- unlist(strsplit(i, '_'))[2]
    resolution <- unlist(strsplit(i, '_'))[3]
    title <- paste0(regions_full[region], "_", year)
    list_name <- paste0(region, "_", year, "_", resolution)
    tmp_df <- read.delim(EP_file, sep = " ", 
                         na.strings = c("NA", "", "Inf"),
                         stringsAsFactors = FALSE)
    print(unique(tmp_df$LUagg2))
    
    output_file <- paste0(output_dir, "EP_", list_name, ".png")
    png(output_file, width = 15, height = 8, units = "cm", res = 600)
    M <- cor(tmp_df[, 4:length(tmp_df)], use = "complete.obs")
    M.sig = cor.mtest(M)
    corrplot(M, method = 'number', type = "lower",
             tl.cex = 0.4, tl.col = "black",
             number.cex = 0.4, cl.cex = 0.4,
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
    png(output_file, width = 16, height = 12, units = "cm", res = 1000)
    
    # output_file <- paste0(output_dir, "buffer_", resolution, "/",title, ".tiff")
    # tiff(output_file, width = 1200, height = 962)
    M <- cor(tmp_df, use = "complete.obs")
    
    corrplot(M, method = 'number', type = "lower",
             tl.cex = 0.4, tl.col = "black",
             number.cex = 0.4, cl.cex = 0.4,
             mar = c(0, 0, 0, 0)
             # title = title,
             # tl.srt=0,
             # p.mat = M.sig$p, sig.level = .05
    )
    title <- paste0(regions_full[region], "_", year)
    mtext(title, at=10, line = 1, cex = 1)
    dev.off()
}

########################### Figure S4 ###########################
# this plot is produced by script "03_stress_selected.R"

########################### Figure S5 ###########################
# helper function
into_group <- function(z0, num_groups = 4){
    z <- z0 %>% filter(!is.na(cover))
    cover <- z$cover
    names(cover) <- z$PlotID # this is within speceis, so each plot would only have one cover value
    sorted_indices <- order(cover, decreasing = TRUE)
    sorted_vec <- na.omit(cover[sorted_indices])
    cumsum_vec <- cumsum(sorted_vec)
    indices_to_use <- sorted_indices[1:which(cumsum_vec >= 0.5)[1]]
    selected_plot_ids <- names(cover)[indices_to_use]
    
    plot_ids <- unique(z0$PlotID)
    grouping_factor <- cut(seq_along(plot_ids), breaks = num_groups, labels = FALSE)
    groups <- split(plot_ids, grouping_factor)
    group_match <- sapply(groups, function(group) {
        all(selected_plot_ids %in% group)
    })
    true_index <- which(group_match)
    if (length(true_index) > 0) {
        return(true_index[1])  # Ensure it returns only one index
    } else {
        return(5)
    }
}

relative_cover <- function(x){
    df <- read_delim(x)
    df <- df %>%
        mutate(across(-PlotID, ~ .x / sum(.x, na.rm = TRUE)))
    return(df)
}

species_cover_dis <- function(envName, env, relativeCover) {
    
    envNames <- c("Grazing", "Mowing", "Fertilization", "LUI",
                  "Patch.Area", "C.Hanski.Gd", "C.NN.dist.Gd",
                  "C.patch.area.w.Gd", "C.patch.no.Gd", "C.edge.dens.Gd")
    formalNames <- c("Grazing", "Mowing", "Fertilization", "LUI",
                     "Patch size", "Hanski connectivity", "NN distance",
                     "Area weighted patch size", "Patch number", "Edge density")
    names(formalNames) <- envNames
    
    regions <- c('AEG', 'HEG', 'SEG')
    regionFormalNames <- c('Alb', 'Hai', 'Sch')
    names(regionFormalNames) <- regions
    
    plotList <- vector(mode = 'list', length = 3)
    
    for (region in regions) {
        # region <- 'AEG'
        # envName <- 'Grazing'
        regionCover <- relativeCover[[region]]
        plotIDs <- regionCover$PlotID
        
        tmp <- env %>%
            dplyr::select(PlotID, all_of(envName)) %>%
            rename(envFactor = 2) %>%
            dplyr::filter(PlotID %in% plotIDs) %>%
            right_join(regionCover) %>%
            pivot_longer(-c(PlotID, envFactor),
                         names_to = 'species',
                         values_to = 'cover') %>%
            arrange(envFactor) %>%
            mutate(PlotID = factor(PlotID, levels = unique(PlotID)))
        
        result <- tmp %>%
            group_by(species) %>%
            summarize(group_num = into_group(pick(everything())), .groups = 'drop')
        # above is to give each species a tag (to order them on y-axis)
        
        tmp2 <- tmp %>% left_join(result, by = "species") %>%
            mutate(species = factor(species, levels = unique(species[order(group_num)])))
        
        plot_ids <- unique(tmp2$PlotID)
        num_groups <- 4
        grouping_factor <- cut(seq_along(plot_ids), breaks = num_groups, labels = FALSE)
        adjusted_cut_points <- seq_along(plot_ids)[tapply(seq_along(plot_ids), 
                                                          grouping_factor, 
                                                          function(x) tail(x, 1))] + 0.5
        
        plotTitle <- paste0(regionFormalNames[region], " ", formalNames[envName])
        
        # Create the plot
        pCover <- ggplot(data = tmp2, aes(x = PlotID, y = species, fill = cover)) +
            geom_tile(color = "lightblue") +
            scale_fill_gradientn(colors = c("yellow", "red"), na.value = "white") +
            geom_vline(xintercept = adjusted_cut_points[1:3]) +
            theme_minimal() +
            theme(axis.text.y = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(fill = "Relative\nCover") +
            xlab("Plot ID (Ordered by Environmental Value)") +
            ggtitle(plotTitle)
        
        # Add the plot to the list
        plotList[[which(regions == region)]] <- pCover
    }
    
    # Use plot_grid from cowplot to combine the plots
    combinedPlot <- cowplot::plot_grid(plotList[[1]], 
                                       plotList[[2]], 
                                       plotList[[3]], ncol = 1)
    
    figureName <- paste0("results/06_manuscript/FigureS5/", 
                         formalNames[envName], ".png")
    ggsave(figureName, combinedPlot, width = 8.40, height = 10)
}


coversPath <- 'results/01_species_cover'
coversType <- 'cover_mean_not_zero_resident_sp.txt'
covers <- list.files(path = coversPath, 
                     pattern = coversType, 
                     full.names = TRUE)

relativeCover <- purrr::map(covers, relative_cover)
regions <- c('AEG', 'HEG', 'SEG')
names(relativeCover) <- regions

env <- read_delim("results/03_proxies/all_proxies_sp.txt")
envNames <- c("Grazing", "Mowing", "Fertilization", "LUI",
              "Patch.Area", "C.Hanski.Gd", "C.NN.dist.Gd",
              "C.patch.area.w.Gd", "C.patch.no.Gd", "C.edge.dens.Gd")

walk(envNames, species_cover_dis, env, relativeCover)

########################### Figure S6 ###########################
### species level power analysis plot
library(cowplot)

rm(list = ls())
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), 
                         c("Alb", "Hai", "Sch"))

new_label <- c("Grazing", "Mowing", "Fertilization", "LUI", 
               "Patch_size", "Area_weighted_patch_size", 
               "Patch_number", "Edge_density", "NN_distance", 
               "Hanski_connectivity")

speciesPowerPlot <- function(trait){
    # trait <- 'GS'
    Level <- "(Species level)"
    inputDir <- "results/04_species_level/resident_sp/"
    inputFile <- "GS_power.txt"
    Title <- "Genome Size"
    
    if (trait != "GS"){
        inputFile <- "PL_power.txt"
        Title <- "Ploidy Level"
    }
    Title <- paste(Title, Level)
    input <- paste0(inputDir, inputFile)
    
    # input <- paste0("/Users/emma/Desktop/analysis_15/", inputDir, inputFile)
    
    df <- read.delim(input) %>% 
        mutate(region = case_when(region == "AEG" ~ "Alb",
                                  region == "HEG" ~ "Hai",
                                  region == "SEG" ~ "Sch"))
    
    a <- setNames(new_label, unique(unique(df$env_proxy)))
    df2 <- df %>% mutate(new_label = a[df$env_proxy])
    
    df2$new_label <- factor(df2$new_label, levels = new_label)
    
    plot_theme <- theme(# legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 16, face = "bold"),
        
        legend.key.size = unit(0.65, 'cm'),
        legend.spacing = unit(0.5, "cm")
    )
    
    
    p <- ggplot(df2, aes(r, power, color = region)) +
        geom_line(linewidth = 1) +
        ggtitle(Title) +
        xlab("Correlation") + ylab("Power") + 
        theme_bw() + 
        facet_wrap(~new_label) +
        scale_color_manual(name = "Region", values = region_color) +  # Corrected the variable name here
        theme(legend.position = "bottom")
    
    print(p) 
}
GS_sp <- speciesPowerPlot('GS')
ggsave("results/06_manuscript/FigureS6/species_gs_power.png", width = 12, height = 8)

PL_sp <- speciesPowerPlot('PL')
ggsave("results/06_manuscript/FigureS6/species_pl_power.png", width = 12, height = 8)

communityPowerPlot <- function(trait){
    # trait <- 'GS'
    Level <- "(Cover-weighted-Mean)"
    inputDir <- "results/05_community_level/all_sp/"
    inputFile <- "GS_power.txt"
    Title <- "Genome Size"
    
    if (trait != "GS"){
        inputFile <- "PL_power.txt"
        Title <- "Ploidy Level"
    }
    Title <- paste(Title, Level)
    input <- paste0(inputDir, inputFile)
    
    
    df <- read.delim(input) %>% 
        mutate(region = case_when(Region == "AEG" ~ "Alb",
                                  Region == "HEG" ~ "Hai",
                                  Region == "SEG" ~ "Sch"))
    
    plot_theme <- theme(# legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 16, face = "bold"),
        
        legend.key.size = unit(0.65, 'cm'),
        legend.spacing = unit(0.5, "cm")
    )
    
    p <- ggplot(df, aes(r, power, color = region)) +
        geom_line(linewidth = 1) +
        ggtitle(Title) +
        xlab("Correlation") + ylab("Power") + 
        theme_bw() + 
        # facet_wrap(~new_label) +
        scale_color_manual(name = "Region", values = region_color) +  # Corrected the variable name here
        theme(legend.position = "bottom")
    
    print(p) 
}
GS_community <- communityPowerPlot("GS")
PL_community <- communityPowerPlot("PL")
final_plot <- plot_grid(GS_community, PL_community, ncol = 2, rel_widths = c(1, 1))
ggsave("results/06_manuscript/FigureS6/community_power.png", width = 12, height = 8)

########################### Figure S7 ##########################
rm(list = ls())
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)

proxy_file <- "results/03_proxies/all_proxies_sp.txt"
sp_type <- "resident_sp"
cover_type <- paste0("cover_mean_not_zero_", sp_type, ".txt")
sp_cover_proxies <- vector(mode = "list")
proxy_df <- read.delim(proxy_file, stringsAsFactors = FALSE)
proxies <- names(proxy_df)[2:length(proxy_df)]
cover_dir <- "results/01_species_cover/"
(cover_file <- list.files(path = cover_dir, pattern = cover_type))

GS_ordered_sp <- read_delim("results/02_traits/regions_traits_resident_sp.txt") %>% 
    filter(!is.na(GS)) %>% 
    arrange(desc(GS)) %>% 
    dplyr::select(species, GS) %>% 
    distinct() %>% 
    pull(species)

for (i in cover_file){
    # i <- "SEG_cover_mean_not_zero_resident_sp.txt"
    region <- unlist(strsplit(i, split = "_"))[1]
    region_proxy <- proxy_df %>% filter(str_detect(PlotID, region)) %>% 
        dplyr::select(PlotID, Fertilization)
    (cover_region <- paste0(cover_dir, i))
    a <- read.delim(cover_region, stringsAsFactors = FALSE, check.names = FALSE)
    a <- a[match(region_proxy$PlotID, a$PlotID), ]
    b <- as.data.frame(apply(a[-1], 2, function(x){x/sum(x, na.rm = TRUE)})) # this gives the percentage of each species cover in each plot
    # a[-1] removed the first column
    rownames(b) <- a$PlotID
    # b <- b[, colSums(is.na(b)) < nrow(b)]
    b <- b[, colSums(is.na(b)) < 45]
    region_sp <- names(b)
    # region_target_sp <- intersect(target_sp, region_sp)[1:10]
    region_target_sp <- intersect(GS_ordered_sp, region_sp)[1:10]
    
    z <- b %>% 
        dplyr::select(all_of(region_target_sp)) %>% 
        rownames_to_column("PlotID") %>% 
        inner_join(region_proxy) %>% 
        pivot_longer(-c(PlotID, Fertilization), names_to = 'species', values_to = 'value') %>% 
        mutate(species = factor(species, levels = region_target_sp))
    
    p <- ggplot(z, aes(y=value, x=Fertilization)) + 
        geom_point(size = 1.6, color = region_color[region]) +
        facet_wrap(~species, scales = 'free_y', nrow = 1) +
        xlab('Standardized Fertilization') +
        ylab(paste0("Relative cover")) +
        ggtitle(regions2[region]) +
        theme_bw() + 
        theme(plot.title = element_text(size = 16, face = 'bold'),
              strip.text = element_text(face = 'italic', size = 10),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(linetype = "dashed"))
    
    output_dir <- "results/06_manuscript/FigureS7/"
    output <- paste0(output_dir, "FigureS7_", region, ".png")
    ggsave(output, width = 50, height = 8.5, units = "cm")
    # print(p)
}

########################## Figure S8 ###########################
library(stringr)
df <- read_csv("data/BE/LUI_not_standardized/LuiData/LUI_default components set__2024-12-19.txt") %>% 
    filter(!Year %in% c("1/1/2021 12:00:00 AM", "1/1/2022 12:00:00 AM")) %>% 
    mutate(year = as.character(sub(".*/(\\d{4}).*", "\\1", Year))) %>% 
    mutate(Exploratory=case_when(str_detect(EP_PlotID, "AEG") ~ "Alb", 
                                 str_detect(EP_PlotID, "HEG") ~ "Hai",
                                 str_detect(EP_PlotID, "SEG") ~ "Sch")) %>% 
    dplyr::select(year, Exploratory, EP_PlotID, TotalFertilization)

df2 <- df %>% group_by(EP_PlotID) %>% 
    summarise(yearlyMean = mean(TotalFertilization)) %>% 
    mutate(Exploratory=case_when(str_detect(EP_PlotID, "AEG") ~ "Alb", 
                                 str_detect(EP_PlotID, "HEG") ~ "Hai",
                                 str_detect(EP_PlotID, "SEG") ~ "Sch"))

df2$year <- "Mean"
df2$TotalFertilization <- df2$yearlyMean
df3 <- df2 %>% dplyr::select(year, Exploratory, EP_PlotID, TotalFertilization)
df4 <- rbind(df, df3)
df4$Exploratory <- factor(df4$Exploratory, levels = c("Alb", "Hai", "Sch"))

df4 %>% ggplot(aes(x = Exploratory, y = TotalFertilization, color = Exploratory)) +
    geom_boxplot(size = 0.8) + geom_jitter(size = 1.4, alpha = 0.5) + 
    scale_color_manual(values = c("#f8766d", "#7cae00", "#00bfc4")) +
    facet_wrap(~year, scales = 'free_y')  + 
    theme_bw() + xlab("") + ylab("Total Fertilization (Nitrogen kg/ha/year)") + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 18),
          strip.text = element_text(size = 16),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"))
ggsave("results/06_manuscript/FigureS8/yearlyNitrogen.png", width = 12, height = 14)








########################## Figure S9 ###########################
library(ggrepel)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)

forest_plot <- function(input_simple, input_quadratic = NA, 
                        global_simple, global_complex = NA,
                        trait = "PL", level = "sp"){
    
    if (!file.exists(input_simple)) stop(paste("File not found:", input_simple))
    if (!file.exists(global_simple)) stop(paste("File not found:", global_simple))
    if (!is.na(input_quadratic) && !file.exists(input_quadratic)) stop(paste("File not found:", input_quadratic))
    if (!is.na(global_complex) && !file.exists(global_complex)) stop(paste("File not found:", global_complex))
    
    
    # trait <- "GS"; level = "sp"
    # input_simple <- "results/04_species_level/resident_sp/GS_phylolm_resident_sp.txt"
    # global_simple <- "results/04_species_level/resident_sp/pooled_GS_phylolm.txt"
    
    # trait <- "PL"; level = "sp"
    # input_simple <- "results/04_species_level/resident_sp/PL_phylolm_resident_sp.txt"
    # global_simple <- "results/04_species_level/resident_sp/pooled_PL_phylolm.txt"
    
    # trait <- "GS"; level = "community"
    # input_simple <- "results/05_community_level/all_sp/simple_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt"
    # input_quadratic <- "results/05_community_level/all_sp/complex_phylogeny_controled_GS_lm_all_sp_yearly_mean.txt"
    
    # trait <- "PL"; level = "community"
    # input_simple <- "results/05_community_level/all_sp/simple_PL_lm_all_sp_yearly_mean.txt"
    # input_quadratic <- "results/05_community_level/all_sp/complex_PL_lm_all_sp_yearly_mean.txt"
    df_global <- read_delim(global_simple)
    if (level != "sp"){
        df_quadratic_global <- read_delim(global_complex)
        df_global$Quadratic <- df_quadratic_global$complexed
    }else{
        df_global$Quadratic <- 'N'
    }
    df_global$region <- 'All'
    df_global <- df_global %>% dplyr::select(c(region, env_proxy, estimate, pvalue, upper, lower, Quadratic))
    
    df <- read_delim(input_simple)
    if (level != "sp"){
        df_quadratic <- read_delim(input_quadratic)
        df$Quadratic <- df_quadratic$complexed
    }else{
        df$Quadratic <- 'N'
    }
    df <- df %>% dplyr::select(c(region, env_proxy, estimate, pvalue, upper, lower, Quadratic))
    
    df <- rbind(df, df_global)
    df$region <- factor(df$region, levels = c('AEG', 'HEG', 'SEG', 'All'))
    previous_names <- c("Grazing", "Mowing", "Fertilization", "LUI", "Patch.Area", 
                        "C.patch.area.w.Gd", "C.patch.no.Gd", "C.edge.dens.Gd", 
                        "C.NN.dist.Gd", "C.Hanski.Gd")
    new_names <- c("Grazing", "Mowing", "Fertilization", "LUI", "Patch_size", 
                   "  Area_weighted\npatch_size  ", "Patch_number", "Edge_density",
                   "NN_distance", "Hanski_index")
    
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
    
    add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.06
    if (level == 'sp' & trait == 'GS'){
        add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.02
    }
    if (level == 'sp' & trait == 'PL'){
        add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.045
    }
    if (level == 'community' & trait == 'PL'){
        add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.15
    }
    
    add_range <-  0.2
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
    
    if (trait == "GS" & level == "sp"){df3$sig_pos <- 0.06}
    if (trait == "GS" & level == "community"){df3$sig_pos <- 0.06}
    if (trait == "PL" & level == "sp"){df3$sig_pos <- 0.04}
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
sp_inputfiles <- c("results/04_species_level/resident_sp/GS_phylolm_resident_sp.txt",
                   "results/04_species_level/resident_sp/pooled_GS_phylolm.txt",
                   "results/04_species_level/resident_sp/PL_phylolm_resident_sp.txt",
                   "results/04_species_level/resident_sp/pooled_PL_phylolm.txt")

inputfiles <- sp_inputfiles
level <- "sp"
sp_type <- "residenet_sp"
GS_plot <- forest_plot(input_simple=sp_inputfiles[1], 
                       global_simple=sp_inputfiles[2],
                       trait = "GS", level = "sp")

output_dir <- "results/06_manuscript/FigureS9_pooled/"
output_gs <- paste0(output_dir, sp_type, "_pooled_GS_forest_plot.png")
ggsave(output_gs, width = 18, height = 22, units = "cm")

PL_plot <- forest_plot(input_simple=sp_inputfiles[3], 
                       global_simple=sp_inputfiles[4],
                       trait = "PL", level="sp")

output_pl <- paste0(output_dir, sp_type, "_pooled_PL_forest_plot.png")
ggsave(output_pl, width = 17, height = 22, units = "cm")


##### community_lmer 
community_inputfiles <- c("results/05_community_level/all_sp/simple_phylogeny_controled_GS_lm_yearly_mean.txt",
                          "results/05_community_level/all_sp/complex_phylogeny_controled_GS_lm_yearly_mean.txt",
                          "results/05_community_level/all_sp/pool_GS_lmer_simple.txt",
                          "results/05_community_level/all_sp/pool_GS_lmer_complex.txt",
                          "results/05_community_level/all_sp/simple_PL_lm_yearly_mean.txt",
                          "results/05_community_level/all_sp/complex_PL_lm_yearly_mean.txt",
                          "results/05_community_level/all_sp/pool_PL_lmer_simple.txt",
                          "results/05_community_level/all_sp/pool_PL_lmer_complex.txt")
inputfiles <- community_inputfiles
level <- "community"
sp_type <- "all_sp"
output_dir <- "results/06_manuscript/pooled/"
GS_plot <- forest_plot(inputfiles[1], inputfiles[2], 
                       inputfiles[3], inputfiles[4],
                       trait = "GS", level)
output_gs <- paste0(output_dir, sp_type, "_community_GS_forest_plot_lmer.png")
ggsave(output_gs, width = 18, height = 22, units = "cm")

PL_plot <- forest_plot(inputfiles[5], inputfiles[6], 
                       inputfiles[7], inputfiles[8],
                       trait = "PL", level)
output_pl <- paste0(output_dir, sp_type, "_community_PL_forest_plot_lmer.png")
ggsave(output_pl, width = 17, height = 22, units = "cm")



########################### Table S1 ###########################
## this table is produced by 03_repeatability.R

########################### Table S2 ###########################
## this table is produced by 03_proxies_selected.R

########################### Table S3 ###########################
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

output_file <- "results/06_manuscript/TableS3.txt"
write.table(all_info, file = output_file, sep = "\t", 
            row.names = FALSE, col.names = TRUE)


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
trait_boxplot <- function(sp_type){
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
            geom_boxplot(size = 0.8) + geom_jitter(size = 1.4, alpha = 0.5) +
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
    output <- paste0("results/06_manuscript/Figure1", "/Figure1c_", sp_type, "_SP_level_trait_boxplot.png")
    png(output, width = 20, height = 12, units = "cm", res = 600)
    print(ggarrange(gs_plot, pl_plot))
    dev.off()
}
walk(c("resident_sp"), trait_boxplot)
# walk(c("all_sp", "resident_sp"), trait_boxplot)

##### Figure 1c: violin plots for trait distribution (species level)

##### Figure 1d: box plots for trait distribution (community level) 
# this plot is produced by script "05_community_analysis.R"
##### Figure 1d: box plots for trait distribution (community level)



########################### Figure 2&3 ###########################
rm(list = ls())
regions <- c("AEG", "HEG", "SEG")
regions2 <- setNames(c("Alb", "Hai", "Sch"), regions)
region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), regions)


forest_plot <- function(input_simple, input_quadratic = NA, trait = "PL", level = "sp"){
    
    # input_simple <- "results/04_species_level/resident_sp/PL_phylolm_resident_sp.txt"; 
    # level <- "sp"; trait <- "PL"
    # input_simple <- "results/05_community_level/all_sp/simple_phylogeny_controled_GS_lm_yearly_mean.txt"
    # input_quadratic <- "results/05_community_level/all_sp/complex_phylogeny_controled_GS_lm_yearly_mean.txt"
    # level <- "community"; trait <- "GS"
    
    df <- read_delim(input_simple)
    # df_quadratic <- read_delim(input_quadratic)
    if (level != "sp"){
        complex <- read_delim(input_quadratic)
        df$Quadratic <- complex$complexed
    }else{
        df$Quadratic <- 'N'
    }
    
    previous_names <- c("Grazing", "Mowing", "Fertilization", "LUI", "Patch.Area", 
                        "C.patch.area.w.Gd", "C.patch.no.Gd", "C.edge.dens.Gd", 
                        "C.NN.dist.Gd", "C.Hanski.Gd")
    new_names <- c("Grazing", "Mowing", "Fertilization", "LUI", "Patch_size", 
                   "  Area_weighted\npatch_size  ", "Patch_number", "Edge_density",
                   "NN_distance", "Hanski_index")
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
    
    
    add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.05
    if (level == 'sp' & trait == 'GS'){
        add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.02
    }
    if (level == 'sp' & trait == 'PL'){
        add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.045
    }
    if (level == 'community' & trait == 'PL'){
        add_text <- max(abs(range(df2$lower, na.rm = TRUE))) + 0.15
    }
    
    
    add_range <-  0.2
    if (trait == "GS" &level != "sp"){
        add_range <- 0.3
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
    
    if (trait == "GS" & level == "sp"){df3$sig_pos <- 0.10}
    if (trait == "GS" & level == "community"){df3$sig_pos <- 0.11}
    if (trait == "PL" & level == "sp"){df3$sig_pos <- 0.06}
    if (trait == "PL" & level == "community"){df3$sig_pos <- 0.12}
    
    plot_theme <- theme(axis.text.y = element_text(size = 35, margin = margin(r = 8, unit = "pt")),
                        axis.text.x = element_text(size = 35),
                        axis.title = element_text(size = 39),
                        plot.title = element_text(size = 42, face = "bold"),
                        axis.ticks.y = element_line(linewidth = 2),
                        legend.text = element_text(size = 28), 
                        legend.title = element_text(size = 32), 
                        legend.position = "bottom")
    # legend.position = "none")
    
    if (trait == "PL") {
        plot_theme <- theme(axis.text.y = element_text(size = 35, margin = margin(r = 8, unit = "pt")),
                            axis.text.x = element_text(size = 35),
                            axis.title = element_text(size = 39),
                            plot.title = element_text(size = 42, face = "bold"),
                            axis.ticks.y = element_blank(),
                            legend.text = element_text(size = 28), 
                            legend.title = element_text(size = 32), 
                            legend.position = "bottom")
        # legend.position = "none")
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
                  size = 8, show.legend = FALSE) +
        # geom_text(aes(x = add_text + 0.01, 
        geom_text(aes(x = add_text + sig_pos, y = index - 0.5, label = sig),
                  hjust = 0, size = 10.5, show.legend = FALSE) +
        geom_line(data = df4, mapping = aes(linetype = Quadratic), linewidth = 3.5) +
        geom_point(size = 1, show.legend = FALSE) +
        geom_point(aes(x = estimate), size = 4.8) + 
        labs(x = expression(paste("Standardized coefficient (", beta, ")")),
             y = y_lab) + 
        guides(fill = "none") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_color_manual(name = "Region", 
                           values = c("#f8766d", "#7cae00", "#00bfc4"),
                           labels = c("Alb", "Hainich", "Schorfheide"),
                           breaks = c("AEG", "HEG", "SEG")) +
        plot_theme +
        ggtitle(Title) 
    # p
    # ggsave("test.png", width = 17, height = 19) 
    return(p)
}


output_dir <- "results/06_manuscript/"

##### sp
sp_type <- "residenet_sp"
level <- "sp"

sp_gs_inputfiles <- "results/04_species_level/resident_sp/GS_phylolm_resident_sp.txt"
sp_GS_plot <- forest_plot(sp_gs_inputfiles, trait = "GS", level = "sp")
output_gs <- paste0(output_dir, "/Figure2/Figure2_", level, "_", sp_type, "_GS_forest_plot.png")
ggsave(output_gs, width = 39, height = 49, units = "cm")

sp_pl_inputfiles <- "results/04_species_level/resident_sp/PL_phylolm_resident_sp.txt"
sp_PL_plot <- forest_plot(sp_pl_inputfiles, trait = "PL", level = "sp")
output_pl <- paste0(output_dir, "/Figure2/Figure2_", level, "_", sp_type, "_PL_forest_plot.png")
ggsave(output_pl, width = 39, height = 49, units = "cm")



##### community 
level <- "community"

community_gs_inputfiles <- c("results/05_community_level/all_sp/simple_phylogeny_controled_GS_lm_yearly_mean.txt",
                             "results/05_community_level/all_sp/complex_phylogeny_controled_GS_lm_yearly_mean.txt")
community_GS_plot <- forest_plot(community_gs_inputfiles[1], community_gs_inputfiles[2], 
                                 trait = "GS", level= "community")
output_gs_community <- paste0(output_dir, "/Figure3/Figure3_", level, "_", "phylo_control_GS_forest_plot.png")
ggsave(output_gs_community, width = 39, height = 49, units = "cm")



community_pl_inputfiles <- c("results/05_community_level/all_sp/simple_PL_lm_yearly_mean.txt",
                             "results/05_community_level/all_sp/complex_PL_lm_yearly_mean.txt")
community_PL_plot <- forest_plot(community_pl_inputfiles[1], community_pl_inputfiles[2], 
                                 trait = "PL", level= "community")
output_pl_community <- paste0(output_dir, "/Figure3/Figure3_", level, "_", "PL_forest_plot.png")
ggsave(output_pl_community, width = 39, height = 49, units = "cm")






