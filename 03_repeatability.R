## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: Test the repeatability of the land use practices across years 2006-2020

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
library(lme4)
rm(list = ls())
df <- read.csv("data/BE/LUI/separate_2006-2020/LUI_default components set_regional_separately_2022-12-02.txt") %>% 
    mutate(Region = case_when(str_detect(PLOTID, "AEG") ~ "Alb", 
                              str_detect(PLOTID, "HEG") ~ "Hai",
                              str_detect(PLOTID, "SEG") ~ "Sch")) %>% 
    mutate(years = as.character(gsub("\\D", "", YEAR))) %>% 
    rename(Fertilization = F_STD, Grazing = G_STD, Mowing = M_STD) %>% 
    pivot_longer(-c(PLOTID,YEAR,EXPLO, Region, years), values_to = 'value', names_to = 'proxy') 

df$EXPLO <- NULL
df$YEAR <- NULL
df$proxy <- factor(df$proxy, levels = c("Grazing", "Mowing", "Fertilization", "LUI"))


# repeatability
results <- list()
data <- df
# Loop through each region (EXPLO) and proxy
for (region in unique(data$Region)) {
    for (management in unique(data$proxy)) {
        # Subset data for this region and proxy
        subset_data <- data %>%
            filter(Region == region, proxy == management)
        
        # Fit a mixed-effects model
        model <- lmer(value ~ (1 | PLOTID), data = subset_data)
        
        # Extract variance components
        var_components <- as.data.frame(VarCorr(model))
        sigma_among <- var_components$vcov[1] # Among-plot variance
        sigma_within <- attr(VarCorr(model), "sc")^2 # Within-plot variance
        
        # Calculate repeatability
        repeatability <- sigma_among / (sigma_among + sigma_within)
        
        # Save results
        results[[paste(region, management, sep = "_")]] <- data.frame(
            Region = region,
            Proxy = management,
            Repeatability = repeatability
        )
    }
}

# Combine results into a single data frame
results_df <- bind_rows(results)

# View results
print(results_df)


region_color <- setNames(c("#f8766d", "#7cae00", "#00bfc4"), 
                         c("Alb", "Hai", "Sch"))

ggplot(df, aes(x = years, y = value, group =PLOTID, color = Region)) + 
    geom_line() +
    scale_color_manual(values = region_color) + 
    facet_grid(proxy ~ Region, scales = "free") +
    theme_bw() + xlab("Years") + ylab("Intensity") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 12),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"))
ggsave("results/06_manuscript/FigureS_timeSeries.png", width = 15, height = 9)    





