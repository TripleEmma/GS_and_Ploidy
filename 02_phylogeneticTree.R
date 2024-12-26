#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: check if all species with trait could be found in the tree

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
#     [1] Taxonstand_2.4  pbapply_1.7-2   ape_5.7-1       lubridate_1.9.3 forcats_1.0.0  
# [6] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0    
# [11] tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.11       pillar_1.9.0      compiler_4.2.0    tools_4.2.0      
# [5] digest_0.6.33     lifecycle_1.0.4   gtable_0.3.4      nlme_3.1-164     
# [9] timechange_0.2.0  lattice_0.22-5    pkgconfig_2.0.3   rlang_1.1.2      
# [13] cli_3.6.2         rstudioapi_0.15.0 parallel_4.2.0    withr_2.5.2      
# [17] generics_0.1.3    vctrs_0.6.5       hms_1.1.3         grid_4.2.0       
# [21] tidyselect_1.2.0  glue_1.6.2        R6_2.5.1          fansi_1.0.6      
# [25] tzdb_0.4.0        magrittr_2.0.3    scales_1.3.0      colorspace_2.1-0 
# [29] utf8_1.2.4        stringi_1.8.3     munsell_0.5.0 

library(tidyverse)
library(ape)
library(Taxonstand)
rm(list = ls())
tree_daphne <- read.tree(file = "data/phylogeny/DaPhnE_01.tre")
tree_sp <- str_replace(tree_daphne$tip.label, "_", " ")
# tpl_treeSP <- TPL(tree_sp)
# output <- paste0("data/phylogeny/tree_Taxonstand.txt")
# write.table(tpl_treeSP,
#             col.names = TRUE, row.names = FALSE, sep = "\t",
#             file = output)
tpl_treeSP <- read.delim("data/phylogeny/tree_Taxonstand.txt", 
                         stringsAsFactors = FALSE, check.names = FALSE)
tree_ori <- str_trim(paste(tpl_treeSP$Genus, tpl_treeSP$Species))
tree_tpl <- str_trim(paste(tpl_treeSP$New.Genus, tpl_treeSP$New.Species))

keepBE <- read.delim("results/01_species_cover/all_sp_mean_cover_each_plot_long.txt",
                     stringsAsFactors = FALSE) %>% 
    filter(cover_mean > 0, 
           taxon_status != "Uncertain", 
           taxon_status != "Accepted", 
           taxon_status != "keepTPL") %>% 
    distinct(useful_name) %>% 
    dplyr::pull(useful_name)

tree_corrected <- data.frame(original_name = tree_ori, 
                             tpl_name = tree_tpl, 
                             status = tpl_treeSP$Taxonomic.status) %>% 
    mutate(final = if_else(original_name %in% keepBE, original_name, tpl_name))
tree_corrected$original_name[tree_corrected$original_name == "Elytrigia repens"] <- "Elymus repens"

BE_sp <- read.delim("results/02_traits/regions_traits_all_sp.txt", 
                    stringsAsFactors = FALSE) %>% 
    filter(!is.na(PL) | !is.na(GS)) %>% pull(species)
if(length(setdiff(BE_sp, tree_corrected$final)) == 0){
    print("All species found in the tree!")
    output <- paste0("data/phylogeny/tree_Taxonstand_corrected.txt")
    write.table(tree_corrected,
                col.names = TRUE, row.names = FALSE, sep = "\t",
                file = output)
}  

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
