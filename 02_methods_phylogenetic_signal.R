#---------------------------------------------------------------------#
## Title ----
# Context-Dependent Relationships Between Genomic Traits and Plant Performance in Temperate Grasslands

## Purposes: access phylogenetic signal of genome size and ploidy level

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
# [1] geiger_2.0.11   caper_1.0.3     mvtnorm_1.2-4   MASS_7.3-60     phylolm_2.6.2  
# [6] phytools_2.0-3  maps_3.4.2      ape_5.7-1       lubridate_1.9.3 forcats_1.0.0  
# [11] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0    
# [16] tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] subplex_1.8             Rcpp_1.0.11             lattice_0.22-5         
# [4] listenv_0.9.0           digest_0.6.33           foreach_1.5.2          
# [7] utf8_1.2.4              parallelly_1.36.0       R6_2.5.1               
# [10] coda_0.19-4             pillar_1.9.0            rlang_1.1.2            
# [13] rstudioapi_0.15.0       phangorn_2.11.1         Matrix_1.6-4           
# [16] combinat_0.0-8          igraph_1.6.0            munsell_0.5.0          
# [19] compiler_4.2.0          numDeriv_2016.8-1.1     pkgconfig_2.0.3        
# [22] mnormt_2.1.1            optimParallel_1.0-2     globals_0.16.2         
# [25] tidyselect_1.2.0        expm_0.999-8            quadprog_1.5-8         
# [28] codetools_0.2-19        fansi_1.0.6             future_1.33.1          
# [31] tzdb_0.4.0              withr_2.5.2             grid_4.2.0             
# [34] nlme_3.1-164            gtable_0.3.4            lifecycle_1.0.4        
# [37] magrittr_2.0.3          scales_1.3.0            future.apply_1.11.1    
# [40] cli_3.6.2               stringi_1.8.3           scatterplot3d_0.3-44   
# [43] doParallel_1.0.17       generics_0.1.3          vctrs_0.6.5            
# [46] fastmatch_1.1-4         deSolve_1.40            iterators_1.0.14       
# [49] tools_4.2.0             glue_1.6.2              hms_1.1.3              
# [52] parallel_4.2.0          timechange_0.2.0        colorspace_2.1-0       
# [55] clusterGeneration_1.3.8


library(tidyverse)
library(phytools)
library(phylolm)
library(caper)
library(ape)
library(geiger)

rm(list = ls())
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
trait_df <- read.delim("results/02_traits/BE_species_traits_both.txt")
# View(trait_df)

gs_df <- trait_df %>% filter(!is.na(usedGS))
rownames(gs_df) <- gs_df$BEname
signal_tree <- phylo_trees(gs_df$BEname)
name_Ok <- name.check(signal_tree, gs_df)
if (name_Ok != 'OK'){stop("The species name is not matched")}
gs_df <- gs_df[match(signal_tree$tip.label, gs_df$BEname), ]
signal_trait <- setNames(gs_df$usedGS, gs_df$BEname)
# signal_trait <- setNames(log(gs_df$usedGS), gs_df$BEname)

signal_lambda <- phylosig(signal_tree, signal_trait, method = "lambda", test = TRUE)
# Phylogenetic signal lambda : 0.93461 
# logL(lambda) : -770.377 
# LR(lambda=0) : 179.933 
# P-value (based on LR test) : 5.01291e-41 
signal_K <- phylosig(signal_tree, signal_trait, test = TRUE, nsim = 10000)
# Phylogenetic signal K : 0.960292 
# P-value (based on 10000 randomizations) : 1e-04 

pl_df <- trait_df %>% filter(!is.na(usedPloidy)) %>% mutate(ploidy = if_else(usedPloidy > 2, 1, 0))
rownames(pl_df) <- pl_df$BEname
signal_tree <- phylo_trees(pl_df$BEname)
name_Ok <- name.check(signal_tree, pl_df)
if (name_Ok != 'OK'){stop("The species name is not matched")}
pl_df <- pl_df[match(signal_tree$tip.label, pl_df$BEname), ]
D <- phylo.d(pl_df, signal_tree, BEname, ploidy, permut = 10000, rnd.bias = NULL)
D
# Estimated D :  0.9066864
# Probability of E(D) resulting from no (random) phylogenetic structure :  0.1027
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0

