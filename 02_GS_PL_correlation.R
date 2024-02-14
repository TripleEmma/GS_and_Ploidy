library(tidyverse)

df <- read_delim("results/02_traits/BE_species_traits_both.txt") %>% 
    filter(!is.na(usedGS), !is.na(usedPloidy)) %>% 
    mutate(PL = case_when(usedPloidy<=2 ~ 2,
                          usedPloidy > 2 & usedPloidy<=4 ~ 4,
                          usedPloidy > 4 ~ 6)) %>% 
    select(BEname, PL, usedGS, usedPloidy) %>% 
    distinct()

nrow(df) # 229
cor.test(df$usedGS, df$usedPloidy)
# Pearson's product-moment correlation
# 
# data:  df$usedGS and df$usedPloidy
# t = 3.5877, df = 227, p-value = 0.0004085
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1051662 0.3507566
# sample estimates:
#       cor 
# 0.2316491 

summary(lm(df$usedGS~df$usedPloidy))
# lm(formula = df$usedGS ~ df$usedPloidy)
# 
# Residuals:
#     Min     1Q Median     3Q    Max 
# -9.082 -4.546 -2.755  1.764 41.680 
# 
# Coefficients:
#     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     3.3076     1.0349   3.196 0.001591 ** 
#     df$usedPloidy   0.9794     0.2730   3.588 0.000409 ***
#     ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 7.98 on 227 degrees of freedom
# Multiple R-squared:  0.05366,	Adjusted R-squared:  0.04949 
# F-statistic: 12.87 on 1 and 227 DF,  p-value: 0.0004085