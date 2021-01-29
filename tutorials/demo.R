#!/usr/bin/env Rscript

# imports

library(dplyr)
library(tidyProt) # twesleyb/tidyProt


# load the iPSD BioID data

data(ipsd_bioid)


# fit a simple linear model
fx <- log2(Intensity) ~ 0 + Condition
fm <- lm(fx, ipsd_bioid %>% subset(Symbol == "Gphn"))

# create a contrast for Gephyrin BioID versus control
LT <- getContrast(fm,"Gephyrin","ConditionControl")

# assess contrast
lmTestContrast(fm, LT) %>% knitr::kable()

# Call:
# lm(formula = fx, data = ipsd_bioid %>% subset(Symbol == "Gphn"))
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.14781 -0.08432  0.01901  0.06472  0.19906 
# 
# Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)    
# ConditionCollybistin-BioID 35.96282    0.06039   595.5   <2e-16 ***
# ConditionControl           29.89757    0.06039   495.0   <2e-16 ***
# ConditionGephyrin-BioID    38.13291    0.06039   631.4   <2e-16 ***
# ConditionInSyn1-BioID      35.53459    0.06039   588.4   <2e-16 ***
# ConditionQC                36.61036    0.04678   782.6   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1046 on 12 degrees of freedom
# Multiple R-squared:      1,	Adjusted R-squared:      1 
# F-statistic: 3.914e+05 on 5 and 12 DF,  p-value: < 2.2e-16


# |Contrast                                 |   log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|
# |:----------------------------------------|--------:|--------------:|---------:|----------:|------:|--:|
# |ConditionGephyrin-BioID-ConditionControl | 8.235345|       301.3601| 0.0854091|   96.42231|      0| 12|
