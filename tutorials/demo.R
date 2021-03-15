#!/usr/bin/env Rscript

# imports

library(dplyr)
library(tidyProt) # soderling-lab/tidyProt


# load the iPSD BioID data

data(ipsd_bioid, package="Uezu2016")


# fit a simple linear model
fx <- log2(Intensity) ~ 0 + Condition
fm <- lm(fx, ipsd_bioid %>% subset(Symbol == "Gphn"))

# create a contrast for Gephyrin BioID versus control
LT <- getContrast(fm,"Gephyrin","ConditionControl")

# assess contrast
lmTestContrast(fm, LT) %>% knitr::kable()

# |Contrast                                 |   log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|
# |:----------------------------------------|--------:|--------------:|---------:|----------:|------:|--:|
# |ConditionGephyrin-BioID-ConditionControl | 8.235345|       301.3601| 0.0854091|   96.42231|      0| 12|
