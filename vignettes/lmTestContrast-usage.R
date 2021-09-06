#!/usr/bin/env Rscript

## using lmTestContrast to assess differential abundance

library(dplyr)
library(data.table)

library(tidyProt)

data(gphn, package="Uezu2016")
data(ipsd_bioid, package="Uezu2016")

# drop QC data before fitting models
tidy_prot <- ipsd_bioid %>% filter(Condition != "QC")

fx <- log2(Intensity) ~ 0 + Condition

# NOTE: by setting the intercept to 0,
# we explicitly estimate all levels of Condition

# fit the model to protein data for gephyrin
fm <- lm(fx, data = tidy_prot %>% subset(Protein == gphn))

# create a contrast
LT <- tidyProt::getContrast(fm, "Gephyrin","Control")
LT

# assess contrast for gphn 
lmTestContrast(fm, LT) 
