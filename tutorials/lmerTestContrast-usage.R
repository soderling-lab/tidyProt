#!/usr/bin/env Rscript

## fit WASH Complex

# load the WASHC UniProt identifiers
data(washc, package="SwipProteomics")

# module-level model includes ranef term Protein
fx <- log2(Intensity) ~ 0 + Condition + (1|Protein)

# fit the model
fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Protein %in% washc))

# assess overall 'Mutant-Control' comparison
res <- lmerTestContrast(fm, LT8) %>% 
	mutate(Contrast='Mutant-Control') %>% unique()
res


## Module-leve comparison for WASHC4 module

library(dplyr)
library(data.table)

library(tidyProt)

data(swip, package="SwipProteomics")
data(swip_tmt, package="SwipProteomics")
data(swip_partition, package="SwipProteomics")

fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Protein)

# Swip is in M38
partition[swip]

m38 <- names(which(partition==38))

fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Protein %in% m38))

LT = getContrast(fm, "Mutant", "Control")

# assess the contrast
lmerTestContrast(fm, LT) %>% 
	mutate(Contrast="Mutant-Control") %>% unique()
