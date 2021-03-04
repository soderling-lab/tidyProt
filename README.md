# tidyProt

## Installation
To install the `tidyProt` package, type the following into an R console:
`devtools::install_github("soderling-lab/tidyProt")`.


## Key Dependencies
Insure you have installed the following R packages:

* lmerTest `install.packages("lmerTest")`
* dplyr `install.packages("dplyr")`
* data.table `install.packages("data.table")`
* reshape2 `install.packages("reshape2")`


## Usage

```R

# work with data from soderling-lab/SwipProteomics
# devtools::install_github("soderling-lab/SwipProteomics")

library(dplyr)
library(tidyProt)

library(SwipProteomics)  

data(washc) # wash complex proteins
data(wash_bioid) # the normalized proximity proteomics data

## linear model-based comparison
fx <- log2(Intensity) ~ 0 + Condition
fm <- lm(fx, data=wash_bioid %>% subset(Accession==washc[1]))

fm

# create a contrast
LT <- getContrast(fm, "WASH","Control")

# assess comparison
lmTestContrast(fm, LT)

## linear mixed-model based comparison
data(swip_tmt)

fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Protein)

fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Protein %in% washc))

LT <- getContrast(fm, "Mutant", "Control")

lmerTestContrast(fm, LT)

```
