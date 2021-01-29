<<<<<<< HEAD
## 12/21/2020

separate dependencies:
tidyProt: 
* lmTestContrast
* lmerTestContrast

getPPIs:
* data(musInteractome)

neten:
* neten

geneLists
* mapIDs==getIDs (from getPPIs)
* other mapping functions
* data('geneLists') # gene lists

=======
# tidyProt
_tidy statistical inference in protein mass spectrometry experiments_

## Installation
To install the `tidyProt` package, for example use `devtools::install_github("twesleyb/tidyProt")`.


## Key Dependencies
Insure you have installed the following R packages:

* dplyr `install.packages("dplyr")`
* impute `BiocManager::install("impute")`
* data.table `install.packages("data.table")`
* reshape2 `install.packages("reshape2")`
* lmerTest `install.packages("lmerTest")`


## Usage

```R

library(tidyProt)

# Linear Model Analysis
data(wash_bioid)

lmTestContrast()


# Mixed-Model Analysis
data(swip_tmt)
lmerTestConstrast()



```
>>>>>>> fa8284bc407769da82e663b3ba503fe495d237dc
