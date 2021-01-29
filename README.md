# tidyProt

## Installation
To install the `tidyProt` package, for example use 
`devtools::install_github("soderling-lab/tidyProt")`.


## Key Dependencies
Insure you have installed the following R packages:

* dplyr `install.packages("dplyr")`
* data.table `install.packages("data.table")`
* reshape2 `install.packages("reshape2")`
* lmerTest `install.packages("lmerTest")`


## Usage

```R

# work with data from soderling-lab/SwipProteomics
# devtools::install_github("soderling-lab/SwipProteomics")

library(tidyProt)
library(SwipProteomics)  

# Linear Model Analysis
data(wash_bioid)

lmTestContrast()


# Mixed-Model Analysis
data(swip_tmt)

lmerTestConstrast()
```
