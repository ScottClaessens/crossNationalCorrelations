# The non-independence of nations and why it matters

This repository contains the R code for the manuscript "The non-independence of nations and why it matters".

## Getting Started

### Installing

To run this code, you will need to [install R](https://www.r-project.org/) and the following R packages:

```
install.packages(c("brms", "cowplot", "conleyreg", "countrycode", 
                   "dagitty", "geosphere", "ggdag", "ggrepel", "ggtext", 
                   "haven", "lmtest", "papaja", "psych", "readxl", 
                   "rstan", "sjlabelled", "targets", "tarchetypes",
				   "tidybayes", "tidyverse"))
```

### Executing code

1. Set the working directory to this code repository
2. Load the `targets` package with `library(targets)`
3. To run all analyses, run `tar_make()`
4. To load individual targets into your environment, run `tar_load(targetName)`

## Help

Any issues, please email scott.claessens@gmail.com.

## Authors

Scott Claessens, scott.claessens@gmail.com
