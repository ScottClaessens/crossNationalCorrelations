# Cross-national analyses require additional controls to account for the non-independence of nations

This repository contains the data and R code for the manuscript "Cross-national analyses require additional controls to account for the non-independence of nations".

## Getting Started

### Cloning this repository

To clone this repository using Git, you will need to first install [Git LFS](https://git-lfs.github.com/) to download large data files. Once Git LFS is installed, run:

```
git clone https://github.com/ScottClaessens/crossNationalCorrelations
```

### Installing required packages

To run this code, you will need to [install R](https://www.r-project.org/) and the following R packages:

```
install.packages(c("brms", "cowplot", "conleyreg", "countrycode", 
                   "dagitty", "geosphere", "ggcorrplot", "ggdag", "ggrepel", 
                   "ggtext", "haven", "lmtest", "papaja", "psych", "readxl", 
                   "rstan", "sf", "sjlabelled", "targets", "tarchetypes",
                   "tidybayes", "tidyverse"))
```

If you run into issues with the pipeline, try installing the specific R version or package versions outlined in `sessionInfo.txt`.

### Executing code

1. Set the working directory to this code repository
2. Load the `targets` package with `library(targets)`
3. To run all analyses, run `tar_make()`
4. To load individual targets into your environment, run `tar_load(targetName)`

## Help

Any issues, please email scott.claessens@gmail.com.

## Authors

Scott Claessens, scott.claessens@gmail.com
