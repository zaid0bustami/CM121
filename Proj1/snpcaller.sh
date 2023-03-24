#!/bin/bash

#install dependencies
sudo apt install r-cran-dplyr
sudo apt install r-cran-stringr
sudo apt install r-cran-tidyr
sudo apt install r-cran-tibble
sudo apt-get install liblzma-dev
sudo apt-get install libbz2-dev

# assign command-line arguments to variables
bamName=$1
metaName=$2

Rscript Proj1.R $bamName $metaName