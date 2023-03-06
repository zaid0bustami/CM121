library(Rsamtools)
library(tidyverse)

bamObj <- BamFile("data/align_sort.bam")
pileupObj <- pileup(bamObj)
pileupDf <- as.data.frame(pileupObj)

snps <- read.table("data/putatative_snps.tsv", sep = "\t", header = TRUE) %>%
  as.data.frame() %>% 
  left_join(pileupDf, by = c("chr" = "seqnames", "pos"))
