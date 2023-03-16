library(Rsamtools)
library(tidyverse)
library(Dict)

complement <- function(base){
  bases <- Dict$new(
    "A" = "T",
    "T" = "A",
    "C" = "G",
    "G" = "C"
  )
  return(bases[base])
}



#finish translating to R
posteriorProbs <- function(df, A, B, maf){
  pAA = (1 - maf) * (1 - maf)
  pBB = maf * maf
  pAB = 1 - (pAA + pBB)
  # generate log likelihoods for each genotype
  llAA <- 0
  llBB <- 0
  llAB <- 0
  for (i in 1:length(df)){
    O = df.at[i, 'observations']
    E = df.at[i, 'probability_of_error']
    if (O == A){
      E = 1 - E
      }
    llAA = llAA + math.log(E)
    llBB = llBB + math.log(1 - E) # all error probabilities are just the complement of the AA genotypes
    llAB = llAB + math.log(0.5) # read's likelihood simplifies to 0.5 under assumption that P(S = A) = P(S = B)
    }
  # calculate denominator of Bayes' theorem equation
  d = math.exp(llAA + math.log(pAA)) + math.exp(llBB + math.log(pBB)) + math.exp(llAB + math.log(pAB))
  # calculate posterior probabilities for each genotype
  ppAA = math.exp(llAA + math.log(pAA) - math.log(d))
  ppBB = math.exp(llBB + math.log(pBB) - math.log(d))
  ppAB = math.exp(llAB + math.log(pAB) - math.log(d))
  return (c(ppAA, ppBB, ppAB))
}

bamObj <- BamFile("Proj1/data/align_sort.bam")
pileupObj <- pileup(bamObj)
pileupDf <- as.data.frame(pileupObj)
qualityDf <- data.frame(scanBam(bamObj)[[1]]$rname, scanBam(bamObj)[[1]]$pos, scanBam(bamObj)[[1]]$qual)


snps <- read.table("Proj1/data/putatative_snps.tsv", sep = "\t", header = TRUE) %>%
  as_tibble() %>% 
  left_join(pileupDf, by = c("chr" = "seqnames", "pos"), multiple = "all") %>% 
  mutate(strand = as.character(strand), nucleotide = as.character(nucleotide))

snpsOriented <- snps
snpsOriented$nucleotide <- apply(snps, 1, function(r){
  if (r["strand"] == "-"){
    r["nucleotide"] = complement(r["nucleotide"])
  }
  return(r["nucleotide"])
})
snpsOriented <- snpsOriented %>% 
  mutate(strand = NULL)

readList <- apply(snpsOriented, 1, function(r){
  tibble("observations" = rep(r["nucleotide"], r["count"]), "probability_of_error" = 0.1)
})
