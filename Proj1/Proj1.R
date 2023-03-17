####Packages####
library(Rsamtools)
library(tidyverse)
library(Dict)

####Functions####
# get the complement of a given base
complement <- function(base){
  bases <- Dict$new(
    "A" = "T",
    "T" = "A",
    "C" = "G",
    "G" = "C"
  )
  return(bases[base])
}

####Inputs####
#generate posterior probabilities for possible genotypes, assuming biallelic model
#finish translating to R
#supply a dataframe with the reference and alternate alleles, minor allele freq, observed reads, and probability of error
posteriorProbs <- function(df){
  # calculate prior probabilities
  A = df$ref %>% unique
  B = df$alt %>% unique
  maf = df$maf %>% unique
  pAA = (1 - maf) * (1 - maf)
  pBB = maf * maf
  pAB = 1 - (pAA + pBB)
  # generate log likelihoods for each genotype
  llAA <- 0
  llBB <- 0
  llAB <- 0
  for (i in seq(nrow(df))){
    O = df$observations[i]
    E = df$probability_of_error[i]
    if (O == A){
      E = 1 - E
      }
    llAA = llAA + log(E)
    llBB = llBB + log(1 - E) # all error probabilities are just the complement of the AA genotypes
    llAB = llAB + log(0.5) # read's likelihood simplifies to 0.5 under assumption that P(S = A) = P(S = B)
    }
  # calculate denominator of Bayes' theorem equation
  d = exp(llAA + log(pAA)) + exp(llBB + log(pBB)) + exp(llAB + log(pAB))
  # calculate posterior probabilities for each genotype
  ppAA = exp(llAA + log(pAA) - log(d))
  ppBB = exp(llBB + log(pBB) - log(d))
  ppAB = exp(llAB + log(pAB) - log(d))
  ef <- tibble(
    "AA" = ppAA,
    "BB" = ppBB,
    "AB" = ppAB
  )
  return (ef)
}

#read in bam as a bam object
bamObj <- BamFile("Proj1/data/align_sort.bam")

#run pileup on bam object to get number of each read at the putatuve snp positions
pileupObj <- pileup(bamObj)
pileupDf <- as.data.frame(pileupObj)

#pull base qualities from bam
qualityDf <- data.frame(chr = scanBam(bamObj)[[1]]$rname, 
                        pos = scanBam(bamObj)[[1]]$pos,
                        qual = scanBam(bamObj)[[1]]$qual) %>% 
  mutate(qual = utf8ToInt(str_sub(qual, 1, 1)) - 33, chr = as.character(chr)) %>% 
  mutate(probability_of_error = 10**(- qual / 10)) %>% 
  as_tibble()

#create table giving read counts for each putative snp
snps <- read.table("Proj1/data/putatative_snps.tsv", sep = "\t", header = TRUE) %>%
  as_tibble()

snpsReads <- snps %>% 
  left_join(pileupDf, by = c("chr" = "seqnames", "pos"), multiple = "all") %>% 
  left_join(qualityDf, by = c("chr", "pos"), multiple = "all") %>% 
  mutate(strand = as.character(strand), nucleotide = as.character(nucleotide))

#modify snps so we find the complement of the bases in the reverse strand
snpsOriented <- snpsReads
snpsOriented$nucleotide <- apply(snpsReads, 1, function(r){
  if (r["strand"] == "-"){
    r["nucleotide"] = complement(r["nucleotide"])
  }
  return(r["nucleotide"])
})
snpsOriented <- snpsOriented %>% 
  mutate(strand = NULL)

#generate the complete table of reads
reads <- apply(snpsOriented, 1, function(r){
  tibble("pos" = r["pos"] %>% as.numeric(),
         "ref" = r["ref"],
         "alt" = r["alt"],
         "maf" = r["maf"] %>% as.numeric(),
         "observations" = rep(r["nucleotide"], r["count"]), 
         "probability_of_error" = qualityDf$probability_of_error[1])#r["probability_of_error"]) #FIX THIS
}) %>% 
  enframe(name = NULL) %>% 
  unnest(cols = value)

# output the posterior probabilities of each genotype for each position
result <- reads %>% 
  split(f = reads$pos) %>% 
  lapply(mutate, pos = NULL) %>% 
  lapply(posteriorProbs) %>% 
  enframe(name = "pos") %>% 
  mutate(pos = as.numeric(pos)) %>% 
  unnest(value) %>% 
  left_join(snps, by = "pos")

####Analysis####
# readGAlignments(get seq and qual)
# Granges(go to a particular pos)
