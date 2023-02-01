library(tidyverse)

b <- exp(1)

df <- read.table("HW1/Data/reads.tsv", sep = "\t", header = TRUE) %>% 
  #calculate the likelihood for each read given a certain genotype
  mutate(l.AA = ifelse(observations == "A", 1 - probability_of_error, probability_of_error),
         l.TT = ifelse(observations == "T", 1 - probability_of_error, probability_of_error),
         l.AT = 0.5) %>% #likelihood simplifies to 0.5 under the assumption that P(S = A) = P(S = B)
  mutate(ll.AA = log(l.AA, base = b),
         ll.TT = log(l.TT, base = b),
         ll.AT = log(l.AT, base = b))

#calculate the likelihood of all reads given a certain genotype (adding the log liklihoods)
ll.AA <- sum(df$ll.AA)
ll.TT <- sum(df$ll.TT)
ll.AT <- sum(df$ll.AT)

#calculate the sum of these with the log of the genotype proportions

denom <- log(
  b^(ll.AA + log(0.95 ^ 2, base = b)) + 
  b^(ll.TT + log(0.05 ^ 2, base = b)) + 
  b^(ll.AA + log(1 - (0.95 ^ 2 + 0.05 ^ 2), base = b)),
  base = b
)

#abaondoned R becuase of denom underflow :/