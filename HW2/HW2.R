#### libraries ####
library(tidyverse)

#
#### 1(a) ####
pj <- read.table("HW2/p.tsv")$V1
Nj <- c(15000, 30000, 150000)
j <- c(1, 2, 3)
# matrix containing the simulations
set.seed(123)
raw.counts <- sapply(j, function(sample){
  Xj = rmultinom(1, Nj[sample], pj)
  return(Xj)
})
colnames(raw.counts) <- paste0("X", j)
# scatterplot
raw.counts %>% 
  as_tibble() %>% 
  rownames_to_column() %>% 
  rename(gene = rowname) %>% 
  ggplot(aes(x = X1, y = X3)) +
  geom_point() +
  labs(title = "Raw Counts for X3 vs X1") +
  geom_abline(aes(intercept = 0, slope = 1), color = "red") +
  geom_abline(aes(intercept = 5, slope = Nj[3]/Nj[1]), color = "blue") +
  geom_text(x=200, y=250, label="y = x", angle=8, color = "red") +
  geom_text(x=200, y=2100, label="y = 10x", angle=30, color = "blue")


#
#### 1(b) ####
normalizers <- function(mat){
  sjg <- apply((mat + 1), 1, function(row){ #mat + 1 adds a psuedocount of 1
    N <- length(row)
    denom <- prod(row) ^ (1 / N) #calculate the denominator
    x <- sapply(row, function(Xjg){ #divide each jth entry in the row by denom
      return(Xjg / denom)
    })
    return(x)
  }) %>% t
  colnames(sjg) <- paste0("s", c(1:ncol(sjg)), "g") #rename columns
  sj <- apply(sjg, 2, median) #calculate the median s value for each j
  return(list(sjg, sj)) #returns a list with the sj values for each gene and sample specific normalization factors
  }
z <- normalizers(raw.counts)

#
#### 1(c) ####
z[[1]] %>% 
  as_tibble() %>% 
  mutate(s2g = NULL) %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "sjg_value") %>% 
  ggplot(aes(x = sjg_value, fill = sample))+
  geom_histogram(bins = 75) +
  # facet_wrap("sample") + # if you want two separate histograms
  scale_fill_viridis_d(end = 0.5) +
  labs(title = "Histograms of s3g vs s1g Values")
z[[1]] %>% 
  as_tibble() %>% 
  ggplot(aes(x = s1g, y = s3g))+
  geom_point() +
  labs(title = "Scatterplot of s3g vs s1g Values") +
  scale_x_continuous(limits = c(0, 6))


#
#### 1(d) ####
normalized.counts <- apply(raw.counts, 1, function(row){
  x = row / z[[2]]
  return(x)
}) %>% t
colnames(normalized.counts) <- paste0("Y", j)
normalized.counts %>% 
  as_tibble() %>% 
  rownames_to_column() %>% 
  rename(gene = rowname) %>% 
  ggplot(aes(x = Y1, y = Y3)) +
  geom_point() +
  labs(title = "Normalized Counts for Y3 vs Y1") +
  geom_abline(aes(intercept = 0, slope = 1), color = "red") +
  geom_text(x=550, y=580, label="y = x", angle=30, color = "red")

#
#### 1(e) ####
Nj.2 <- c(1e6, 1e6, 1e6)
# matrix containing the simulations
set.seed(123)
raw.counts.2 <- sapply(j, function(sample){
  Xj = rmultinom(1, Nj.2[sample], pj)
  return(Xj)
})
colnames(raw.counts.2) <- paste0("X", j)
# estimating sjg and sj
z.2 <- normalizers(raw.counts.2)

#
#### 1(f) ####
qj <- read.table("HW2/q.tsv")$V1
Nj.3 <- c(1e6, 1e6, 1e6)
# matrix containing the simulations
set.seed(123)
raw.counts.3 <- sapply(j, function(sample){
  Xj = rmultinom(1, Nj.3[sample], qj)
  return(Xj)
})
colnames(raw.counts.3) <- paste0("X", j)
# estimating sjg and sj
z.3 <- normalizers(raw.counts.3)

#
#### 1(g) ####
raw.counts.4 <- cbind(raw.counts.2, raw.counts.3)
j.4 <- 1:6
colnames(raw.counts.4) <- paste0("X", j.4)
z.4 <- normalizers(raw.counts.4)
normalized.counts.4 <- apply(raw.counts.4, 1, function(row){
  x = row / z.4[[2]]
  return(x)
}) %>% t
colnames(normalized.counts.4) <- paste0("Y", j.4)
normalized.counts.4 %>% 
  as_tibble() %>% 
  rownames_to_column() %>% 
  rename(gene = rowname) %>% 
  ggplot(aes(x = Y1, y = Y6)) +
  geom_point() +
  labs(title = "Normalized Counts for Y6 vs Y1") +
  geom_abline(aes(intercept = 0, slope = 1), color = "red") +
  geom_text(x=10500, y=12500, label="y = x", angle=10, color = "red")
cbind(pj, qj) %>% 
  as_tibble() %>% 
  ggplot(aes(x = pj, y = qj)) +
  geom_point() +
  labs(title = "Scatterplot of qj vs pj")
cbind(pj, qj) %>% 
  as_tibble() %>% 
  pivot_longer(cols = c("pj", "qj"), names_to = "distribution", values_to = "probability") %>% 
  ggplot(aes(x = probability, fill = distribution)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis_d() +
  labs(title = "Distributions of pj and qj") +
  scale_x_sqrt()
