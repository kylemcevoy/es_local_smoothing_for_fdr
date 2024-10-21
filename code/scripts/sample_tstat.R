# Script to generate null distribution of t-statistics.

source("code/functions/tstat_sampling_dist.R")

# for fwrite
library(data.table)

set.seed(96273)

# sample size of multivariate normal
n1 <- 50
n2 <- 75
n3 <- 100

# upper limit for interval null hypothesis
delta <- 1e-3

# true coefficient between Y and the selected gridbox
beta <- 1

#upper limit for rho
rho <- delta / beta

# t statistic generation
chol_var <- delta_chol(rho)
t_sample_list <- vector(mode = "list", length = 10L)

for (i in 1:10) {
  t_sample_list[[i]] <- replicate(1e6, draw_tstat(chol_S = chol_var,
                                                  n = n1,
                                                  beta = beta))

}
t_sample <- unlist(t_sample_list)

## write the csv for later use
fwrite(list(t_stat = t_sample), file = "data/t_sample_delta_1e-3_n50.csv")

t_sample_list2 <- vector(mode = "list", length = 10L)

for (i in 1:10) {
  t_sample_list2[[i]] <- replicate(1e6, draw_tstat(chol_S = chol_var,
                                                  n = n2,
                                                  beta = beta))
  
}

t_sample2 <- unlist(t_sample_list2)

fwrite(list(t_stat = t_sample2), file = "data/t_sample_delta_1e-3_n75.csv")

t_sample_list3 <- vector(mode = "list", length = 10L)

for (i in 1:10) {
  t_sample_list3[[i]] <- replicate(1e6, draw_tstat(chol_S = chol_var,
                                                  n = n3,
                                                  beta = beta))
  
}

t_sample3 <- unlist(t_sample_list3)

fwrite(list(t_stat = t_sample3), file = "data/t_sample_delta_1e-3_n100.csv")
