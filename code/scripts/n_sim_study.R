# Script for generating simulation results investigating the effect of varying
# the correlation length scale parameter of the matern-3/2 covariance on the
# results.

# The sample_tstat.R script must have already been run to generate a sampling
# distribution of the t statistic under the null hypothesis.
# FIXTHIS set.seed()

source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")

set.seed(1615656706)

t_sample_n100 <- fread("data/t_sample_delta_1e-3_n100.csv")$t_stat
t_ecdf_n100 <- ecdf(t_sample_n100)

# create lat, lon grid on which to simulate data. Choose lat. and lon. ranges
# and choose units to create resolution of grid.
lat_lon_grid <- build_lat_lon_grid(lat_range = c(-60, 60),
                                   lon_range = c(0, 358),
                                   lat_unit = 2,
                                   lon_unit = 2)

# approximate distances are found by calculating an EW distance and a NS distance
# between each gridbox. See paper/function code for details. The same distance
# is calculated in anisotropic and isotropic cases for consistency.

szon_dist <- find_signed_zonal_dist(lat_lon_grid)
smer_dist <- find_signed_merid_dist(lat_lon_grid)

#isotropic distance measure
sim_dist <- sqrt(szon_dist^2 + smer_dist^2)

# matern-3/2 with variance paramter = 1 and correlation length parameters 1.54
# 
sim_cov <- matern_1.5_cov(sim_dist, phi = 1)
sqrt_cov <- sqrt_eigen_cov(sim_cov)

loc_list <- get_loc_list(lat_lon_grid)

# see documentation of run_sim() in functions/sim.R for more details.
# iters controls the number of independent spatial fields generated
# n controls number of independent samples per spatial field.
# delta gives the max value for rho * beta under the interval null
# beta is the true linear relationship between Y and the random location
# noise sd is the amount of random normal noise added to the random location


results_out_n100_sigma1 <- run_sim(iter = 10,
                            lat_lon_grid = lat_lon_grid,
                            t_ecdf = t_ecdf_n100,
                            szon_dist = szon_dist,
                            smer_dist = smer_dist,
                            loc_list = loc_list,
                            cov_mat = sim_cov,
                            sqrt_cov = sqrt_cov,
                            n = 100,
                            delta = 1e-3,
                            beta = 1,
                            noise_sd = 1,
                            verbose = TRUE)

fwrite(results_out_n100_sigma1, file = "data/iso_sim_study_n100_sigma1.csv")

results_out_n100_sigma2 <- run_sim(iter = 10,
                                  lat_lon_grid = lat_lon_grid,
                                  t_ecdf = t_ecdf_n100,
                                  szon_dist = szon_dist,
                                  smer_dist = smer_dist,
                                  loc_list = loc_list,
                                  cov_mat = sim_cov,
                                  sqrt_cov = sqrt_cov,
                                  n = 100,
                                  delta = 1e-3,
                                  beta = 1,
                                  noise_sd = 2,
                                  verbose = TRUE)

fwrite(results_out_n100_sigma2, file = "data/iso_sim_study_n100_sigma2.csv")

results_out_n100_sigma3 <- run_sim(iter = 10,
                                   lat_lon_grid = lat_lon_grid,
                                   t_ecdf = t_ecdf_n100,
                                   szon_dist = szon_dist,
                                   smer_dist = smer_dist,
                                   loc_list = loc_list,
                                   cov_mat = sim_cov,
                                   sqrt_cov = sqrt_cov,
                                   n = 100,
                                   delta = 1e-3,
                                   beta = 1,
                                   noise_sd = 3,
                                   verbose = TRUE)

fwrite(results_out_n100_sigma3, file = "data/iso_sim_study_n100_sigma3.csv")

results_out_n100_sigma1$n <- 100
results_out_n100_sigma1$sigma <- 1
results_out_n100_sigma1$phi <- 1
results_out_n100_sigma1$theta <- 0
results_out_n100_sigma1$psi <- 1

results_out_n100_sigma2$n <- 100
results_out_n100_sigma2$sigma <- 2
results_out_n100_sigma2$phi <- 1
results_out_n100_sigma2$theta <- 0
results_out_n100_sigma2$psi <- 1

results_out_n100_sigma3$n <- 100
results_out_n100_sigma3$sigma <- 3
results_out_n100_sigma3$phi <- 1
results_out_n100_sigma3$theta <- 0
results_out_n100_sigma3$psi <- 1

results_out_n100 <- rbind(results_out_n100_sigma1,
                          results_out_n100_sigma2,
                          results_out_n100_sigma3)

results_out_n100$base_FDX <- results_out_n100$base_FDP > 0.1
results_out_n100$smooth_FDX <- results_out_n100$smooth_FDP > 0.1

n100_mean_results <- results_out_n100[, lapply(.SD, mean), by = .(sigma,
                                                                  n,
                                                                  phi,
                                                                  theta,
                                                                  psi)]

n100_mean_results[, iteration := NULL]

n_sim_study_mean <- melt(n100_mean_results, 
                         id.vars = c("sigma", "n", "phi", "theta", "psi"),
                         measure.vars = measure(method,
                                                value.name,
                                                pattern = "^([[:alpha:]]*)_(.*)$"))

n_sim_study_mean <- n_sim_study_mean[order(sigma)]

fwrite(results_out_n100, file = "data/n_sim_study.csv")
fwrite(n_sim_study_mean, file = "data/n_sim_study_mean.csv")
