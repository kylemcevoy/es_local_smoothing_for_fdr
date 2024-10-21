# Script for generating simulation results for matern model misspecification.

# The sample_tstat.R script must have already been run to generate a sampling
# distribution of the t statistic under the null hypothesis.
set.seed(193259264)

source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")

t_sample <- fread("data/t_sample_delta_1e-3_n50.csv")$t_stat
t_ecdf <- ecdf(t_sample)

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

sim_dist <- sqrt(szon_dist^2 + smer_dist^2)

#matern-1/2 with variance paramter = 1 and correlation length scale = 1
sim_cov_0.5 <- matern_0.5_cov(sim_dist)
sqrt_cov_0.5 <- sqrt_eigen_cov(sim_cov_0.5)

#matern-5/2 with variance paramter = 1 and correlation length scale = 1
sim_cov_2.5 <- matern_2.5_cov(sim_dist)
sqrt_cov_2.5 <- sqrt_eigen_cov(sim_cov_2.5)

loc_list <- get_loc_list(lat_lon_grid)

# see documentation of run_sim for more details.
# iters controls the number of independent spatial fields generated
# n controls number of independent samples per spatial field.
# delta gives the max value for rho * beta under the interval null
# beta is the true linear relationship between Y and the random location
# noise sd is the amount of random normal noise added to the random location

results_out1 <- run_sim(iter = 10,
                        lat_lon_grid = lat_lon_grid,
                        t_ecdf = t_ecdf,
                        szon_dist = szon_dist,
                        smer_dist = smer_dist,
                        loc_list = loc_list,
                        cov_mat = sim_cov_0.5,
                        sqrt_cov = sqrt_cov_0.5,
                        n = 50,
                        delta = 1e-3,
                        beta = 1,
                        noise_sd = 1,
                        verbose = TRUE)

fwrite(results_out1, file = "data/iso_sim_study_n50_sigma1_matern0_5.csv")

results_out2 <- run_sim(iter = 10,
                        lat_lon_grid = lat_lon_grid,
                        t_ecdf = t_ecdf,
                        szon_dist = szon_dist,
                        smer_dist = smer_dist,
                        loc_list = loc_list,
                        cov_mat = sim_cov_2.5,
                        sqrt_cov = sqrt_cov_2.5,
                        n = 50,
                        delta = 1e-3,
                        beta = 1,
                        noise_sd = 1,
                        verbose = TRUE)

fwrite(results_out2, file = "data/iso_sim_study_n50_sigma1_matern2_5.csv")
