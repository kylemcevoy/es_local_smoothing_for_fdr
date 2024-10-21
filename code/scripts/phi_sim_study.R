# Script for generating simulation results investigating the effect of varying
# the correlation length scale parameter of the matern-3/2 covariance on the
# results.

# The sample_tstat.R script must have already been run to generate a sampling
# distribution of the t statistic under the null hypothesis.

source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")

set.seed(0816153725)

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

#isotropic distance measure
sim_dist <- sqrt(szon_dist^2 + smer_dist^2)

# matern-3/2 with variance paramter = 1 and correlation length parameters 1.54
# 
sim_cov_phi_1.54 <- matern_1.5_cov(sim_dist, phi = 1.54)
sqrt_cov_phi_1.54 <- sqrt_eigen_cov(sim_cov_phi_1.54)

loc_list <- get_loc_list(lat_lon_grid)

# see documentation of run_sim() in functions/sim.R for more details.
# iters controls the number of independent spatial fields generated
# n controls number of independent samples per spatial field.
# delta gives the max value for rho * beta under the interval null
# beta is the true linear relationship between Y and the random location
# noise sd is the amount of random normal noise added to the random location

results_out_phi_1.54 <- run_sim(iter = 10,
                               lat_lon_grid = lat_lon_grid,
                               t_ecdf = t_ecdf,
                               szon_dist = szon_dist,
                               smer_dist = smer_dist,
                               loc_list = loc_list,
                               cov_mat = sim_cov_phi_1.54,
                               sqrt_cov = sqrt_cov_phi_1.54,
                               n = 50,
                               delta = 1e-3,
                               beta = 1,
                               noise_sd = 1,
                               verbose = TRUE)

fwrite(results_out_phi_1.54, file = "data/iso_sim_study_n50_sigma1_phi1_54.csv")

results_out_phi_1.54$sigma <- 1
results_out_phi_1.54$n <- 50
results_out_phi_1.54$phi <- 1.54
results_out_phi_1.54$theta <- 0
results_out_phi_1.54$psi <- 1

results_out_phi_1.54$base_FDX <- results_out_phi_1.54$base_FDP > 0.1
results_out_phi_1.54$smooth_FDX <- results_out_phi_1.54$smooth_FDP > 0.1

phi_results_mean <- results_out_phi_1.54[, lapply(.SD, mean), by = .(sigma,
                                                                     n,
                                                                     phi,
                                                                     theta,
                                                                     psi)]
phi_results_mean[, iteration := NULL]

phi_results_mean <- melt(phi_results_mean, 
                         id.vars = c("sigma", "n", "phi", "theta", "psi"),
                         measure.vars = measure(method,
                                                value.name,
                                                pattern = "^([[:alpha:]]*)_(.*)$"))

fwrite(results_out_phi_1.54, file = "data/phi_sim_study.csv")

fwrite(phi_results_mean, file = "data/phi_sim_study_mean.csv")



