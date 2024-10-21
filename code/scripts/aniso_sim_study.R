source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")

set.seed(567629070)

t_sample <- fread("data/t_sample_delta_1e-3_n50.csv")$t_stat
t_ecdf <- ecdf(t_sample)

# create lat, lon grid on which to simulate data. Choose lat. and lon. ranges
# and choose units to create resolution of grid.
lat_lon_grid <- build_lat_lon_grid(lat_range = c(-60, 60),
                                   lon_range = c(0, 358),
                                   lat_unit = 2,
                                   lon_unit = 2)

# number of gridboxes
p <- nrow(lat_lon_grid)

# approximate distances are found by calculating an EW distance and a NS distance
# between each gridbox. See paper/function code for details. The same distance
# is calculated in anisotropic and isotropic cases for consistency.

aniso_mat <- find_aniso_mat(aniso_angle = pi / 4, aniso_ratio = 2)
aniso_mat2 <- find_aniso_mat(aniso_angle = pi / 4, aniso_ratio = 1.5)

szon_dist <- find_signed_zonal_dist(lat_lon_grid)
smer_dist <- find_signed_merid_dist(lat_lon_grid)

sim_dist_pairs <- cbind(c(szon_dist), c(smer_dist))

aniso_dist <- sqrt(rowSums((sim_dist_pairs %*% aniso_mat) * sim_dist_pairs))
aniso_dist2 <- sqrt(rowSums((sim_dist_pairs %*% aniso_mat2) * sim_dist_pairs))

aniso_cov <- matrix(matern_1.5_cov(aniso_dist, phi = 1, sigma2 = 1),
                    nrow = p,
                    ncol = p)

aniso_cov2 <- matrix(matern_1.5_cov(aniso_dist2, phi = 1, sigma2 = 1),
                     nrow = p,
                     ncol = p)



sqrt_cov <- sqrt_eigen_cov(aniso_cov)
sqrt_cov2 <- sqrt_eigen_cov(aniso_cov2)

loc_list <- get_loc_list(lat_lon_grid)

# see documentation of run_sim for more details.
# iters controls the number of independent spatial fields generated
# n controls number of independent samples per spatial field.
# delta gives the max value for rho * beta under the interval null
# beta is the true linear relationship between Y and the random location
# noise sd is the amount of random normal noise added to the random location

aniso_results1 <- run_sim(iter = 10,
                          lat_lon_grid = lat_lon_grid,
                          t_ecdf = t_ecdf,
                          szon_dist = szon_dist,
                          smer_dist = smer_dist,
                          loc_list = loc_list,
                          cov_mat = aniso_cov,
                          sqrt_cov = sqrt_cov,
                          n = 50,
                          delta = 1e-3,
                          beta = 1,
                          noise_sd = 1,
                          verbose = TRUE)

fwrite(aniso_results1, file = "data/aniso_sim_study_0_25pi_psi2.csv")

aniso_results2 <- run_sim(iter = 10,
                          lat_lon_grid = lat_lon_grid,
                          t_ecdf = t_ecdf,
                          szon_dist = szon_dist,
                          smer_dist = smer_dist,
                          loc_list = loc_list,
                          cov_mat = aniso_cov2,
                          sqrt_cov = sqrt_cov2,
                          n = 50,
                          delta = 1e-3,
                          beta = 1,
                          noise_sd = 1,
                          verbose = TRUE)

fwrite(aniso_results2, file = "data/aniso_sim_study_0_25pi_psi1_5.csv")

aniso_results1$sigma <- 1
aniso_results1$n <- 50
aniso_results1$phi <- 1
aniso_results1$theta <- pi / 4
aniso_results1$psi <- 2

aniso_results2$sigma <- 1
aniso_results2$n <- 50
aniso_results2$phi <- 1
aniso_results2$theta <- pi / 4
aniso_results2$psi <- 1.5

aniso_results <- rbind(aniso_results2,
                       aniso_results1)

aniso_results$base_FDX <- aniso_results$base_FDP > 0.1
aniso_results$smooth_FDX <- aniso_results$smooth_FDP > 0.1

mean_aniso_results <- aniso_results[, lapply(.SD, mean), by = .(sigma,
                                                                n,
                                                                phi,
                                                                theta,
                                                                psi)]
aniso_results_mean <- melt(mean_aniso_results, 
                           id.vars = c("sigma", "n", "phi", "theta", "psi"),
                           measure.vars = measure(method,
                                                  value.name,
                                                  pattern = "^([[:alpha:]]*)_(.*)$"))

aniso_results_mean <- aniso_results_mean[order(psi)]

fwrite(aniso_results, file = "data/aniso_sim_study.csv")
fwrite(aniso_results_mean, file = "data/aniso_sim_study_mean.csv")
