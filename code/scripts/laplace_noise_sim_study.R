# Script for generating simulation results for laplace heavier tailed noise
# rather than normal noise.

# The sample_tstat.R script must have already been run to generate a sampling
# distribution of the t statistic under the null hypothesis.
set.seed(580103789)


library(extraDistr)

source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")

run_sim_laplace <- function(iter,
                    lat_lon_grid,
                    t_ecdf,
                    szon_dist,
                    smer_dist,
                    loc_list,
                    cov_mat,
                    beta,
                    delta,
                    n,
                    sqrt_cov,
                    noise_scale_param,
                    verbose = FALSE) {
  
  p <- nrow(lat_lon_grid)
  
  results_mat <- matrix(NA, nrow = iter * 1000, ncol = 13)
  
  init_cond_mat <- matrix(c(1, 1, 1,
                            0.2 * pi, 0.5 * pi, 0.8 * pi,
                            1.5, 1.5, 1.5),
                          nrow = 3)
  
  MLE_df <- data.frame(phi = rep(NA, p),
                       ani_angle = rep(NA, p),
                       ani_ratio = rep(NA, p),
                       convergence = rep(NA, p))
  
  rnd_indx <- c(replicate(iter, sample(p, size = 1000, replace = FALSE)))
  
  
  for (i in 1:iter) {
    if (verbose) print(i)
    
    sim_data <- gen_data_eigen(n = n, sqrt_eigen_cov = sqrt_cov)
    sim_data_std <- scale(sim_data)
    
    for (j in 1:p) {
      MLE_df[j, ] <- find_local_MLE(loc = j,
                                    coord_dt = lat_lon_grid,
                                    zonal_dist_mat = szon_dist,
                                    merid_dist_mat = smer_dist,
                                    data_mat_std = sim_data_std,
                                    lik_func = m1.5_lik,
                                    return_best = TRUE,
                                    init_cond = init_cond_mat,
                                    return_likelihood = FALSE)
      
      if (verbose && j %% 1000 == 0) print(j)
    }
    MLE_mat <- data.matrix(MLE_df[, 1:3])
    
    cov_weights <- find_smooth_weights(loc_list = loc_list,
                                       zonal_sdist = szon_dist,
                                       merid_sdist = smer_dist,
                                       cov_func = matern_1.5_cov,
                                       mle_params = MLE_mat,
                                       progress = FALSE)
    
    smooth_data <- sim_data_std %*% cov_weights
    
    for (k in 1:1000) {
      if (verbose && k %% 100 == 0) print(k)
      
      true_cov <- cov_mat[rnd_indx[(i - 1) * 1000 + k], ]
      true_null <- true_cov * beta < delta
      
      Y <- sim_data[, rnd_indx[(i - 1) * 1000 + k]] + 
        rlaplace(n, sigma = noise_scale_param)
      
      base_lm_fits <- fit_lm_tstat(y = Y, data_mat = sim_data_std)
      smooth_lm_fits <- fit_lm_tstat(y = Y, data_mat = smooth_data)
      
      base_pvals <- 1 - t_ecdf(base_lm_fits)
      smooth_pvals <- 1 - t_ecdf(smooth_lm_fits)
      
      base_reject <- find_BH_reject(base_pvals, alpha = 0.1)
      smooth_reject <- find_BH_reject(smooth_pvals, alpha = 0.1)
      
      base_metrics <- estim_metrics(p_dt = base_reject,
                                    true_null = true_null)
      smooth_metrics <- estim_metrics(p_dt = smooth_reject,
                                      true_null = true_null)
      
      results_mat[(i - 1) * 1000 + k, ] <- c(i,
                                             base_metrics,
                                             smooth_metrics)
    }
  }
  
  results_frame <- data.table(results_mat)
  
  metric_names <- c("reject",
                    "true_reject",
                    "FDP",
                    "FNP",
                    "sens",
                    "spec")
  
  names(results_frame) <- c("iteration",
                            paste0("base_", metric_names),
                            paste0("smooth_", metric_names))
  
  return(results_frame)
}

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

#matern-3/2 with variance paramter = 1 and correlation length scale = 1
sim_cov <- matern_1.5_cov(sim_dist)
sqrt_cov <- sqrt_eigen_cov(sim_cov)

loc_list <- get_loc_list(lat_lon_grid)

# see documentation of run_sim for more details.
# iters controls the number of independent spatial fields generated
# n controls number of independent samples per spatial field.
# delta gives the max value for rho * beta under the interval null
# beta is the true linear relationship between Y and the random location
# noise sd is the amount of random normal noise added to the random location

results_out1 <- run_sim_laplace(iter = 10,
                        lat_lon_grid = lat_lon_grid,
                        t_ecdf = t_ecdf,
                        szon_dist = szon_dist,
                        smer_dist = smer_dist,
                        loc_list = loc_list,
                        cov_mat = sim_cov,
                        sqrt_cov = sqrt_cov,
                        n = 50,
                        delta = 1e-3,
                        beta = 1,
                        noise_scale_param = 1,
                        verbose = TRUE)

fwrite(results_out1, file = "data/iso_sim_study_n50_laplace_sigma1.csv")
