####### Selection Stability #######
library(data.table)
library(viridis)
library(colorspace)
library(ggplot2)

source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")

set.seed(3519350)

find_chull <- function(coord_dt, condition) {
  map_back <- which(condition)
  convex_hull <- chull(coord_dt[condition == TRUE, ])
  convex_hull <- c(convex_hull, convex_hull[1])
  orig_indices <- map_back[convex_hull]
  chull_orig_data <- coord_dt[orig_indices, ]
  
  return(chull_orig_data)
}

axis_labels_theme <- theme(axis.title.x = element_text(size = 13),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 13),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 12),
                           legend.text = element_text(size = 11),
                           legend.key.height = unit(1, "cm"),
                           legend.key.width = unit(0.75, "cm")) 

n <- 50

t_sample <- fread("data/t_sample_delta_1e-3_n50.csv")$t_stat
t_ecdf <- ecdf(t_sample)

# create lat, lon grid on which to simulate data. Choose lat. and lon. ranges
# and choose units to create resolution of grid.
lat_lon_grid <- build_lat_lon_grid(lat_range = c(-60, 60),
                                   lon_range = c(0, 358),
                                   lat_unit = 2,
                                   lon_unit = 2)

p <- nrow(lat_lon_grid)

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

tmp_MLE <- data.frame(phi = rep(0, p),
                      ani_angle = rep(0, p),
                      ani_ratio = rep(0, p),
                      convergence = rep(0, p))

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

sim_data <- gen_data_eigen(n = n, sqrt_cov)

sim_data_std <- scale(sim_data)

for (j in 1:p) {
  
  tmp_MLE[j, ] <- find_local_MLE(loc = j,
                                 coord_dt = lat_lon_grid,
                                 zonal_dist_mat = szon_dist,
                                 merid_dist_mat = smer_dist,
                                 data_mat_std = sim_data_std,
                                 return_best = TRUE,
                                 init_cond = init_cond_mat,
                                 return_likelihood = FALSE)
  
  if (j %% 1000 == 0) print(j)
  
}

MLE_mat <- data.matrix(tmp_MLE[, 1:3])

cov_weights <- find_smooth_weights(loc_list = loc_list,
                                   zonal_sdist = szon_dist,
                                   merid_sdist = smer_dist,
                                   cov_func = matern_1.5_cov,
                                   mle_params = MLE_mat,
                                   progress = FALSE)

smooth_data <- sim_data_std %*% cov_weights

base_reject_loc <- matrix(NA, nrow = 1000, ncol = 10980)
smooth_reject_loc <- matrix(NA, nrow = 1000, ncol = 10980)

for (i in 1:1000) {
  if (i %% 100 == 0) message(i)
  
  tmp_y <- sim_data[, 5475] + rnorm(n = 50)
  base_tstat <- fit_lm_tstat(tmp_y, data_mat = sim_data_std)
  smooth_tstat <- fit_lm_tstat(tmp_y, data_mat = smooth_data)
  
  base_pvals <- 1 - t_ecdf(base_tstat)
  smooth_pvals <- 1 - t_ecdf(smooth_tstat)
  
  BH_base <- find_BH_reject(base_pvals, alpha = 0.1)
  BH_smooth <- find_BH_reject(smooth_pvals, alpha = 0.1)
  
  base_reject_loc[i, ] <- BH_base$BH_reject
  smooth_reject_loc[i, ] <- BH_smooth$BH_reject
}

base_reject_freq <- apply(base_reject_loc, 2, mean)
smooth_reject_freq <- apply(smooth_reject_loc, 2, mean)
reject_freq_diff <- smooth_reject_freq - base_reject_freq

true_cov <- sim_cov[5475, ]
beta <- 1
test_chull <- find_chull(lat_lon_grid, !(true_cov * beta < 1e-3))
test_chull2 <- find_chull(lat_lon_grid, !(true_cov * beta < 1e-2))
test_chull3 <- find_chull(lat_lon_grid, !(true_cov * beta < 1e-1))

test_chull$corr <- 1e-3
test_chull2$corr <- 1e-2
test_chull3$corr <- 1e-3

test_chull_comb <- rbind(test_chull, test_chull2, test_chull3)

selection_stability_plot <- ggplot(data = cbind(lat_lon_grid,
                                            reject_freq_diff = reject_freq_diff),
       aes(x = lon, y = lat)) + 
  geom_tile(aes(fill = reject_freq_diff)) +
  geom_point(data = lat_lon_grid[5475, ], color = "red", alpha = 0.3) +
  scale_fill_stepsn(colors = rev(diverging_hcl(20,
                                               h = c(180, 50),
                                               c = 80,
                                               l = c(20, 95),
                                               power = c(0.7, 1.3))),
                    n.breaks = 10,
                    limits = c(-0.5, 0.5),
                    name = "reject freq. diff.") +
  geom_path(data = test_chull, color = "red") +
  geom_path(data = test_chull2, color = "purple") +
  geom_path(data = test_chull3, color = "blue")+ 
  coord_sf(xlim = c(100, 196), ylim = c(-45, 45)) + 
  labs(x = "longitude", y = "latitude") +
  theme_bw() +
  axis_labels_theme

selection_stability_plot

ggsave('images/selection_stability_plot.png', selection_stability_plot)

### Results: Top 10 average rejection freq and top 50 average rejection freq.

mean(sort(base_reject_freq, decreasing = TRUE)[1:10])
mean(sort(smooth_reject_freq, decreasing = TRUE)[1:10])

mean(sort(base_reject_freq, decreasing = TRUE)[1:50])
mean(sort(smooth_reject_freq, decreasing = TRUE)[1:50])
  
