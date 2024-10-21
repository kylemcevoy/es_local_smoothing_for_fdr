###### BH smooth illustration figure #####

library(sf)
library(ggplot2)
library(maps)
library(data.table)

source("code/functions/dist_cov.R")
source("code/functions/process_data.R")
source("code/functions/inference.R")
source("code/functions/cov_smooth.R")
source("code/functions/blockwise_MLE.R")

set.seed(41501501)

axis_labels_theme <- theme(axis.title.x = element_text(size = 12),
                            axis.text.x = element_text(size = 11),
                            axis.title.y = element_text(size = 12),
                            axis.text.y = element_text(size = 11),
                            legend.title = element_text(size = 10),
                            legend.text = element_text(size = 11),
                            legend.key.height = unit(1.4, "cm")) 

sst_sub_mle <- fread("data/sst_jan_50yrsub_mle_fit.csv")

sst_jan <- load_obs_data(path = "data/sst_midlat_anom_jan.csv",
                         max_const_prop = 0.25)


sst_jan_coord_dt <- extract_coord_dt(sst_jan)

param_dt <- cbind(sst_jan_coord_dt, sst_sub_mle)

names(param_dt)[c(5, 6)] <- c("theta", "psi")

param_dt$phi_psi <- param_dt$phi * sqrt(param_dt$psi)

land_map <- sf::st_as_sf(map(wrap = c(0, 360), plot = FALSE, fill = TRUE))

sst_jan_loc_list <- get_loc_list(sst_jan_coord_dt)

sst1_zonal_sdist <- find_signed_zonal_dist(sst_jan_coord_dt)
sst1_merid_sdist <- find_signed_merid_dist(sst_jan_coord_dt)

sst_jan_subset <- sst_jan[year > 1920 & year <= 1970]

sst_mat_sub <- extract_data_matrix(sst_jan_subset, var = "sst", scale = TRUE)

mle_param_mat <- data.matrix(param_dt[, 4:6])

sst_cov_weights <- find_smooth_weights(loc_list = sst_jan_loc_list,
                                       zonal_sdist = sst1_zonal_sdist,
                                       merid_sdist = sst1_merid_sdist,
                                       cov_func = matern_1.5_cov,
                                       mle_params = mle_param_mat,
                                       progress = FALSE)

sst_jan_sub_smooth <- sst_mat_sub %*% sst_cov_weights

y <- sst_mat_sub[, 3453] + rnorm(50)

base_tstat <- fit_lm_tstat(y, data_mat = sst_mat_sub)
smooth_tstat <- fit_lm_tstat(y, data_mat = sst_jan_sub_smooth)

t_sample <- fread("data/t_sample_delta_1e-3_n50.csv")$t_stat
t_ecdf <- ecdf(t_sample)

base_pvals <- 1 - t_ecdf(base_tstat)
smooth_pvals <- 1 - t_ecdf(smooth_tstat)

BH_base <- find_BH_reject(base_pvals, alpha = 0.1)
BH_smooth <- find_BH_reject(smooth_pvals, alpha = 0.1)

BH_data_base <- cbind(sst_jan_coord_dt, 
                      reject = BH_base$BH_reject,
                      t_stat = base_tstat,
                      method = "base")

BH_data_smooth <- cbind(sst_jan_coord_dt,
                        reject = BH_smooth$BH_reject,
                        t_stat = smooth_tstat,
                        method = "smooth")

BH_data <- rbind(BH_data_base, BH_data_smooth)

sst_BH_example_plot <- ggplot(data = land_map) + 
  geom_tile(data = BH_data,
            aes(x = lon,
                y = lat,
                fill = t_stat)) + 
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  scale_fill_continuous(type = 'viridis') +
  geom_point(data = BH_data[reject == TRUE],
             aes(x = lon,
                 y = lat),
             color = "red",
             alpha = 0.3,
             shape = 3) +
  geom_point(data = BH_data[c(3453,
                              8146 + 3453), ],
             aes(x = lon,
                 y = lat),
             color = "red",
             alpha = 0.5) + 
  labs(x = "longitude", y = "latitude",
       fill = "t stat.") + 
  facet_wrap(vars(method), nrow = 2) + 
  theme_bw() + 
  axis_labels_theme

sst_BH_example_plot 

ggsave("images/sst_BH_example_plot.png",
       plot = sst_BH_example_plot,
       units = "in",
       width = 8,
       height = 5.5)

sum(BH_data_base$reject)

sum(BH_data_smooth$reject)  
