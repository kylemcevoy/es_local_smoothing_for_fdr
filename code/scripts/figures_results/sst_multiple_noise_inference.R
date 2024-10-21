library(sf)
library(ggplot2)
library(maps)
library(data.table)
library(gridExtra)

source("code/functions/dist_cov.R")
source("code/functions/process_data.R")
source("code/functions/inference.R")
source("code/functions/cov_smooth.R")
source("code/functions/blockwise_MLE.R")

set.seed(513384462)

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

t_sample <- fread("data/t_sample_delta_1e-3_n50.csv")$t_stat
t_ecdf <- ecdf(t_sample)

sim_rejections <- function(indx, reps) {
  p <- ncol(sst_mat_sub)
  output_data_list <- vector(mode = "list", length = reps)
  
  for (i in 1:reps) {
    y <- sst_mat_sub[, indx] + rnorm(50)
    
    base_tstat <- fit_lm_tstat(y, data_mat = sst_mat_sub)
    smooth_tstat <- fit_lm_tstat(y, data_mat = sst_jan_sub_smooth)
    
    base_pvals <- 1 - t_ecdf(base_tstat)
    smooth_pvals <- 1 - t_ecdf(smooth_tstat)
    
    BH_base <- find_BH_reject(base_pvals, alpha = 0.1)
    BH_smooth <- find_BH_reject(smooth_pvals, alpha = 0.1)
    
    reject_base_only <- BH_base$BH_reject & !BH_smooth$BH_reject
    reject_smooth_only <- !BH_base$BH_reject & BH_smooth$BH_reject
    reject_both <- BH_base$BH_reject & BH_smooth$BH_reject
    reject_any <- BH_base$BH_reject | BH_smooth$BH_reject
    
    reject_type_vec <- vector(mode = "character", length = p)
    reject_type_vec[reject_base_only] <- "base"
    reject_type_vec[reject_smooth_only] <- "smooth"
    reject_type_vec[reject_both] <- "both"
    reject_type_vec[!reject_any] <- "none"
    reject_type_vec <- factor(reject_type_vec, levels = c("none", "base", "smooth", "both"))
    
    BH_data <- cbind(sst_jan_coord_dt,
                     rep = i, 
                     reject_type = reject_type_vec,
                     reject_any = reject_any)
    
    output_data_list[[i]] <- BH_data
  }
  
  return(output_data_list)
}

output_data <- sim_rejections(indx = 3453, reps = 8)

BH_data <- rbindlist(output_data)

color_pal <- palette.colors(palette = "ggplot")[2:4][c(1, 3, 2)]


sst_BH_example_plot <- ggplot(data = land_map) + 
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  geom_point(data = BH_data[reject_any == TRUE],
             aes(x = lon,
                 y = lat,
                 color = reject_type),
             alpha = 0.5,
             size = 0.8,
             shape = 3) + 
  scale_x_continuous(breaks = seq(0, 340, by = 30),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  geom_point(data = BH_data[location == 3453],
             aes(x = lon, y = lat),
             color = "red") + 
  labs(x = "longitude",
       y = "latitude",
       color = "rejections") +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1")) + 
  facet_wrap(vars(rep), nrow = 4) + 
  guides(color = guide_legend(override.aes = list(size=4)))

ggsave(file = "images/plot_grid.png", sst_BH_example_plot)
