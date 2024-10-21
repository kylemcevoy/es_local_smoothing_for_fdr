library(sf)
library(tidyr)
library(ggplot2)
library(maps)
library(pals)
library(scales)
library(data.table)
library(gridExtra)
library(gtable)
library(grid)

source("code/functions/dist_cov.R")
source("code/functions/process_data.R")
source("code/functions/inference.R")
source("code/functions/cov_smooth.R")
source("code/functions/blockwise_MLE.R")

## Function from StilRiv answer at https://stackoverflow.com/questions/26159495/align-multiple-ggplot-graphs-with-and-without-legends?noredirect=1&lq=1
## for aligning multiple ggplots.

AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}

axis_labels_theme <- theme(axis.title.x = element_text(size = 12),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 12),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 14),
                           legend.text = element_text(size = 11),
                           legend.key.height = unit(0.6, "cm"))

axis_labels_theme_phi <- theme(axis.title.x = element_text(size = 12),
                               axis.text.x = element_text(size = 11),
                               axis.title.y = element_text(size = 12),
                               axis.text.y = element_text(size = 11),
                               legend.title = element_text(size = 14),
                               legend.text = element_text(size = 11),
                               legend.key.height = unit(0.62, "cm"))

sst_sub_mle <- fread("data/sst_jan_50yrsub_mle_fit.csv")

sst_jan <- load_obs_data(path = "data/sst_midlat_anom_jan.csv",
                         max_const_prop = 0.25)


sst_jan_coord_dt <- extract_coord_dt(sst_jan)

param_dt <- cbind(sst_jan_coord_dt, sst_sub_mle)

names(param_dt)[c(5, 6)] <- c("theta", "psi")

param_dt$phi_psi <- param_dt$phi * sqrt(param_dt$psi)

land_map <- sf::st_as_sf(map(wrap = c(0, 360), plot = FALSE, fill = TRUE))

phi_plot <- ggplot(data = land_map) +
  geom_tile(data = param_dt, aes(x = lon, y = lat, fill = phi)) +
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  scale_fill_binned(type = "viridis",
                    name = expression(paste(~~phi, " [", 1 %.% 10^3, " km]")),
                    breaks = seq(0, 3, by = 1/2),
                    limits = c(0, 3)
  ) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  labs(x = NULL, y = "latitude") + 
  theme_bw()

phi_plot

phase_palette <- ocean.phase(5)

theta_plot <- ggplot(data = land_map) +
  geom_tile(data = param_dt, aes(x = lon, y = lat, fill = theta)) +
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  scale_fill_stepsn(colors = phase_palette,
                    name = expression(paste(theta, " [radians]")),
                    breaks = seq(0, pi, by = pi / 5),
                    labels = as.list(c(0,
                               expression(pi / 5), 
                               expression(2 * pi / 5),
                               expression(3 * pi / 5),
                               expression(4 * pi / 5)))
  ) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  labs(x = NULL, y = "latitude") + 
  theme_bw()

theta_plot

phipsi_plot <- ggplot(data = land_map) +
  geom_tile(data = param_dt, aes(x = lon, y = lat, fill = phi_psi)) +
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  scale_fill_binned(type = "viridis",
                    name = expression(paste(phi * sqrt(psi), " [", 1 %.% 10^3, " km]")),
                    breaks = seq(0, 5, by = 1),
                    limits = c(0, 5)
  )  + 
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  labs(x = NULL, y = "latitude")

phipsi_plot

aligned_plots <- AlignPlots(phi_plot, theta_plot, phipsi_plot)

gridded_mle <- grid.arrange(aligned_plots[[1]], aligned_plots[[2]], aligned_plots[[3]],
                            nrow = 3,
                            bottom = "longitude")

ggsave("images/stacked_sst_mle2.png",
       gridded_mle,
       units = 'in',
       width = 8,
       height = 8)


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

sst_smooth_plot_dt <- cbind(sst_jan_coord_dt, 
                        unsmooth = sst_mat_sub[33, ],
                        smooth = sst_jan_sub_smooth[33, ])


sst_smooth_plot_dt_pl <- pivot_longer(sst_smooth_plot_dt,
                                      cols = c(unsmooth, smooth),
                                      values_to = "anomaly",
                                      names_to = "smooth")

sst_smooth_plot_dt_pl$smooth <- factor(sst_smooth_plot_dt_pl$smooth,
                                       levels = c("unsmooth", "smooth"))

axis_labels_theme2 <- theme(axis.title.x = element_text(size = 12),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 12),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 10),
                           legend.text = element_text(size = 11),
                           legend.key.height = unit(1.4, "cm")) 

smooth_vs_unsmooth_ex_plot <- ggplot(data = land_map) +
  geom_tile(data = sst_smooth_plot_dt_pl,
            aes(x = lon,
                y = lat, 
                fill = anomaly)) +
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  scale_fill_continuous(name = "std. anom. [ ]",
                        type = "viridis") +
  facet_wrap(vars(smooth), nrow = 2) +
  theme_bw() +
  labs(x = "longitude", y = "latitude") + 
  axis_labels_theme2

smooth_vs_unsmooth_ex_plot

ggsave("images/sst_1953_example2.png",
       plot = smooth_vs_unsmooth_ex_plot,
       units = "in",
       width = 8, height = 5.5)
