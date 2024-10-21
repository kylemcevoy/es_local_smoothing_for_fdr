library(tidyr)
library(ggplot2)

axis_labels_theme <- theme(axis.title.x = element_text(size = 12),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 12),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 11),
                           legend.text = element_text(size = 11),
                           legend.key.height = unit(0.67, "cm")) 

# n = 50
noise_sim_study_results <- fread("data/noise_sim_study.csv")
# n = 100
n_sim_study_results <- fread("data/n_sim_study.csv")

filter_nonzero <- function(data) {
  filter_data <- data[data$base_reject > 0 | data$smooth_reject > 0]
  filter_data$true_reject_diff <- (filter_data$smooth_true_reject -
                                     filter_data$base_true_reject)
  return(filter_data)
}

noise_sim_filtered <- filter_nonzero(noise_sim_study_results)
n_sim_filtered <- filter_nonzero(n_sim_study_results)

noise_sim_filtered$sigma <- factor(noise_sim_filtered$sigma, levels = c(1, 2, 3))
noise_sim_filtered$n <- factor(noise_sim_filtered$n, levels = c(50, 100))

n_sim_filtered$sigma <- factor(n_sim_filtered$sigma, levels = c(1, 2, 3))
n_sim_filtered$n <- factor(n_sim_filtered$n, levels = c(50, 100))

sim_filtered <- rbind(noise_sim_filtered, n_sim_filtered)

mean_frame <- sim_filtered[, .(mean = mean(true_reject_diff)), by = .(sigma, n)]
mean_frame$zero <- 0

facet_labels <- c("50" = "n = 50",
                  "100" = "n = 100")

reject_hist_plot <- ggplot(data = sim_filtered, aes(x = true_reject_diff,
                                  fill = sigma)) +
  geom_histogram(color = "black", binwidth = 10, position = "identity", alpha = 0.7) +
  geom_vline(data = mean_frame,
             aes(xintercept = zero),
                 color = "black",
             linetype = 2) +
  geom_vline(data = mean_frame,
             aes(xintercept = mean, color = sigma),
             linetype = 2) +
  labs(x = "difference in # of true rejections per sim. (smoothed - base)",
       color = expression(sigma[epsilon]),
       fill = expression(sigma[epsilon])) + 
  facet_wrap(vars(n), nrow = 2, labeller = as_labeller(facet_labels)) + 
  theme_bw() +
  axis_labels_theme

reject_hist_plot

ggsave(filename = "images/reject_hist_plot.png", plot = reject_hist_plot,
       units = "in",
       width = 8,
       height = 6)
