# Results for Table 2 of paper on misspecification.

library(data.table)

matern_05_results_out <- fread("data/iso_sim_study_n50_sigma1_matern0_5.csv")
matern_25_results_out <- fread("data/iso_sim_study_n50_sigma1_matern2_5.csv")
laplace_results_out <- fread("data/iso_sim_study_n50_laplace_sigma1.csv")

# Matern 1/2

## proportion of time base methods rejects something but smooth method doesn't
sum(matern_05_results_out$smooth_reject == 0 & matern_05_results_out$base_reject > 0) / 10000


## mean difference in number of true rejections per simulation
mean(matern_05_results_out$base_true_reject)
mean(matern_05_results_out$smooth_true_reject)

# mean FDP for the base and smooth methods
mean(matern_05_results_out$base_FDP)
mean(matern_05_results_out$smooth_FDP)

# mean sensitivity
mean(matern_05_results_out$base_sens)
mean(matern_05_results_out$smooth_sens)
#proportional increase
100 * (mean(matern_05_results_out$smooth_sens) / mean(matern_05_results_out$base_sens) - 1)

# mean specificity
mean(matern_05_results_out$base_spec)
mean(matern_05_results_out$smooth_spec)

mean(matern_05_results_out$base_FDP > 0.1)
mean(matern_05_results_out$smooth_FDP > 0.1)

# Matern 5/2

sum(matern_25_results_out$smooth_reject == 0 & matern_25_results_out$base_reject > 0) / 10000


## mean difference in number of true rejections per simulation
mean(matern_25_results_out$base_true_reject)
mean(matern_25_results_out$smooth_true_reject)

# mean FDP for the base and smooth methods
mean(matern_25_results_out$base_FDP)
mean(matern_25_results_out$smooth_FDP)

# mean sensitivity
mean(matern_25_results_out$base_sens)
mean(matern_25_results_out$smooth_sens)
#proportional increase
100 * (mean(matern_25_results_out$smooth_sens) / mean(matern_25_results_out$base_sens) - 1)

# mean specificity
mean(matern_25_results_out$base_spec)
mean(matern_25_results_out$smooth_spec)

# FDX
mean(matern_25_results_out$base_FDP > 0.1)
mean(matern_25_results_out$smooth_FDP > 0.1)


## Laplace Noise

## mean difference in number of true rejections per simulation
mean(laplace_results_out$base_true_reject)
mean(laplace_results_out$smooth_true_reject)

# mean FDP for the base and smooth methods
mean(laplace_results_out$base_FDP)
mean(laplace_results_out$smooth_FDP)

# mean sensitivity
mean(laplace_results_out$base_sens)
mean(laplace_results_out$smooth_sens)
#proportional increase
100 * (mean(laplace_results_out$smooth_sens) / mean(laplace_results_out$base_sens) - 1)

# mean specificity
mean(laplace_results_out$base_spec)
mean(laplace_results_out$smooth_spec)

# FDX
mean(laplace_results_out$base_FDP > 0.1)
mean(laplace_results_out$smooth_FDP > 0.1)

