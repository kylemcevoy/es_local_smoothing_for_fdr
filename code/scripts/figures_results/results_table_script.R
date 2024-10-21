## Reported values in results section

noise_sim_study_mean <- fread("data/noise_sim_study_mean.csv")
n_sim_study_mean <- fread("data/n_sim_study_mean.csv")
phi_sim_study_mean <- fread("data/phi_sim_study_mean.csv")
aniso_sim_study_mean <- fread("data/aniso_sim_study_mean.csv")

results_table <- rbind(noise_sim_study_mean,
                       n_sim_study_mean,
                       phi_sim_study_mean,
                       aniso_sim_study_mean)

results_table[, FNP := NULL]
results_table[, reject := NULL]

results_table <- results_table[, .(method, sigma, n, phi, theta, psi, FDP, FDX, true_reject, sens, spec)]

results_table

sens_table <- results_table[, .(method, sigma, n, phi, theta, psi, sens)]

wide_sens_table <- dcast(sens_table,
                         sigma + n + phi + theta + psi ~ method,
                         value.var = "sens")

wide_sens_table[, perc_delta_sens := 100 * ((smooth - base) / base)]

wide_sens_table[order(theta, psi, phi, n, sigma)]


