source("code/functions/blockwise_MLE.R")
source("code/functions/process_data.R")
source("code/functions/dist_cov.R")

# Only gridboxes with fewer than 25% of values being identical are kept.
# This mostly arises at locations that are frozen over for a portion of the
# year.

sst_jan <- load_obs_data(path = "data/sst_midlat_anom_jan.csv",
                         max_const_prop = 0.25)

sst_jan_mat_std <- extract_data_matrix(sst_jan, "sst", scale = TRUE)

sst_jan_coord_dt <- extract_coord_dt(sst_jan)

sst1_zonal_sdist <- find_signed_zonal_dist(sst_jan_coord_dt)
sst1_merid_sdist <- find_signed_merid_dist(sst_jan_coord_dt)

MLE_out <- data.frame(phi = rep(0, 8146),
                      ani_angle = rep(0, 8146),
                      ani_ratio = rep(0, 8146),
                      convergence = rep(0, 8146))

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

for (i in 1:8146) {
  MLE_out[i, ] <- find_local_MLE(loc = i,
                                 coord_dt = sst_jan_coord_dt,
                                 zonal_dist_mat = sst1_zonal_sdist,
                                 merid_dist_mat = sst1_merid_sdist,
                                 data_mat_std = sst_jan_mat_std,
                                 init_cond = init_cond_mat,
                                 return_best = TRUE,
                                 return_likelihood = FALSE)
  
}

fwrite(MLE_out, file = "data/sst_jan_mle_fit.csv")

# 6 locations have convergence issues in at least one init. cond., but not in the
# init. cond. which had the max likelihood. In some locations the best location has an
# abnormal termination of line search convergence status in L-BFGS, but the second
# best likelihood is within a relative tolerance of the max likelihood and has
# convergence. 

# Subset to 50 years

sst_jan_subset <- sst_jan[year > 1920 & year <= 1970]

sst_mat_sub <- extract_data_matrix(sst_jan_subset, var = "sst", scale = TRUE)

MLE_sub <- data.frame(phi = rep(0, 8146),
                      ani_angle = rep(0, 8146),
                      ani_ratio = rep(0, 8146),
                      convergence = rep(0, 8146))

for (i in 1:8146) {
  MLE_sub[i, ] <- find_local_MLE(loc = i,
                                 coord_dt = sst_jan_coord_dt,
                                 zonal_dist_mat = sst1_zonal_sdist,
                                 merid_dist_mat = sst1_merid_sdist,
                                 data_mat_std = sst_mat_sub,
                                 init_cond = init_cond_mat,
                                 return_best = TRUE,
                                 return_likelihood = FALSE)
  
  
}

# 9 initial condition warnings: 455, 532, 533, 1379, 1507, 1550, 1551, 1553, 6838

fwrite(MLE_sub, file = "data/sst_jan_50yrsub_mle_fit.csv")
