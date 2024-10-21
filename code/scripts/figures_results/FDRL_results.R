# Results the FDR_L method of Zhang et al. 2011. See paper for full citation.

FDRL_small_sigma1 <- fread('data/fdrl_small_iso_sim_study_n50_sigma1.csv')
FDRL_small_sigma2 <- fread('data/fdrl_small_iso_sim_study_n50_sigma2.csv')
FDRL_small_sigma3 <- fread('data/fdrl_small_iso_sim_study_n50_sigma3.csv')

FDRL_big_sigma1 <- fread('data/fdrl_iso_sim_study_n50_sigma1.csv')
FDRL_big_sigma2 <- fread('data/fdrl_iso_sim_study_n50_sigma2.csv')
FDRL_big_sigma3 <- fread('data/fdrl_iso_sim_study_n50_sigma3.csv')

base_sigma1 <- fread('data/iso_sim_study_n50_sigma1.csv')
base_sigma2 <- fread('data/iso_sim_study_n50_sigma2.csv')
base_sigma3 <- fread('data/iso_sim_study_n50_sigma3.csv')

# results sigma = 1

mean_frame_base1 <- apply(base_sigma1, 2, mean)
mean_frame_FDRL_small1 <- apply(FDRL_small_sigma1, 2, mean)
mean_frame_FDRL_big1 <- apply(FDRL_big_sigma1, 2, mean)

mean_frame_FDRL_small1
(mean_frame_FDRL_small1[6] - mean_frame_base1[6]) / mean_frame_base1[6]

mean_frame_FDRL_big1
(mean_frame_FDRL_big1[6] - mean_frame_base1[6]) / mean_frame_base1[6]


# results sigma = 2

mean_frame_base2 <- apply(base_sigma2, 2, mean)
mean_frame_FDRL_small2 <- apply(FDRL_small_sigma2, 2, mean)
mean_frame_FDRL_big2 <- apply(FDRL_big_sigma2, 2, mean)

mean_frame_FDRL_small2
(mean_frame_FDRL_small2[6] - mean_frame_base2[6]) / mean_frame_base2[6]

mean_frame_FDRL_big2
(mean_frame_FDRL_big2[6] - mean_frame_base2[6]) / mean_frame_base2[6]

# results sigma = 3

mean_frame_base3 <- apply(base_sigma3, 2, mean)
mean_frame_FDRL_small3 <- apply(FDRL_small_sigma3, 2, mean)
mean_frame_FDRL_big3 <- apply(FDRL_big_sigma3, 2, mean)

mean_frame_FDRL_small3
(mean_frame_FDRL_small3[6] - mean_frame_base3[6]) / mean_frame_base3[6]

mean_frame_FDRL_big3
(mean_frame_FDRL_big3[6] - mean_frame_base3[6]) / mean_frame_base3[6]

