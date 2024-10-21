Code for implementing the local smoothing method proposed in "Improving the Power of Multiple Hypothesis Testing in Large Environmental Datasets Using Spatial Smoothing" by Kyle R. McEvoy and Karen A. McKinnon, and code for generating the simulation results and figures presented in that paper.

All code for results and figures in the paper was written and executed in the R programming language using a 2021 M1 Macbook Pro. macOS versions 14.4-14.6 were used with the Apple Accelerate framework BLAS linked. R session and package info below. 

R sessionInfo():
R version 4.3.0 (2023-04-21)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 14.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] colorspace_2.1-0   viridis_0.6.3      viridisLite_0.4.2  gtable_0.3.3       gridExtra_2.3     
 [6] tidyr_1.3.0        scales_1.3.0       pals_1.8           maps_3.4.1         sf_1.0-14         
[11] ggplot2_3.5.1      data.table_1.15.99

loaded via a namespace (and not attached):
 [1] dplyr_1.1.2        compiler_4.3.0     tidyselect_1.2.0   Rcpp_1.0.10        dichromat_2.0-0.1 
 [6] R6_2.5.1           generics_0.1.3     classInt_0.4-9     mapproj_1.2.11     tibble_3.2.1      
[11] units_0.8-3        munsell_0.5.0      DBI_1.1.3          pillar_1.9.0       RColorBrewer_1.1-3
[16] rlang_1.1.1        utf8_1.2.3         cli_3.6.1          withr_2.5.0        magrittr_2.0.3    
[21] class_7.3-22       rstudioapi_0.14    lifecycle_1.0.3    vctrs_0.6.2        KernSmooth_2.23-21
[26] proxy_0.4-27       glue_1.6.2         fansi_1.0.4        e1071_1.7-13       purrr_1.0.1       
[31] tools_4.3.0        pkgconfig_2.0.3   