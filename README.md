Code for implementing the local smoothing method proposed in "Improving the Power of Multiple Hypothesis Testing in Large Environmental Datasets Using Spatial Smoothing" by Kyle R. McEvoy and Karen A. McKinnon, and code for generating the simulation results and figures presented in that paper.

All code for results and figures in the paper was written and executed in the R programming language using a 2021 M1 Macbook Pro. macOS versions 14.4-14.6 were used with the Apple Accelerate framework BLAS linked. R session and package info below. 

R sessionInfo():
R version 4.3.0 (2023-04-21)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 14.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] data.table_1.15.99

loaded via a namespace (and not attached):
 [1] utf8_1.2.3         R6_2.5.1           tidyselect_1.2.0   e1071_1.7-13      
 [5] magrittr_2.0.3     glue_1.6.2         tibble_3.2.1       KernSmooth_2.23-21
 [9] pkgconfig_2.0.3    generics_0.1.3     dplyr_1.1.2        lifecycle_1.0.3   
[13] classInt_0.4-9     sf_1.0-14          cli_3.6.1          fansi_1.0.4       
[17] grid_4.3.0         vctrs_0.6.2        DBI_1.1.3          proxy_0.4-27      
[21] class_7.3-22       compiler_4.3.0     rstudioapi_0.14    tools_4.3.0       
[25] pillar_1.9.0       Rcpp_1.0.10        rlang_1.1.1        units_0.8-3 
