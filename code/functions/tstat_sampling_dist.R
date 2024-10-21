# t statistic monte carlo estimated sampling distribution generation. The 
# sampling distribution is used with ecdf() to generate p-values.

#' Cholesky Decomp. of 2 x 2 symmetric matrix.
#' 
#' Finds the cholesky decomposition of a 2 x 2 symmetric correlation matrix
#' with correlation parameter rho.
#' @param rho the correlation between the two variables
#' @returns the 2 x 2 upper triangular cholesky decomposition of the matrix.
delta_chol <- function(rho) {
  chol(matrix(c(1, rho, rho, 1), nrow = 2))
}

#' Monte carlo draw of the t statistic.
#' 
#' Draw a simulated t statistic under the specified null.
#' See the paper for more details.
#' @param chol_S the upper triangular cholesky decomposition of the desired
#' correlation matrix
#' @param n the sample size of the multivariate normal
#' @param beta the true beta param between Y and X1.
#' @returns a t statistic sampled under the given null.
draw_tstat <- function(chol_S, n, beta) {
  Z <- matrix(rnorm(2 * n), ncol = 2)
  X <- Z %*% chol_S
  Y <- beta * X[, 1] + rnorm(n)

  fitted_lm <- lm.fit(y = Y, x = cbind(1, X[, 2]))
  beta_x <- unname(fitted_lm$coefficients[2])
  MSE <- sum(fitted_lm$residuals^2) / (n - 2)
  SS_x <- (n - 1) * var(X[, 2])
  se_beta <- sqrt(MSE / SS_x)
  t_stat <- beta_x / se_beta

  return(t_stat)
}
