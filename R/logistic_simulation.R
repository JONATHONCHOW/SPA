#' @title Simulations in logistic regression
#'
#' @description Numerical simulations in logistic regression model.
#'
#' @param cases
#' @param controls
#' @param maf
#' @param X1prob1_case
#' @param X1prob1_control
#' @param samplerX2_given_Y0
#' @param samplerX2_given_Y1
#'
#' @return A list.
#'
#' @export
#'
#' @examples
#' @references Dey, Rounak, et al. "A fast and accurate algorithm to test for binary phenotypes and its application to PheWAS." The American Journal of Human Genetics 101.1 (2017): 37-49.
logistic_simulation <- function(cases, controls, maf, X1prob1_case = 0.7260159, X1prob1_control = 0.4974979, samplerX2_given_Y0 = arscpp::ars(f = Log_Dens_X2_given_Y0, f_prime = Log_Dens_X2_given_Y0_prime, xlb = -Inf, xrb = Inf, x = c(-2, 0, 2), beta0 = -5.6, prevalence = 0.0109492), samplerX2_given_Y1 = arscpp::ars(f = Log_Dens_X2_given_Y1, f_prime = Log_Dens_X2_given_Y1_prime, xlb = -Inf, xrb = Inf, x = c(-2, 0, 2), beta0 = -5.6, prevalence = 0.0109492)) {
  N <- cases + controls
  ccratio <- cases / N

  genos <- rbinom(N, 2, maf)
  pheno <- c(rep(0, controls), rep(1, cases))
  X1 <- c(rbinom(n = controls, size = 1, prob = X1prob1_control), rbinom(n = cases, size = 1, prob = X1prob1_case))
  X2 <- c(samplerX2_given_Y0$sample(controls), samplerX2_given_Y1$sample(cases))
  cov <- as.matrix(data.frame(intercept = rep(1, N), X1 = X1, X2 = X2))
  return(list(ccratio = ccratio, genos = genos, pheno = pheno, cov = cov))
}

# Log density of X2 given Y1 = 1
Log_Dens_X2_given_Y1 <- function(x, beta0, prevalence) {
  return(log((0.5 * (1 / (1 + exp(-beta0 - x))) * (sqrt(1 / (2 * pi))) * exp(-0.5 * x^2) + 0.5 * (1 / (1 + exp(-beta0 - 1 - x))) * (sqrt(1 / (2 * pi))) * exp(-0.5 * x^2)) / (prevalence)))
}
Log_Dens_X2_given_Y1_prime <- function(x, beta0, prevalence) {
  return(-x + (1 / (1 / (1 + exp(-beta0 - x)) + 1 / (1 + exp(-beta0 - x - 1)))) * (exp(-beta0 - x) / (1 + exp(-beta0 - x))^2 + exp(-beta0 - 1 - x) / (1 + exp(-beta0 - 1 - x))^2))
}

# Log density of X2 given Y1 = 0
Log_Dens_X2_given_Y0 <- function(x, beta0, prevalence) {
  return(log((dnorm(x) - 0.5 * (1 / (1 + exp(-beta0 - x))) * (sqrt(1 / (2 * pi))) * exp(-0.5 * x^2) - 0.5 * (1 / (1 + exp(-beta0 - 1 - x))) * (sqrt(1 / (2 * pi))) * exp(-0.5 * x^2)) / (1 - prevalence)))
}
Log_Dens_X2_given_Y0_prime <- function(x, beta0, prevalence) {
  return(-x + (1 / (1 - 0.5 / (1 + exp(-beta0 - x)) - 0.5 / (1 + exp(-beta0 - x - 1)))) * (-0.5 * exp(-beta0 - x) / (1 + exp(-beta0 - x))^2 - 0.5 * exp(-beta0 - 1 - x) / (1 + exp(-beta0 - x - 1))^2))
}
