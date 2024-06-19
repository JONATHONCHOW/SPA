# Modified from https://github.com/saigegit/SAIGE

#' @title Fast SPA in logistic regression
#'
#' @description Fast saddlepoint approximation method in logistic regression model.
#'
#' @param genos A vector or matrix containing the genotypes. If matrix is provided then rows should correspond to SNPs and columns should correspond to subjects. Optional, but needed if \code{obj.null} is missing.
#' @param pheno A vector containing the phenotypes. Optional, but needed if \code{obj.null} is missing.
#' @param cov A matrix or data frame containing the covariates. Optional, but needed if \code{obj.null} is missing.
#' @param obj.null An object of class "\code{logistic_null}". (Optional)
#' @param minmac An integer denoting the minimum minor allele count threshold to run SPA test. Default value is 5.
#' @param cutoff An integer denoting the standard deviation cutoff to be used. If the test statistic lies within the standard deviation cutoff of the mean, p-value based on traditional score test is returned. Default value is 2.
#' @param log.p A Boolean value determines whether to return natural log-transformed p-values. Default value is \code{FALSE}.
#'
#' @return A list of p values.
#'
#' @export
#'
#' @examples
#' # See vignettes
#' @references Zhou, Wei, et al. "Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies." Nature genetics 50.9 (2018): 1335-1341.
logistic_spa_fast <- function(genos, pheno, cov, obj.null, minmac = 5, cutoff = 2, log.p = FALSE) {
  if (missing(obj.null)) {
    if (missing(cov) || is.null(cov)) {
      cov <- rep(1, length(pheno))
    }
    obj.null <- fit_null(as.matrix(pheno) ~ as.matrix(cov))
  }
  # cov <- obj.null$X1
  # pheno <- obj.null$y
  cov <- obj.null$X1
  pheno <- obj.null$y

  genos <- as.matrix(genos)
  if (ncol(genos) == 1) {
    m <- 1
    genos <- t(genos)
  } else {
    m <- nrow(genos)
  }
  p.value <- rep(NA, m)
  p.value.NA <- rep(NA, m)
  # Is.converge <- rep(NA, m)
  for (i in 1:m)
  {
    try({
      ina <- which(is.na(genos[i, ]))
      if (length(ina) > 0) {
        genos[i, ina] <- mean(genos[i, ], na.rm = TRUE)
      }
      MAC <- min(sum(genos[i, ]), sum(2 - genos[i, ]))
      if (MAC >= minmac) {
        re <- score_test(as.vector(genos[i, , drop = FALSE]), obj.null, cutoff = cutoff, log.p = log.p)

        p.value[i] <- re$p.value
        p.value.NA[i] <- re$p.value.NA
        # Is.converge[i] <- re$Is.converge
      }
    })
    if (i %% 1000 == 0) print(paste("Processed", i, "SNPs", sep = " "))
  }

  return(list(p.value = p.value, p.value.NA = p.value.NA))
}

# Fit null model

#' @importFrom stats model.matrix
#' @importFrom stats glm
compute_XVX_inv <- function(XVX) {
  P_decomp <- eigen(XVX, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus <- 1 / P_decomp$values
  XVX_inv <- U %*% (Lambda_minus * t(U))
  return(XVX_inv)
}

fit_null <- function(formula, data = NULL) {
  X1 <- model.matrix(formula, data = data)
  # X1 <- X1_adj(X1)
  X1 <- X1[, -1]

  glmfit <- glm(formula, data = data, family = "binomial")
  mu <- glmfit$fitted.values

  V <- mu * (1 - mu)
  res <- glmfit$y - mu
  XV <- t(X1 * V)
  XVX_inv <- compute_XVX_inv(t(X1) %*% (X1 * V))
  XXVX_inv <- X1 %*% XVX_inv

  re <- list(y = glmfit$y, mu = mu, res = res, V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv)
  class(re) <- "logistic_null"

  return(re)
}

# Compute p value

#' @importFrom stats pchisq
#' @importFrom Rcpp evalCpp
get_pvalue <- function(q.sum, mu.a, g.sum, cutoff = 2, log.p = FALSE) {
  m1 <- sum(mu.a * g.sum)
  var1 <- sum(mu.a * (1 - mu.a) * g.sum^2)
  qinv <- -sign(q.sum - m1) * abs(q.sum - m1) + m1
  pval_noadj <- pchisq((q.sum - m1)^2 / var1, lower.tail = FALSE, df = 1, log.p = log.p)
  ## cat("pval_noadj ", pval_noadj, "\n")
  if (abs(q.sum - m1) / sqrt(var1) < cutoff) {
    p.value_burden <- pval_noadj
  } else {
    p.value_burden <- SPA_binary(mu.a, g.sum, q.sum, qinv, pval_noadj, .Machine$double.eps^0.25, 100, log.p)
  }
  # p.value_burden<-SPAtest:::Saddle_Prob(q.sum, mu=mu.a, g=g.sum, Cutoff=2,alpha=2.5*10^-6, log.p =TRUE)$p.value
  # cat("p.value_burden ", p.value_burden, "\n")
  return(list(p.value = p.value_burden, p.value.NA = pval_noadj))
}

score_test <- function(G, obj.null, cutoff = 2, log.p = FALSE) {
  if (class(obj.null) != "logistic_null") {
    stop("obj.null should be a returned object from fit_null")
  }

  y <- obj.null$y
  mu <- obj.null$mu
  res <- obj.null$res

  n.g <- sum(G)
  if (n.g / (2 * length(G)) > 0.5) {
    G <- 2 - G
    n.g <- sum(G)
  }
  G1 <- G - obj.null$XXVX_inv %*% (obj.null$XV %*% G)
  q <- sum(G1 * y)
  out <- get_pvalue(q, mu, G1, cutoff = cutoff, log.p = log.p)

  return(out)
}
