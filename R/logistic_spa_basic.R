# Modified from https://github.com/leeshawn/SPAtest

#' @title Basic SPA in logistic regression
#'
#' @description Basic saddlepoint approximation method in logistic regression model.
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
#' @references Dey, Rounak, et al. "A fast and accurate algorithm to test for binary phenotypes and its application to PheWAS." The American Journal of Human Genetics 101.1 (2017): 37-49.
logistic_spa <- function(genos, pheno, cov, obj.null, minmac = 5, cutoff = 2, log.p = FALSE) {
  if (missing(obj.null)) {
    if (missing(cov) || is.null(cov)) {
      cov <- rep(1, length(pheno))
    }
    obj.null <- logistic_fit_null(as.matrix(pheno) ~ as.matrix(cov))
  }
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
  Is.converge <- rep(NA, m)
  for (i in 1:m)
  {
    try({
      ina <- which(is.na(genos[i, ]))
      if (length(ina) > 0) {
        genos[i, ina] <- mean(genos[i, ], na.rm = TRUE)
      }
      MAC <- min(sum(genos[i, ]), sum(2 - genos[i, ]))
      if (MAC >= minmac) {
        re <- logistic_score_test(as.vector(genos[i, , drop = FALSE]), obj.null, cutoff = cutoff, log.p = log.p)

        p.value[i] <- re$p.value
        p.value.NA[i] <- re$p.value.NA
        Is.converge[i] <- re$Is.converge
      }
    })
    if (i %% 1000 == 0) print(paste("Processed", i, "SNPs", sep = " "))
  }

  return(list(p.value = p.value, p.value.NA = p.value.NA, Is.converge = Is.converge))
}

# Fit null model

logistic_X1_adj <- function(X1) {
  q1 <- ncol(X1)
  if (q1 >= 2) {
    if (sum(abs(X1[, 1] - X1[, 2])) == 0) {
      X1 <- X1[, -2]
      q1 <- q1 - 1
    }
  }
  X1.qr <- qr(X1)
  if (X1.qr$rank < q1) {
    X1.svd <- svd(X1)
    X1 <- X1.svd$u[, 1:X1.qr$rank]
  }
  return(X1)
}

#' @importFrom stats model.matrix
#' @importFrom stats glm
logistic_fit_null <- function(formula, data = NULL) {
  X1 <- model.matrix(formula, data = data)
  X1 <- logistic_X1_adj(X1)

  glmfit <- glm(formula, data = data, family = "binomial")
  mu <- glmfit$fitted.values

  V <- mu * (1 - mu)
  res <- glmfit$y - mu
  XV <- t(X1 * V)
  XVX_inv <- solve(t(X1) %*% (X1 * V))
  XXVX_inv <- X1 %*% XVX_inv

  re <- list(y = glmfit$y, mu = mu, res = res, V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv)
  class(re) <- "logistic_null"

  return(re)
}

# Compute CGF function

logistic_K0 <- function(t, mu, g) {
  n.t <- length(t)
  out <- rep(0, n.t)

  for (i in 1:n.t) {
    t1 <- t[i]
    temp <- log(1 - mu + mu * exp(g * t1))
    out[i] <- sum(temp)
  }

  return(out)
}

logistic_K1 <- function(t, mu, g, q) {
  n.t <- length(t)
  out <- rep(0, n.t)

  for (i in 1:n.t) {
    t1 <- t[i]
    temp1 <- (1 - mu) * exp(-g * t1) + mu
    temp2 <- mu * g
    out[i] <- sum(temp2 / temp1) - q
  }

  return(out)
}

logistic_K2 <- function(t, mu, g) {
  n.t <- length(t)
  out <- rep(0, n.t)

  for (i in 1:n.t) {
    t1 <- t[i]
    temp1 <- ((1 - mu) * exp(-g * t1) + mu)^2
    temp2 <- (1 - mu) * mu * g^2 * exp(-g * t1)
    out[i] <- sum(temp2 / temp1, na.rm = TRUE)
  }

  return(out)
}

logistic_get_root <- function(init, mu, g, q, tol = .Machine$double.eps^0.25, maxiter = 1000) {
  g.pos <- sum(g[which(g > 0)])
  g.neg <- sum(g[which(g < 0)])

  if (q >= g.pos || q <= g.neg) {
    return(list(root = Inf, n.iter = 0, Is.converge = TRUE))
  } else {
    t <- init
    K1_eval <- logistic_K1(t, mu, g, q)
    prevJump <- Inf
    rep <- 1
    repeat {
      K2_eval <- logistic_K2(t, mu, g)
      tnew <- t - K1_eval / K2_eval
      if (is.na(tnew)) {
        conv <- FALSE
        break
      }
      if (abs(tnew - t) < tol) {
        conv <- TRUE
        break
      }
      if (rep == maxiter) {
        conv <- FALSE
        break
      }

      newK1 <- logistic_K1(tnew, mu, g, q)
      if (sign(K1_eval) != sign(newK1)) {
        if (abs(tnew - t) > prevJump - tol) {
          tnew <- t + sign(newK1 - K1_eval) * prevJump / 2
          newK1 <- logistic_K1(tnew, mu, g, q)
          prevJump <- prevJump / 2
        } else {
          prevJump <- abs(tnew - t)
        }
      }

      rep <- rep + 1
      t <- tnew
      K1_eval <- newK1
    }
    return(list(root = t, n.iter = rep, Is.converge = conv))
  }
}

#' @importFrom stats pnorm
logistic_get_prob <- function(zeta, mu, g, q, log.p = FALSE) {
  k1 <- logistic_K0(zeta, mu, g)
  k2 <- logistic_K2(zeta, mu, g)

  if (is.finite(k1) && is.finite(k2)) {
    temp1 <- zeta * q - k1

    w <- sign(zeta) * (2 * temp1)^{
      1 / 2
    }
    v <- zeta * (k2)^{
      1 / 2
    }

    Z.test <- w + 1 / w * log(v / w)

    if (Z.test > 0) {
      pval <- pnorm(Z.test, lower.tail = FALSE, log.p = log.p)
    } else {
      pval <- -pnorm(Z.test, lower.tail = TRUE, log.p = log.p)
    }
  } else {
    if (log.p) {
      pval <- -Inf
    } else {
      pval <- 0
    }
  }

  return(pval)
}

# Compute p value

logistic_add_logp <- function(p1, p2) {
  p1 <- -abs(p1)
  p2 <- -abs(p2)
  maxp <- max(p1, p2)
  minp <- min(p1, p2)

  return(maxp + log(1 + exp(minp - maxp)))
}

#' @importFrom stats pchisq
logistic_get_pvalue <- function(q, mu, g, cutoff = 2, log.p = FALSE) {
  m1 <- sum(mu * g)
  var1 <- sum(mu * (1 - mu) * g^2)
  p1 <- NULL
  p2 <- NULL

  score <- q - m1
  qinv <- -sign(q - m1) * abs(q - m1) + m1

  # Noadj
  pval.noadj <- pchisq((q - m1)^2 / var1, lower.tail = FALSE, df = 1, log.p = log.p)
  Is.converge <- TRUE

  if (cutoff < 10^-1) cutoff <- 10^-1

  if (abs(q - m1) / sqrt(var1) < cutoff) {
    pval <- pval.noadj
  } else {
    out.uni1 <- logistic_get_root(0, mu = mu, g = g, q = q)
    out.uni2 <- logistic_get_root(0, mu = mu, g = g, q = qinv)
    if (out.uni1$Is.converge == TRUE && out.uni2$Is.converge == TRUE) {
      p1 <- tryCatch(logistic_get_prob(out.uni1$root, mu, g, q, log.p = log.p), error = function(e) {
        if (log.p) {
          return(pval.noadj - log(2))
        } else {
          return(pval.noadj / 2)
        }
      })
      p2 <- tryCatch(logistic_get_prob(out.uni2$root, mu, g, qinv, log.p = log.p), error = function(e) {
        if (log.p) {
          return(pval.noadj - log(2))
        } else {
          return(pval.noadj / 2)
        }
      })
      if (log.p) {
        pval <- logistic_add_logp(p1, p2)
      } else {
        pval <- abs(p1) + abs(p2)
      }
      Is.converge <- TRUE
    } else {
      print("Error_Converge")
      pval <- pval.noadj
      Is.converge <- FALSE
    }
  }

  if (pval != 0 && pval.noadj / pval > 10^3) {
    return(logistic_get_pvalue(q, mu, g, cutoff = cutoff * 2, log.p = log.p))
  } else {
    return(list(p.value = pval, p.value.NA = pval.noadj, Is.converge = Is.converge, score = score))
  }
}

logistic_score_test <- function(G, obj.null, cutoff = 2, log.p = FALSE) {
  if (class(obj.null) != "logistic_null") {
    stop("obj.null should be a returned object from logistic_fit_null")
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
  out <- logistic_get_pvalue(q, mu = mu, g = G1, cutoff = cutoff, log.p = log.p)

  return(out)
}
