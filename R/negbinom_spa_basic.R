#' @title Basic SPA in NB generalized linear model
#'
#' @description Basic saddlepoint approximation method in NB generalized linear model.
#'
#' @param perturbation A vector or matrix containing the perturbation indicators. Optional, but needed if \code{obj.null} is missing.
#' @param gene A vector containing the gene expressions. Optional, but needed if \code{obj.null} is missing.
#' @param cov A matrix or data frame containing the covariates. Optional, but needed if \code{obj.null} is missing.
#' @param obj.null An object of class "\code{logistic_null}". (Optional)
#' @param minperturb An integer denoting the minimum perturbation number. Default value is 5.
#' @param cutoff An integer denoting the standard deviation cutoff to be used. If the test statistic lies within the standard deviation cutoff of the mean, p-value based on traditional score test is returned. Default value is 2.
#' @param log.p A Boolean value determines whether to return natural log-transformed p-values. Default value is \code{FALSE}.
#'
#' @return A list of p values.
#'
#' @export
#'
#' @examples
#' # See vignettes
negbinom_spa <- function(perturbation, gene, cov, obj.null, minperturb = 5, cutoff = 2, log.p = FALSE) {
  if (missing(obj.null)) {
    if (missing(cov) || is.null(cov)) {
      cov <- rep(1, length(gene))
    }
    obj.null <- negbinom_fit_null(as.matrix(gene) ~ as.matrix(cov))
  }
  cov <- obj.null$Z1
  gene <- obj.null$y

  perturbation <- as.matrix(perturbation)
  if (ncol(perturbation) == 1) {
    m <- 1
    perturbation <- t(perturbation)
  } else {
    m <- nrow(perturbation)
  }
  p.value <- rep(NA, m)
  p.value.NA <- rep(NA, m)
  Is.converge <- rep(NA, m)
  for (i in 1:m)
  {
    try({
      ina <- which(is.na(perturbation[i, ]))
      if (length(ina) > 0) {
        perturbation[i, ina] <- mean(perturbation[i, ], na.rm = TRUE)
      }
      perturb <- min(sum(perturbation[i, ]))
      if (perturb >= minperturb) {
        re <- negbinom_score_test(as.vector(perturbation[i, , drop = FALSE]), obj.null, cutoff = cutoff, log.p = log.p)

        p.value[i] <- re$p.value
        p.value.NA[i] <- re$p.value.NA
        Is.converge[i] <- re$Is.converge
      }
    })
  }

  return(list(p.value = p.value, p.value.NA = p.value.NA, Is.converge = Is.converge))
}

# Fit null model

negbinom_Z1_adj <- function(Z1) {
  q1 <- ncol(Z1)
  if (q1 >= 2) {
    if (sum(abs(Z1[, 1] - Z1[, 2])) == 0) {
      Z1 <- Z1[, -2]
      q1 <- q1 - 1
    }
  }
  Z1.qr <- qr(Z1)
  if (Z1.qr$rank < q1) {
    Z1.svd <- svd(Z1)
    Z1 <- Z1.svd$u[, 1:Z1.qr$rank]
  }
  return(Z1)
}

#' @importFrom stats model.matrix
#' @importFrom MASS glm.nb
negbinom_fit_null <- function(formula, data = NULL) {
  Z1 <- model.matrix(formula, data = data)
  Z1 <- negbinom_Z1_adj(Z1)

  glmfit <- glm.nb(formula, data = data)
  mu <- glmfit$fitted.values
  theta <- glmfit$theta
  pp <- theta / (theta + mu)

  W <- mu / (1 + mu / theta)
  M <- 1 / mu
  res <- glmfit$y - mu
  ZW <- t(Z1 * W)
  ZWZ_inv <- solve(t(Z1) %*% (Z1 * W))
  ZZWZ_inv <- Z1 %*% ZWZ_inv

  re <- list(y = glmfit$y, mu = mu, theta = theta, pp = pp, res = res, Z1 = Z1, W = W, M = M, ZW = ZW, ZZWZ_inv = ZZWZ_inv)
  class(re) <- "negbinom_null"

  return(re)
}

# Compute CGF function

negbinom_K0 <- function(t, mu, theta, pp, xx) {
  n.t <- length(t)
  out <- rep(0, n.t)

  for (i in 1:n.t) {
    t1 <- t[i]
    temp <- log(pp / (1 - (1 - pp) * exp(xx * pp * t1)))
    out[i] <- theta * sum(temp)
  }

  return(out)
}

negbinom_K1 <- function(t, mu, theta, pp, xx, q) {
  n.t <- length(t)
  out <- rep(0, n.t)

  for (i in 1:n.t) {
    t1 <- t[i]
    temp1 <- exp(-xx * pp * t1) - (1 - pp)
    temp2 <- xx * pp * (1 - pp)
    out[i] <- theta * sum(temp2 / temp1) - q
  }

  return(out)
}

negbinom_K2 <- function(t, mu, theta, pp, xx) {
  n.t <- length(t)
  out <- rep(0, n.t)

  for (i in 1:n.t) {
    t1 <- t[i]
    temp1 <- (1 - (1 - pp) * exp(xx * pp * t1))^2
    temp2 <- (1 - pp) * pp^2 * xx^2 * exp(xx * pp * t1)
    out[i] <- theta * sum(temp2 / temp1, na.rm = TRUE)
  }

  return(out)
}

negbinom_get_root <- function(init, mu, theta, pp, xx, q, tol = .Machine$double.eps^0.25, maxiter = 1000) {
  xx.pos <- sum(xx[which(xx > 0)])
  xx.neg <- sum(xx[which(xx < 0)])

  if (q >= xx.pos || q <= xx.neg) {
    return(list(root = Inf, n.iter = 0, Is.converge = TRUE))
  } else {
    t <- init
    K1_eval <- negbinom_K1(t, mu, theta, pp, xx, q)
    prevJump <- Inf
    rep <- 1
    repeat {
      K2_eval <- negbinom_K2(t, mu, theta, pp, xx)
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

      newK1 <- negbinom_K1(tnew, mu, theta, pp, xx, q)
      if (sign(K1_eval) != sign(newK1)) {
        if (abs(tnew - t) > prevJump - tol) {
          tnew <- t + sign(newK1 - K1_eval) * prevJump / 2
          newK1 <- negbinom_K1(tnew, mu, theta, pp, xx, q)
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
negbinom_get_prob <- function(zeta, mu, theta, pp, xx, q, log.p = FALSE) {
  k1 <- negbinom_K0(zeta, mu, theta, pp, xx)
  k2 <- negbinom_K2(zeta, mu, theta, pp, xx)

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

negbinom_add_logp <- function(p1, p2) {
  p1 <- -abs(p1)
  p2 <- -abs(p2)
  maxp <- max(p1, p2)
  minp <- min(p1, p2)

  return(maxp + log(1 + exp(minp - maxp)))
}

#' @importFrom stats pchisq
negbinom_get_pvalue <- function(q, mu, theta, pp, xx, cutoff = 2, log.p = FALSE) {
  m1 <- sum(mu * pp * xx)
  var1 <- sum(mu / (1 + mu / theta) * xx^2)
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
    out.uni1 <- negbinom_get_root(0, mu = mu, theta = theta, pp = pp, xx = xx, q = q)
    out.uni2 <- negbinom_get_root(0, mu = mu, theta = theta, pp = pp, xx = xx, q = qinv)
    if (out.uni1$Is.converge == TRUE && out.uni2$Is.converge == TRUE) {
      p1 <- tryCatch(negbinom_get_prob(out.uni1$root, mu, theta, pp, xx, q, log.p = log.p), error = function(e) {
        if (log.p) {
          return(pval.noadj - log(2))
        } else {
          return(pval.noadj / 2)
        }
      })
      p2 <- tryCatch(negbinom_get_prob(out.uni2$root, mu, theta, pp, xx, qinv, log.p = log.p), error = function(e) {
        if (log.p) {
          return(pval.noadj - log(2))
        } else {
          return(pval.noadj / 2)
        }
      })
      if (log.p) {
        pval <- negbinom_add_logp(p1, p2)
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
    return(negbinom_get_pvalue(q, mu, theta, pp, xx, cutoff = cutoff * 2, log.p = log.p))
  } else {
    return(list(p.value = pval, p.value.NA = pval.noadj, Is.converge = Is.converge, score = score))
  }
}

negbinom_score_test <- function(xx, obj.null, cutoff = 2, log.p = FALSE) {
  if (class(obj.null) != "negbinom_null") {
    stop("obj.null should be a returned object from negbinom_fit_null")
  }

  y <- obj.null$y
  mu <- obj.null$mu
  theta <- obj.null$theta
  pp <- obj.null$pp
  res <- obj.null$res

  xx1 <- xx - obj.null$ZZWZ_inv %*% (obj.null$ZW %*% xx)
  q <- sum(xx1 * pp * y)
  out <- negbinom_get_pvalue(q = q, mu = mu, theta = theta, pp = pp, xx = xx1, cutoff = cutoff, log.p = log.p)

  return(out)
}
