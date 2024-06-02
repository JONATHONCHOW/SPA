#' @title Draw QQ plot
#'
#' @description Draw QQ plot for p values.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom stats ppoints
#' @importFrom stats qbeta
#'
#' @param ps A vector in which each element is an observed p value.
#' @param ci A parameter to control tolerance.
#'
#' @return A QQ plot of p values.
#'
#' @export
#'
#' @examples
#' gg_qqplot(stats::runif(100, min = 0, max = 1))
gg_qqplot <- function(ps, ci = 0.95) {
  n <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw() +
    theme(panel.grid = element_blank())
}
