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
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 geom_hline
#' @importFrom stats ppoints
#' @importFrom stats qbeta
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr n
#' @importFrom dplyr row_number
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#'
#' @param ps1 A vector in which each element is an observed p value.
#' @param ps2 A vector in which each element is an observed p value.
#' @param ps3 A vector in which each element is an observed p value.
#' @param ps4 A vector in which each element is an observed p value.
#' @param ps5 A vector in which each element is an observed p value.
#' @param ps6 A vector in which each element is an observed p value.
#' @param methods A column of strings representing the methods.
#' @param ci A parameter to control tolerance.
#' @param p_hline A parameter to display a reasonable range of p values.
#' @param plot_title A string representing the title.
#' @param plot_colors A column of strings representing the colors.
#'
#' @return A QQ plot of p values.
#'
#' @export
#'
#' @examples
#' ps1 <- stats::runif(100, min = 0, max = 1)
#' ps2 <- stats::runif(200, min = 0, max = 1)
#' ps3 <- stats::runif(300, min = 0, max = 1)
#' ps4 <- stats::runif(100, min = 0, max = 1)
#' ps5 <- stats::runif(200, min = 0, max = 1)
#' ps6 <- stats::runif(300, min = 0, max = 1)
#' gg_qqplot(ps1)
#' gg_qqplot(ps1, ps2, ps3, ps4, ps5, ps6, methods = c("a", "b", "c", "d", "e", "f"), p_hline = 0.01, plot_title = "FUN")
gg_qqplot <- function(ps1, ps2 = NULL, ps3 = NULL, ps4 = NULL, ps5 = NULL, ps6 = NULL, methods = NULL, ci = 0.95, p_hline = NULL, plot_title = NULL, plot_colors = c(aa = "firebrick3", bb = "hotpink2", cc = "darkslategray4", dd = "royalblue4", ee = "darkslategray4", ff = "magenta4")) {
  if (is.null(ps2) && is.null(ps3) && is.null(ps4) && is.null(ps5) && is.null(ps6)) {
    if (is.null(methods)) {
      methods <- c("1")
    }
    if (is.null(plot_title)) {
      plot_title <- ""
    }
    gg_qqplot_1(ps1 = ps1, methods = methods, ci = ci, p_hline = p_hline, plot_title = plot_title, plot_colors = plot_colors)
  } else if (is.null(ps3) && is.null(ps4) && is.null(ps5) && is.null(ps6)) {
    if (is.null(methods)) {
      methods <- c("1", "2")
    }
    if (is.null(plot_title)) {
      plot_title <- ""
    }
    gg_qqplot_2(ps1 = ps1, ps2 = ps2, methods = methods, ci = ci, p_hline = p_hline, plot_title = plot_title, plot_colors = plot_colors)
  } else if (is.null(ps4) && is.null(ps5) && is.null(ps6)) {
    if (is.null(methods)) {
      methods <- c("1", "2", "3")
    }
    if (is.null(plot_title)) {
      plot_title <- ""
    }
    gg_qqplot_3(ps1 = ps1, ps2 = ps2, ps3 = ps3, methods = methods, ci = ci, p_hline = p_hline, plot_title = plot_title, plot_colors = plot_colors)
  } else if (is.null(ps5) && is.null(ps6)) {
    if (is.null(methods)) {
      methods <- c("1", "2", "3", "4")
    }
    if (is.null(plot_title)) {
      plot_title <- ""
    }
    gg_qqplot_4(ps1 = ps1, ps2 = ps2, ps3 = ps3, ps4 = ps4, methods = methods, ci = ci, p_hline = p_hline, plot_title = plot_title, plot_colors = plot_colors)
  } else if (is.null(ps6)) {
    if (is.null(methods)) {
      methods <- c("1", "2", "3", "4", "5")
    }
    if (is.null(plot_title)) {
      plot_title <- ""
    }
    gg_qqplot_5(ps1 = ps1, ps2 = ps2, ps3 = ps3, ps4 = ps4, ps5 = ps5, methods = methods, ci = ci, p_hline = p_hline, plot_title = plot_title, plot_colors = plot_colors)
  } else {
    if (is.null(methods)) {
      methods <- c("1", "2", "3", "4", "5", "6")
    }
    if (is.null(plot_title)) {
      plot_title <- ""
    }
    gg_qqplot_6(ps1 = ps1, ps2 = ps2, ps3 = ps3, ps4 = ps4, ps5 = ps5, ps6 = ps6, methods = methods, ci = ci, p_hline = p_hline, plot_title = plot_title, plot_colors = plot_colors)
  }
}

gg_qqplot_1 <- function(ps1, methods, ci, p_hline, plot_title, plot_colors) {
  dat1 <- data.frame(pvalue = ps1)

  df_NTC <- rbind(
    dat1 %>% mutate(method = methods[1])
  ) %>%
    mutate(method = factor(x = as.character(method), levels = methods, labels = methods)) %>%
    arrange(method)

  df <- df_NTC %>%
    group_by(method) %>%
    mutate(
      r = row_number(pvalue), expected = ppoints(n())[r],
      clower = qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = n() + 1 - r),
      cupper = qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = n() + 1 - r)
    ) %>%
    ungroup() %>%
    mutate()

  p_thresh <- 1e-8
  p_val <- df %>%
    mutate(clower = ifelse(method == df[which.max(df$r), ]$method, clower, NA), cupper = ifelse(method == df[which.max(df$r), ]$method, cupper, NA)) %>%
    mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>%
    ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) +
    geom_point(aes(color = method), size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = setNames(plot_colors[c("aa")], NULL), name = "Method") +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      legend.position = c(0.25, 0.8),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    (if (!is.null(p_hline)) geom_hline(yintercept = p_hline, linetype = "dashed") else NULL)

  p_val
}

gg_qqplot_2 <- function(ps1, ps2, methods, ci, p_hline, plot_title, plot_colors) {
  dat1 <- data.frame(pvalue = ps1)
  dat2 <- data.frame(pvalue = ps2)

  df_NTC <- rbind(
    dat1 %>% mutate(method = methods[1]),
    dat2 %>% mutate(method = methods[2])
  ) %>%
    mutate(method = factor(x = as.character(method), levels = methods, labels = methods)) %>%
    arrange(method)

  df <- df_NTC %>%
    group_by(method) %>%
    mutate(
      r = row_number(pvalue), expected = ppoints(n())[r],
      clower = qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = n() + 1 - r),
      cupper = qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = n() + 1 - r)
    ) %>%
    ungroup() %>%
    mutate()

  p_thresh <- 1e-8
  p_val <- df %>%
    mutate(clower = ifelse(method == df[which.max(df$r), ]$method, clower, NA), cupper = ifelse(method == df[which.max(df$r), ]$method, cupper, NA)) %>%
    mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>%
    ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) +
    geom_point(aes(color = method), size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = setNames(plot_colors[c("aa", "bb")], NULL), name = "Method") +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      legend.position = c(0.25, 0.8),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    (if (!is.null(p_hline)) geom_hline(yintercept = p_hline, linetype = "dashed") else NULL)

  p_val
}

gg_qqplot_3 <- function(ps1, ps2, ps3, methods, ci, p_hline, plot_title, plot_colors) {
  dat1 <- data.frame(pvalue = ps1)
  dat2 <- data.frame(pvalue = ps2)
  dat3 <- data.frame(pvalue = ps3)

  df_NTC <- rbind(
    dat1 %>% mutate(method = methods[1]),
    dat2 %>% mutate(method = methods[2]),
    dat3 %>% mutate(method = methods[3])
  ) %>%
    mutate(method = factor(x = as.character(method), levels = methods, labels = methods)) %>%
    arrange(method)

  df <- df_NTC %>%
    group_by(method) %>%
    mutate(
      r = row_number(pvalue), expected = ppoints(n())[r],
      clower = qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = n() + 1 - r),
      cupper = qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = n() + 1 - r)
    ) %>%
    ungroup() %>%
    mutate()

  p_thresh <- 1e-8
  p_val <- df %>%
    mutate(clower = ifelse(method == df[which.max(df$r), ]$method, clower, NA), cupper = ifelse(method == df[which.max(df$r), ]$method, cupper, NA)) %>%
    mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>%
    ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) +
    geom_point(aes(color = method), size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = setNames(plot_colors[c("aa", "bb", "cc")], NULL), name = "Method") +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      legend.position = c(0.25, 0.8),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    (if (!is.null(p_hline)) geom_hline(yintercept = p_hline, linetype = "dashed") else NULL)

  p_val
}

gg_qqplot_4 <- function(ps1, ps2, ps3, ps4, methods, ci, p_hline, plot_title, plot_colors) {
  dat1 <- data.frame(pvalue = ps1)
  dat2 <- data.frame(pvalue = ps2)
  dat3 <- data.frame(pvalue = ps3)
  dat4 <- data.frame(pvalue = ps4)

  df_NTC <- rbind(
    dat1 %>% mutate(method = methods[1]),
    dat2 %>% mutate(method = methods[2]),
    dat3 %>% mutate(method = methods[3]),
    dat4 %>% mutate(method = methods[4])
  ) %>%
    mutate(method = factor(x = as.character(method), levels = methods, labels = methods)) %>%
    arrange(method)

  df <- df_NTC %>%
    group_by(method) %>%
    mutate(
      r = row_number(pvalue), expected = ppoints(n())[r],
      clower = qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = n() + 1 - r),
      cupper = qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = n() + 1 - r)
    ) %>%
    ungroup() %>%
    mutate()

  p_thresh <- 1e-8
  p_val <- df %>%
    mutate(clower = ifelse(method == df[which.max(df$r), ]$method, clower, NA), cupper = ifelse(method == df[which.max(df$r), ]$method, cupper, NA)) %>%
    mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>%
    ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) +
    geom_point(aes(color = method), size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = setNames(plot_colors[c("aa", "bb", "cc", "dd")], NULL), name = "Method") +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      legend.position = c(0.25, 0.8),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    (if (!is.null(p_hline)) geom_hline(yintercept = p_hline, linetype = "dashed") else NULL)

  p_val
}

gg_qqplot_5 <- function(ps1, ps2, ps3, ps4, ps5, methods, ci, p_hline, plot_title, plot_colors) {
  dat1 <- data.frame(pvalue = ps1)
  dat2 <- data.frame(pvalue = ps2)
  dat3 <- data.frame(pvalue = ps3)
  dat4 <- data.frame(pvalue = ps4)
  dat5 <- data.frame(pvalue = ps5)

  df_NTC <- rbind(
    dat1 %>% mutate(method = methods[1]),
    dat2 %>% mutate(method = methods[2]),
    dat3 %>% mutate(method = methods[3]),
    dat4 %>% mutate(method = methods[4]),
    dat5 %>% mutate(method = methods[5])
  ) %>%
    mutate(method = factor(x = as.character(method), levels = methods, labels = methods)) %>%
    arrange(method)

  df <- df_NTC %>%
    group_by(method) %>%
    mutate(
      r = row_number(pvalue), expected = ppoints(n())[r],
      clower = qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = n() + 1 - r),
      cupper = qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = n() + 1 - r)
    ) %>%
    ungroup() %>%
    mutate()

  p_thresh <- 1e-8
  p_val <- df %>%
    mutate(clower = ifelse(method == df[which.max(df$r), ]$method, clower, NA), cupper = ifelse(method == df[which.max(df$r), ]$method, cupper, NA)) %>%
    mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>%
    ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) +
    geom_point(aes(color = method), size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = setNames(plot_colors[c("aa", "bb", "cc", "dd", "ee")], NULL), name = "Method") +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      legend.position = c(0.25, 0.8),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    (if (!is.null(p_hline)) geom_hline(yintercept = p_hline, linetype = "dashed") else NULL)

  p_val
}

gg_qqplot_6 <- function(ps1, ps2, ps3, ps4, ps5, ps6, methods, ci, p_hline, plot_title, plot_colors) {
  dat1 <- data.frame(pvalue = ps1)
  dat2 <- data.frame(pvalue = ps2)
  dat3 <- data.frame(pvalue = ps3)
  dat4 <- data.frame(pvalue = ps4)
  dat5 <- data.frame(pvalue = ps5)
  dat6 <- data.frame(pvalue = ps6)

  df_NTC <- rbind(
    dat1 %>% mutate(method = methods[1]),
    dat2 %>% mutate(method = methods[2]),
    dat3 %>% mutate(method = methods[3]),
    dat4 %>% mutate(method = methods[4]),
    dat5 %>% mutate(method = methods[5]),
    dat6 %>% mutate(method = methods[6])
  ) %>%
    mutate(method = factor(x = as.character(method), levels = methods, labels = methods)) %>%
    arrange(method)

  df <- df_NTC %>%
    group_by(method) %>%
    mutate(
      r = row_number(pvalue), expected = ppoints(n())[r],
      clower = qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = n() + 1 - r),
      cupper = qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = n() + 1 - r)
    ) %>%
    ungroup() %>%
    mutate()

  p_thresh <- 1e-8
  p_val <- df %>%
    mutate(clower = ifelse(method == df[which.max(df$r), ]$method, clower, NA), cupper = ifelse(method == df[which.max(df$r), ]$method, cupper, NA)) %>%
    mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue)) %>%
    ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) +
    geom_point(aes(color = method), size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = setNames(plot_colors[c("aa", "bb", "cc", "dd", "ee", "ff")], NULL), name = "Method") +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      legend.position = c(0.25, 0.8),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    (if (!is.null(p_hline)) geom_hline(yintercept = p_hline, linetype = "dashed") else NULL)

  p_val
}

revlog_trans <- function(base = exp(1)) {
  trans <- function(x) {
    -log(x, base)
  }
  inv <- function(x) {
    base^(-x)
  }
  scales::trans_new(
    name = paste("revlog-", base, sep = ""),
    transform = trans,
    inverse = inv,
    breaks = scales::log_breaks(base = base),
    domain = c(1e-100, Inf)
  )
}

# #' @title Draw QQ plot
# #'
# #' @description Draw QQ plot for p values.
# #'
# #' @importFrom ggplot2 ggplot
# #' @importFrom ggplot2 geom_ribbon
# #' @importFrom ggplot2 geom_point
# #' @importFrom ggplot2 geom_abline
# #' @importFrom ggplot2 xlab
# #' @importFrom ggplot2 ylab
# #' @importFrom ggplot2 theme_bw
# #' @importFrom ggplot2 theme
# #' @importFrom ggplot2 aes
# #' @importFrom ggplot2 element_blank
# #' @importFrom stats ppoints
# #' @importFrom stats qbeta
# #'
# #' @param ps A vector in which each element is an observed p value.
# #' @param ci A parameter to control tolerance.
# #'
# #' @return A QQ plot of p values.
# #'
# #' @export
# #'
# #' @examples
# #' gg_qqplot(stats::runif(100, min = 0, max = 1))
# gg_qqplot <- function(ps, ci = 0.95) {
#   n <- length(ps)
#   df <- data.frame(
#     observed = -log10(sort(ps)),
#     expected = -log10(ppoints(n)),
#     clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
#     cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
#   )
#   log10Pe <- expression(paste("Expected -log"[10], plain(P)))
#   log10Po <- expression(paste("Observed -log"[10], plain(P)))
#   ggplot(df) +
#     geom_ribbon(
#       mapping = aes(x = expected, ymin = clower, ymax = cupper),
#       alpha = 0.1
#     ) +
#     geom_point(aes(expected, observed), shape = 1, size = 3) +
#     geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
#     xlab(log10Pe) +
#     ylab(log10Po) +
#     theme_bw() +
#     theme(panel.grid = element_blank())
# }
