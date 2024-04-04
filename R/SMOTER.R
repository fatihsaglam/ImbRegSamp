#' @title SMOTE for Regression
#'
#' @description \code{SMOTER} applies undersampling and oversampling to
#' imbalanced regression data sets.
#'
#' @param x feature matrix or dataframe. Only numeric variables for now.
#' @param y target variable.
#' @param phi Relevance values for each sample. Default is NULL. If NULL,
#' inverse of kernel density estimation is used.
#' @param thresh_rel threshold to determine rare values. Default is 0.5.
#' @param perc_ov_lower percentage of oversampling of lower rare samples. If
#' NULL, it is automatically determined by algorithm. Cannot be lower than 1.
#' @param perc_ov_upper percentage of oversampling of upper rare samples. If
#' NULL, it is automatically determined by algorithm. Cannot be lower than 1.
#' @param perc_un percentage of undersampling of non-rare samples. If NULL, it
#' is automatically determined by algorithm. Cannot be higher than 1.
#' @param k number of neighbors for links
#' @param method bandwith estimation method used in relevance value calculation.
#' Default is "ucv".
#'
#' @details
#' Despite its name, it can undersample and oversample imbalanced data sets. It is
#' used to better estimate the rare values in regression. Algorithm is a little
#' different to UBL package and more similar to original paper (Torgo et al., 2013).
#' Relevance values are the inverse of kernel density estimations using unbiased
#' cross-validation bandwidth.There are three classes: lower rare,
#' not rare and upper rare. Lower rare class is samples that satisfy \eqn{\phi > treshold}
#' and \eqn{y < \tilde{y}}. Upper rare class is samples that satisfy \eqn{\phi > treshold}
#' and \eqn{y > \tilde{y}}. Other samples are considered as not rare.
#' \eqn{\tilde{y}} is median of \eqn{y}.
#'
#' @return an list object which includes:
#'  \item{x_new}{SMOTEd feature matrix}
#'  \item{y_new}{SMOTEd target variable}
#'  \item{x_syn}{Synthetic feature matrix}
#'  \item{y_syn}{Synthetic target variable}
#'  \item{phi}{relevance function values for y}
#'  \item{y}{original y to be used to calculate new target values}
#'  \item{h}{bandwidth value used to calculate densities for relevance values}
#'
#' @references
#' Torgo, L., Ribeiro, R. P., Pfahringer, B., & Branco, P. (2013, September).
#' Smote for regression. In Portuguese conference on artificial intelligence
#' (pp. 378-389). Berlin, Heidelberg: Springer Berlin Heidelberg.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @examples
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' m_SMOTER <- SMOTER(x = x, y = y)
#'
#' plot(x, y)
#' plot(m_SMOTER$x_new, m_SMOTER$y_new)
#'
#' @importFrom ks kde
#' @importFrom stats bw.ucv
#' @importFrom stats median
#'
#' @rdname SMOTER
#' @export

SMOTER <-
  function(x,
           y,
           thresh_rel = 0.5,
           phi = NULL,
           perc_ov_lower = NULL,
           perc_ov_upper = NULL,
           perc_un = NULL,
           k = 5,
           method = "ucv") {
    data <- as.matrix(cbind(x, y))
    n <- nrow(data)
    p <- ncol(data) - 1

    h <- NULL
    if (is.null(phi)) {
      m_rel <- relevance_density(y = y, method = method)
      h <- m_rel$h
      phi <- m_rel$rel
    } else {
      if (length(phi) != n) {
        stop("phi must be equal length to y")
      }
      if (!is.vector(phi)) {
        stop("phi must be vector of length y")
      }
      if (any(phi < 0)) {
        stop("phi cannot be negative")
      }
    }

    i_rare_lower <- which(phi > thresh_rel & y < median(y))
    i_rare_upper <- which(phi > thresh_rel & y > median(y))
    i_notRare <- which(phi <= thresh_rel)

    data_notRare <- data[i_notRare, ]
    n_notRare <- nrow(data_notRare)
    data_rare_lower <- data[i_rare_lower,]
    n_rare_lower <- nrow(data_rare_lower)
    data_rare_upper <- data[i_rare_upper,]
    n_rare_upper <- nrow(data_rare_upper)

    k_lower <- min(k, n_rare_lower - 1)
    k_upper <- min(k, n_rare_upper - 1)

    n_effbump <- sum(n_notRare > 0,
                     n_rare_lower > 0,
                     n_rare_upper > 0)

    if (is.null(perc_ov_lower)) {
      perc_ov_lower <- n / n_effbump / n_rare_lower
      if (perc_ov_lower < 1) {
        perc_ov_lower <- 1
      }
    } else {
      if (perc_ov_lower < 1) {
        perc_ov_lower <- 1
      }
    }
    if (is.null(perc_ov_upper)) {
      perc_ov_upper <- n / n_effbump / n_rare_upper
      if (perc_ov_upper < 1) {
        perc_ov_upper <- 1
      }
    } else {
      if (perc_ov_upper < 1) {
        perc_ov_upper <- 1
      }
    }
    if (is.null(perc_un)) {
      perc_un <- n / n_effbump / n_notRare
      if (perc_un > 1) {
        perc_un <- 1
      }
    } else {
      if (perc_un > 1) {
        perc_un <- 1
      }
    }

    cat(
      "perc_ov_lower:",
      perc_ov_lower,
      "| perc_ov_upper:",
      perc_ov_upper,
      "| perc_un:",
      perc_un,
      "\n"
    )
    cat(
      "n_rare_lower:",
      n_rare_lower,
      "| n_rare_upper:",
      n_rare_upper,
      "| n_notRare:",
      n_notRare,
      "\n"
    )

    i_notRare_undersampled <-
      sample(1:n_notRare, round(n_notRare * perc_un))
    data_notRare_undersampled <-
      data_notRare[i_notRare_undersampled,]

    data_syn_lower <-
      generator_SMOTER(data_rare = data_rare_lower,
                       perc_ov = perc_ov_lower,
                       k = k_lower)
    data_syn_upper <-
      generator_SMOTER(data_rare = data_rare_upper,
                       perc_ov = perc_ov_upper,
                       k = k_upper)

    data_syn <- rbind(data_syn_lower,
                      data_syn_upper)
    data_new <- rbind(data_notRare_undersampled,
                      data_syn,
                      data_rare_lower,
                      data_rare_upper)

    results <- list(
      x_new = data_new[, 1:p],
      y_new = data_new[, p + 1],
      x_syn = data_syn[, 1:p],
      y_syn = data_syn[, p + 1],
      phi = phi,
      h = h
    )

    return(results)
  }
