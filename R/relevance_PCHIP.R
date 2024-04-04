#' @title Piecewise Cubic Hermite Interpolating Polynomials
#'
#' @description Relevance function using piecewise cubic Hermite interpolating
#' polynomials (PCHIP)
#'
#' @param y target values.
#' @param y_new new target values for calculating relevance values. Default is y.
#' @param coef coefficient for IQR outlier range. Default is 1.5.
#' @param S_points a numeric vector of length 3. They are control points in y.
#' Lower than first value are lower rare, higher than third value are upper rare
#' area. Middle value is median and have zero relevance. Default is NULL and
#' automatically decided based on IQR statistics.
#' @param S_phi relevance values for respected points given in S_points. Default
#' is NULL and automatically decided.
#' @param S_der derivative of relevance values for respected points given in
#' S_points. If extreme points are given in S_points, all three values should
#' be zero. Default is NULL, which means all values are zero.
#'
#' @details
#' \code{relevance_PCHIP} is the most popular relevance function. It is used to
#' determine rare and not rare values in a regression problem. It calculates
#' interpolation values to determine weights of samples based on given bumps
#' (S_points).
#'
#' @return a vector of relevance values.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @importFrom stats IQR
#'
#' @examples
#'
#' ### regression dataset
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' ### relevance values
#' rel <- relevance_PCHIP(y)
#'
#' @rdname relevance_PCHIP
#' @export

relevance_PCHIP <-
  function(y,
           y_new = y,
           coef = 1.5,
           S_points = NULL,
           S_phi = NULL,
           S_der = NULL) {
    # y <- sort(y)

    if (is.null(S_points)) {
      IQR <- IQR(y)
      S_points <-
        quantile(as.numeric(y), c(0.25, 0.5, 0.75)) + c(-IQR * coef, 0, IQR * coef)
    }
    if (is.null(S_phi)) {
      S_phi <-
        c(as.numeric(any(y < S_points[1])), 0, as.numeric(any(y > S_points[3])))
    }
    if (is.null(S_der)) {
      S_der <- c(0, 0, 0)
    }

    s <- length(S_points)
    h <- diff(S_points)
    gamma <- diff(S_phi) / h
    a <- S_phi[-s]
    b <- S_der


    ### check_slopes ###
    for (k in 1:(s - 1)) {
      if (gamma[k] == 0) {
        b[k] <- b[k + 1] <- 0
      } else {
        alpha <- b[k] / gamma[k]
        beta <- b[k + 1] / gamma[k]
        if (b[k] != 0 & alpha < 0) {
          b[k] <- -b[k]
          alpha <- b[k] / gamma[k]
        }
        if (b[k + 1] != 0 & beta < 0) {
          b[k + 1] <- -b[k + 1]
          beta <- b[k + 1] / gamma[k]
        }
        tau1 <- 2 * alpha + beta - 3
        tau2 <- alpha + 2 * beta - 3
        if (tau1 > 0 & tau2 > 0 & alpha * (tau1 + tau2) < tau1 * tau2) {
          tau <- 3 * gamma[k] / sqrt(alpha ^ 2 + beta ^ 2)
          b[k] <- alpha * tau
          b[k + 1] <- beta * tau
        }
      }
    }
    ### end checking, new derivatives ###

    c <- (3 * gamma - 2 * b[-s] + b[-1]) / h
    d <- (b[-s] - 2 * gamma + b[-1]) / h ^ 2

    ii <- rep(1, length(y_new))
    for (j in 2:(s - 1)) {
      ii[S_points[j] <= y_new] <- j
    }

    diffs <- y_new - S_points[ii]

    phi <- a[ii] + b[ii] * diffs + c[ii] * diffs ^ 2 + d[ii] * diffs ^ 3
    phi[y_new < S_points[1]] <- 1
    phi[y_new > S_points[3]] <- 1

    return(list(
      rel = phi,
      rel_model = list(
        method = "PCHI",
        pars = list(
          coef = coef,
          S_points = S_points,
          S_phi = S_phi,
          S_der = S_der
        )
      )
    ))
  }
