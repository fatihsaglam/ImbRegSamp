#' @title Relevance function based on sigmoid
#'
#' @description Relevance value generation using sigmoid function
#'
#' @param y target values.
#' @param y_new new target values for calculating relevance values. Default is y.
#' @param coef coefficient for IQR outlier range. Default is \eqn{1.5}.
#' @param adj_points a numeric vector of length 2. They are control points which
#' represents the points where relevance equals \eqn{0.5}. First value is the lower,
#' second is the upper \eqn{0.5}. Default is NULL and automatically decided based on
#' IQR statistics.
#' @param delta precision value. Default is \eqn{1e-4}.
#' @param k  decay factor. Default is \eqn{0.1}.
#'
#' @details
#' \code{relevance_sigmoid} is a relevance function. It is used to determine
#' rare and not rare values in a regression problem. It calculates relevance
#' based on a modified sigmoid function (Torgo & Ribeiro, 2021).
#'
#' @return a vector of relevance values.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @import stats
#'
#' @references
#' Torgo, L., & Ribeiro, R. (2009). Precision and recall for regression. In
#' Discovery Science: 12th International Conference, DS 2009, Porto, Portugal,
#' October 3-5, 2009 12 (pp. 332-346). Springer Berlin Heidelberg.
#'
#' @examples
#' ### regression dataset
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' ### relevance values
#' rel <- relevance_sigmoid(y)
#'
#' @rdname relevance_sigmoid
#' @export

relevance_sigmoid <-
  function(y,
           y_new = y,
           coef = 1.5,
           adj_points = NULL,
           delta = 1e-4,
           k = 0.5) {

    if (is.null(adj_points)) {
      IQR <- IQR(y)
      adj_points <-
        quantile(y, c(0.25, 0.75)) + c(-IQR * coef, IQR * coef)
    }

    adj_L <- adj_points[1]
    adj_H <- adj_points[2]

    s_L <- -log(1/delta - 1)/abs(adj_L*k)
    s_H <- log(1/delta - 1)/abs(adj_H*k)

    rel_L <- 1/(1 + exp(-s_L*(y_new - adj_L)))
    rel_H <- 1/(1 + exp(-s_H*(y_new - adj_H)))

    rel <- rep(1, length(y_new))
    rel[y_new > median(y)] <- rel_H[y_new > median(y)]
    rel[y_new <= median(y)] <- rel_L[y_new <= median(y)]

    return(list(
      rel = rel,
      rel_model = list(method = "sigmoid",
                       pars = list(y = y,
                                   adj_points = adj_points,
                                   delta = delta,
                                   k = k))
    ))
  }
