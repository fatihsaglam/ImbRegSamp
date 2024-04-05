#' @title Relevance function based on density
#'
#' @description Relevance value generation using kernel density estimation
#'
#' @param y target values.
#' @param y_new new target values for calculating relevance values. Default is y.
#' @param h bandwith value for kernel density estimation. Default is NULL. If
#' NULL, it is estimated based on \code{method}.
#' @param bw_method method to be used in estimating bandwidth. Default is unbiased
#' cross validation.
#'
#' @details
#' \code{relevance_density} is a relevance function. It is used to determine
#' rare and not rare values in a regression problem. It calculates kernel densities
#' and inverse of densities are considered relevance values (Steininger et al., 2021).
#'
#' @return a vector of relevance values.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @importFrom ks kde
#' @import stats
#' @importFrom methods getFunction
#'
#' @references
#' Steininger, M., Kobs, K., Davidson, P., Krause, A., & Hotho, A. (2021).
#' Density-based weighting for imbalanced regression. Machine Learning, 110,
#' 2187-2211.
#'
#' @examples
#'
#' ### regression dataset
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' ### relevance values
#' rel <- relevance_density(y)
#'
#' @rdname relevance_density
#' @export

relevance_density <-
  function(y,
           y_new = y,
           h = NULL,
           alpha = 1,
           bw_method = "ucv",
           epsilon = 1e-10) {
    match.arg(bw_method, c("bcv", "nrd", "nrd0", "SJ", "ucv"))
    if (is.null(h)) {
      f_bw <- get(paste0("bw.", bw_method))
      h <- f_bw(y)
    }
    dens <- kde(y, h = h, eval.points = y_new)$estimate
    dens_scaled <- (dens - min(dens)) / (max(dens) - min(dens))
    rel <- max(1 - alpha*dens, epsilon)

    return(list(
      rel = rel,
      rel_model = list(method = "density",
                       pars = list(h = h))
    ))
  }
