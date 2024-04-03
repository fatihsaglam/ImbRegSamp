#' @title Relevance function based on density
#'
#' @description Relevance value generation using kernel density estimation
#'
#' @param y target values.
#' @param y_new new target values for calculating relevance values. Default is y.
#' @param h bandwith value for kernel density estimation. Default is NULL. If
#' NULL, it is estimated based on \code{method}.
#' @param method method to be used in estimating bandwidth. Default is unbiased
#' cross validation.
#'
#' @details
#' \code{relevance_density} is a relevance function. It is used to determine
#' rare and not rare values in a regression problem. It calculates kernel densities
#' and inverse of densities are considered relevance values.
#'
#' @return a vector of relevance values.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @importFrom ks kde
#' @import stats
#' @importFrom methods getFunction
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

relevance_density <- function(y, y_new = y, h = NULL, method = "ucv") {
  match.arg(method, c("bcv", "nrd", "nrd0", "SJ", "ucv"))
  if (is.null(h)) {
    f_bw <- get(paste0("bw.", method))
    h <- f_bw(y)
  }
  dens <- kde(y, h = h, eval.points = y_new)$estimate
  dens[dens <= 1e-10] <- 1e-10
  rel <- 1/dens
  return(list(rel = rel,
              h = h))
}
