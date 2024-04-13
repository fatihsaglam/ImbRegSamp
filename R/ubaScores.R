#' @title ubaScores
#'
#' @description Calculation of ubaScores.
#'
#' @param u ...
#' @param z ...
#'
#' @details
#' ubaScores
#'
#' @return a vector of ubaScores
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#' @references
#' Torgo, L., & Ribeiro, R. (2007). Utility-based regression. In Knowledge
#' Discovery in Databases: PKDD 2007: 11th European Conference on Principles
#' and Practice of Knowledge Discovery in Databases, Warsaw, Poland,
#' September 17-21, 2007. Proceedings 11 (pp. 597-604). Springer Berlin
#' Heidelberg.
#'
#' @importFrom utils tail
#'
#' @noRd

ubaScores <- function(u, z) {
  o <- order(u, decreasing = FALSE)
  s <- z[o]

  n <- length(s)

  repeat{
    i_dec <- c(which(diff(s) < 0), n)
    if (length(i_dec) == 1) {
      break
    }
    where_dec <- i_dec[1]:(i_dec[2] - 1)
    s[where_dec] <-
      sum(s[where_dec] / (tail(where_dec, 1) - where_dec[1] + 1))
  }
  return(s)
}
