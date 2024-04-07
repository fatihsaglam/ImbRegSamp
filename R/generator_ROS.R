#' @title Data Generation Algorithm of Random Oversampling
#'
#' @description Generation of data in ROS algorithm.
#'
#' @param data_rare rare dataset.
#' @param perc_ov percentage of oversampling.
#'
#' @details
#' To be used inside ROUS to generate data when oversampling.
#'
#' @return oversampled data.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @noRd


generator_ROS <- function(data_rare, perc_ov) {
  n_rare <- nrow(data_rare)
  p_rare <- ncol(data_rare) - 1
  if (n_rare == 0) {
    return(matrix(NA, nrow = 0, ncol = p_rare + 1))
  }

  n_syn <- round((perc_ov - 1) * n_rare)

  C <- rep(floor(n_syn / n_rare), n_rare)
  n_needed <- n_syn - sum(C)
  ii <- sample(1:n_rare, n_needed)
  C[ii] <- C[ii] + 1

  i_selected <- unlist(sapply(1:n_rare, function(m) {
    rep(m, C[m])
  }))

  data_syn <- data_rare[i_selected,,drop = FALSE]
  return(data_syn)
}
