#' @title Data Generation Algorithm of GNOR
#'
#' @description Generation of data in GNOR algorithm.
#'
#' @param data_rare rare dataset.
#' @param perc_ov percentage of oversampling.
#' @param pert coefficient for noise variances.
#'
#' @details
#' To be used inside GNOR to generate data when oversampling.
#'
#' @return oversampled data.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @noRd

generator_GNOR <- function(data_rare, perc_ov, pert = 0.02) {
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

  data_syn <- matrix(NA, nrow = 0, ncol = p_rare + 1)
  sds <- apply(data_rare, 2, sd) +
    1e-6 ### for zero or near zero variances

  for (i in 1:n_rare) {
    for (j in 1:C[i]) {
      center <- data_rare[i, , drop = FALSE]
      syn <- center + rnorm(n = p_rare + 1, mean = 0, sd = sds*pert)

      data_syn <- rbind(data_syn,
                        syn)
    }
  }
  return(data_syn)
}
