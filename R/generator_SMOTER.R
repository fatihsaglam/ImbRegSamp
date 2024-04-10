#' @title Data Generation Algorithm of SMOTER
#'
#' @description Generation of data in SMOTER algorithm.
#'
#' @param data_rare rare dataset.
#' @param perc_ov percentage of oversampling.
#' @param k number of nearest neighbors.
#'
#' @details
#' To be used inside SMOTER to generate data when oversampling.
#'
#' @return oversampled data.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @importFrom FNN get.knn
#' @importFrom Rfast dista
#' @importFrom stats runif
#'
#' @noRd

generator_SMOTER <- function(data_rare, perc_ov, k) {
  n_rare <- nrow(data_rare)
  p_rare <- ncol(data_rare) - 1
  if (n_rare == 0) {
    return(matrix(NA, nrow = 0, ncol = p_rare + 1))
  }

  n_syn <- round((perc_ov - 1) * n_rare)

  NN_rare2rare <- get.knn(data_rare[,1:p_rare], k = k)$nn.index
  C <- rep(floor(n_syn / n_rare), n_rare)
  n_needed <- n_syn - sum(C)
  ii <- sample(1:n_rare, n_needed)
  C[ii] <- C[ii] + 1

  data_syn <- matrix(NA, nrow = 0, ncol = p_rare + 1)

  for (i in 1:n_rare) {
    center <- data_rare[rep(i, C[i]), , drop = FALSE]
    target <-
      data_rare[NN_rare2rare[i,sample(1:k, C[i], replace = TRUE)],, drop = FALSE]
    r <- runif(C[i])
    syn <- center + r * (target - center)

    d1 <- diag(dista(syn, center, trans = FALSE))
    d2 <- diag(dista(syn, target, trans = FALSE))

    y_syn_old <- syn[, p_rare + 1, drop = FALSE]

    y_center <- center[, p_rare + 1]
    y_target <- target[, p_rare + 1]

    y_syn_new <- (d1 * y_target + d2 * y_center) / (d1 + d2)

    y_syn_new[(d1 + d2) == 0] <- y_syn_old[(d1 + d2) == 0]

    syn[, p_rare + 1] <- y_syn_new

    data_syn <- rbind(data_syn,
                      syn)
  }

  return(data_syn)
}
