#' @title Data Generation Algorithm of SMOGN
#'
#' @description Generation of data in SMOGN algorithm.
#'
#' @param data_rare rare dataset.
#' @param perc_ov percentage of oversampling.
#' @param k number of nearest neighbors.
#' @param pert coefficient for noise variances.
#'
#' @details
#' To be used inside SMOGN to generate data when oversampling.
#'
#' @return a matrix of oversampled data.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @importFrom FNN get.knn
#' @importFrom Rfast dista
#' @importFrom Rfast Dist
#' @importFrom stats runif
#'
#' @noRd

generator_SMOGN <- function(data_rare, perc_ov, k, pert = 0.02) {
  n_rare <- nrow(data_rare)
  p_rare <- ncol(data_rare) - 1
  if (n_rare == 0) {
    return(matrix(NA, nrow = 0, ncol = p_rare + 1))
  }

  n_syn <- round((perc_ov - 1) * n_rare)

  m_NN_rare2rare <- get.knn(data = data_rare[,1:p_rare], k = k)
  i_NN_rare2rare <- m_NN_rare2rare$nn.index
  d_NN_rare2rare <- m_NN_rare2rare$nn.dist
  d_rare2rare <- Dist(data_rare)

  maxD <- sapply(1:n_rare, function(i) {
    median(d_rare2rare[i,-i])
  })

  C <- rep(floor(n_syn / n_rare), n_rare)
  n_needed <- n_syn - sum(C)
  ii <- sample(1:n_rare, n_needed)
  C[ii] <- C[ii] + 1

  data_syn <- matrix(NA, nrow = 0, ncol = p_rare + 1)
  sds <- apply(data_rare, 2, sd) +
    1e-6 ### for zero or near zero variances

  for (i in 1:n_rare) {
    for (j in 1:C[i]) {
      i_nn_selected <- sample(1:k, 1)
      DistM_x <- d_NN_rare2rare[i,i_nn_selected]
      center <- data_rare[i, , drop = FALSE]
      if (DistM_x < maxD[i]) {
        target <- data_rare[i_NN_rare2rare[i,i_nn_selected], , drop = FALSE]
        r <- runif(1)
        syn <- center + r * (target - center)
      } else {
        pert <- min(maxD[i], pert)
        syn <- center + rnorm(n = p_rare + 1, mean = 0, sd = sds*pert)
      }

      data_syn <- rbind(data_syn,
                        syn)
    }
  }
  return(data_syn)
}
