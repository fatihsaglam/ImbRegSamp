#' @title Data Generation Algorithm of GSMOTER
#'
#' @description Generation of data in GSMOTER algorithm.
#'
#' @param data_rare rare dataset.
#' @param perc_ov percentage of oversampling.
#' @param k number of nearest neighbors.
#' @param alpha_sel selection method. Can be "minority", "majority" or "combined".
#' Default is "combined".
#' @param alpha_trunc truncation factor. A numeric value in \eqn{[-1,1]}.
#' Default is 0.5.
#' @param alpha_def deformation factor. A numeric value in \eqn{[0,1]}.
#' Default is 0.5
#'
#' @details
#' To be used inside GSMOTER to generate data when oversampling.
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

generator_GSMOTER <- function(data_rare, data_notRare, perc_ov, k,
                              alpha_sel, alpha_trunc, alpha_def) {
  n_rare <- nrow(data_rare)
  p_rare <- ncol(data_rare) - 1
  if (n_rare == 0) {
    return(matrix(NA, nrow = 0, ncol = p_rare + 1))
  }

  n_syn <- round((perc_ov - 1) * n_rare)

  m_rare2notRare <- FNN::get.knnx(data = data_notRare,
                             query = data_rare,
                             k = 1 + 1)
  NN_rare2notRare <- m_rare2notRare$nn.index[, -1, drop = FALSE]
  D_rare2notRare <- m_rare2notRare$nn.dist[, -1, drop = FALSE]

  m_rare2rare <- FNN::get.knnx(data = data_rare,
                             query = data_rare,
                             k = k + 1)
  NN_rare2rare <- m_rare2rare$nn.index[, -1, drop = FALSE]
  D_rare2rare <- m_rare2rare$nn.dist[, -1, drop = FALSE]

  C <- rep(floor(n_syn / n_rare), n_rare)
  n_needed <- n_syn - sum(C)
  ii <- sample(1:n_rare, n_needed)
  C[ii] <- C[ii] + 1

  data_syn <- matrix(NA, nrow = 0, ncol = p_rare + 1)

  for (i in 1:n_rare) {

    data_center <- data_rare[i, , drop = FALSE]

    for (j in 1:C[i]) {
      ### Surface ###
      if (alpha_sel == "minority") {
        i_selected_neighbor_rare <- sample(1:k, 1)
        data_surface <-
          data_rare[NN_rare2rare[i, i_selected_neighbor_rare], , drop = FALSE]
      } else if (alpha_sel == "majority") {
        i_selected_neighbor_notRare <- sample(1:1, 1)
        data_surface <-
          data_notRare[NN_rare2notRare[i, i_selected_neighbor_notRare], , drop = FALSE]
      } else {
        i_selected_neighbor_rare <- sample(1:k, 1)
        i_selected_neighbor_notRare <- sample(1:1, 1)

        if (D_rare2rare[i, i_selected_neighbor_rare] < D_rare2notRare[i, i_selected_neighbor_notRare]) {
          data_surface <-
            data_rare[NN_rare2rare[i, i_selected_neighbor_rare], , drop = FALSE]
        } else {
          data_surface <-
            data_notRare[NN_rare2notRare[i, i_selected_neighbor_notRare], , drop = FALSE]
        }
      }
      ###############

      ### Hyperball ###
      v_normal <- rnorm(p_rare + 1)
      e_sphere <- v_normal / norm(v_normal, type = "2")
      r <- runif(1)
      data_gen <- matrix(r ^ (1 / (p_rare + 1)) * e_sphere, nrow = 1)
      #################

      ### Vectors ###
      R <- norm(data_surface - data_center, type = "2")

      if (R == 0) {
        data_syn <- rbind(data_syn, data_center)
        next
      }

      e_parallel <- (data_surface - data_center) / R
      data_parallel_proj <- c(tcrossprod(data_gen, e_parallel))
      data_parallel <- data_parallel_proj * e_parallel
      data_perpendicular <- data_gen - data_parallel
      ###############

      ### Truncate ###
      if (abs(alpha_trunc - data_parallel_proj) > 1) {
        data_gen <- data_gen - 2 * data_parallel
      }
      ################

      ### Deform ###
      data_gen <- data_gen - alpha_def * data_perpendicular
      ##############

      ### Translate ###
      data_syn_step <- data_center + R * data_gen
      #################

      data_syn <- rbind(data_syn, data_syn_step)
    }

  }

  return(data_syn)
}
