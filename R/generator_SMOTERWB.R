#' @title Data Generation Algorithm of SMOTERWB
#'
#' @description Generation of data in SMOTERWB algorithm.
#'
#' @param data_rare rare dataset.
#' @param data_notRare non-rare dataset.
#' @param perc_ov percentage of oversampling.
#' @param k_max to increase maximum number of neighbors. Determined automatically.
#' @param n_weak_learner number of weak learners for boosting.
#'
#' @details
#' To be used inside SMOTERWB to generate data when oversampling.
#'
#' @return oversampled data.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @importFrom FNN knnx.index
#' @importFrom stats runif
#'
#' @noRd

generator_SMOTERWB <- function(data_rare,
                               data_notRare,
                               perc_ov,
                               k_max,
                               n_weak_learner,
                               type = 2,
                               regressor = "lm",
                               lr = 0.1,
                               loss = "linear") {
  n_rare <- nrow(data_rare)
  p_rare <- ncol(data_rare) - 1
  n_notRare <- nrow(data_notRare)
  n_all <- n_rare + n_notRare
  if (n_rare == 0) {
    return(matrix(NA, nrow = 0, ncol = p_rare + 1))
  }

  x <-
    rbind(data_rare[, 1:p_rare, drop = FALSE], data_notRare[, 1:p_rare, drop = FALSE])
  y <- c(data_rare[, p_rare + 1], data_notRare[, p_rare + 1])

  class_rare <- "rare"
  class_notRare <- "notRare"
  z <-
    as.factor(c(rep(class_rare, n_rare), rep(class_notRare, n_notRare)))


  if (type == 1) {
    w <-
      boosted_weights_regression(
        x = x,
        y = y,
        n_iter = n_weak_learner,
        regressor = regressor,
        lr = lr,
        loss = loss
      )
  } else if (type == 2) {
    w <-
      boosted_weights(x = x, y = z, n_iter = n_weak_learner)
  }

  w_rare <- w[z == class_rare]
  w_notRare <- w[z == class_notRare]

  wclass_rare <- n_all / n_rare * 0.5
  wclass_notRare <- n_all / n_notRare * 0.5

  T_rare <- (1 / n_all) * wclass_rare
  T_notRare <- (1 / n_all) * wclass_notRare

  scl <- T_rare * n_rare + T_notRare * n_notRare

  T_rare <- T_rare / scl
  T_notRare <- T_notRare / scl

  nl_notRare <- ifelse(w_notRare > T_notRare, "noise", "notnoise")
  nl_rare <- ifelse(w_rare > T_rare, "noise", "notnoise")

  n_notRare_noise <- sum(nl_notRare == "noise")
  n_rare_noise <- sum(nl_rare == "noise")
  n_notRare_notnoise <- sum(nl_notRare == "notnoise")
  n_rare_notnoise <- sum(nl_rare == "notnoise")

  data_notRare_noise <-
    data_notRare[nl_notRare == "noise", , drop = FALSE]
  data_rare_noise <- data_rare[nl_rare == "noise", , drop = FALSE]
  data_notRare_notnoise <-
    data_notRare[nl_notRare == "notnoise", , drop = FALSE]
  data_rare_notnoise <-
    data_rare[nl_rare == "notnoise", , drop = FALSE]

  data_notnoise <- rbind(data_rare_notnoise, data_notRare_notnoise)
  y_notnoise <- c(rep(class_rare, n_rare_notnoise),
                  rep(class_notRare, n_notRare_notnoise))

  k_max <- min(k_max, nrow(data_notnoise) - 1)
  NN <-
    knnx.index(data = data_notnoise[, 1:p_rare, drop = FALSE],
               query = data_rare[, 1:p_rare, drop = FALSE],
               k = k_max + 1)
  NN_temp <- matrix(data = NA,
                    nrow = n_rare,
                    ncol = k_max)
  NN_temp[nl_rare == "noise",] <-
    NN[nl_rare == "noise", -(k_max + 1)]
  NN_temp[nl_rare == "notnoise",] <- NN[nl_rare == "notnoise", -1]
  NN <- NN_temp

  k <- c()
  fl <- c()

  for (i in 1:n_rare) {
    cls <- y_notnoise[NN[i,]]

    if (all(cls == class_rare)) {
      k[i] <- k_max
    } else {
      k[i] <- which(cls == class_notRare)[1] - 1
    }

    if (k[i] == 0 & nl_rare[i] == "noise") {
      fl[i] <- "bad"
    }
    if (k[i] == 0 & nl_rare[i] == "notnoise") {
      fl[i] <- "lonely"
    }
    if (k[i] > 0) {
      fl[i] <- "good"
    }
  }

  n_syn <- round((perc_ov - 1) * n_rare)
  C <- numeric(n_rare)
  n_good_and_lonely <- sum((fl == "good") + (fl == "lonely"))
  for (i in 1:n_rare) {
    if (fl[i] == "good" | fl[i] == "lonely") {
      C[i] <- ceiling(n_syn / n_good_and_lonely)
    }
  }
  n_diff <- n_syn - sum(C)
  ii <-
    sample(which(fl == "good" | fl == "lonely"), size = abs(n_diff))
  C[ii] <- C[ii] + n_diff / abs(n_diff)

  data_syn <- matrix(nrow = 0, ncol = p_rare + 1)
  for (i in 1:n_rare) {
    if (fl[i] == "lonely") {
      i_step <- rep(i, C[i])
      data_syn_step <- data_rare[i_step,]
      data_syn <- rbind(data_syn, data_syn_step)
    }
    if (fl[i] == "good") {
      if (C[i] == 0) {
        next
      }
      NN_i <- NN[i, 1:k[i]]
      i_k <- sample(1:k[i], C[i], replace = TRUE)
      lambda <- runif(C[i])
      kk <- data_notnoise[NN_i, , drop = FALSE]
      kk <- kk[i_k,]
      data_rare_i_temp <- data_rare[rep(i, C[i]), , drop = FALSE]
      data_syn_step <-
        data_rare_i_temp + (kk - data_rare_i_temp) * lambda
      data_syn <- rbind(data_syn, data_syn_step)
    }
  }

  return(data_syn)
}
