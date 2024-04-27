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
                               nl_rare,
                               nl_notRare,
                               class_rare,
                               class_notRare,
                               phi) {
  n_rare <- nrow(data_rare)
  p <- ncol(data_rare) - 1
  n_notRare <- nrow(data_notRare)
  n_all <- n_rare + n_notRare
  if (n_rare == 0) {
    return(matrix(NA, nrow = 0, ncol = p + 1))
  }

  x <-
    rbind(data_rare[, 1:p, drop = FALSE], data_notRare[, 1:p, drop = FALSE])
  y <- c(data_rare[, p + 1], data_notRare[, p + 1])
  z <- as.factor(c(rep(class_rare, nrow(data_rare)), rep(class_notRare, nrow(data_notRare))))
  nl <- as.factor(c(nl_rare, nl_notRare))

  n_notRare_noise <- sum(nl_notRare == "noise")
  n_rare_noise <- sum(nl_rare == "noise")
  n_notRare_notNoise <- sum(nl_notRare == "notNoise")
  n_rare_notNoise <- sum(nl_rare == "notNoise")

  data_notRare_noise <-
    data_notRare[nl_notRare == "noise", , drop = FALSE]
  data_rare_noise <- data_rare[nl_rare == "noise", , drop = FALSE]
  data_notRare_notNoise <-
    data_notRare[nl_notRare == "notNoise", , drop = FALSE]
  data_rare_notNoise <-
    data_rare[nl_rare == "notNoise", , drop = FALSE]

  data_notNoise <- rbind(data_rare_notNoise, data_notRare_notNoise)
  z_notNoise <- as.factor(c(rep(class_rare, nrow(data_rare_notNoise)),
                            rep(class_notRare, nrow(data_notRare_notNoise))))

  k_max <- min(k_max, nrow(data_notNoise) - 1)

  NN <-
    knnx.index(data = data_notNoise[, 1:p, drop = FALSE],
               query = data_rare[, 1:p, drop = FALSE],
               k = k_max + 1)
  NN_temp <- matrix(data = NA,
                    nrow = n_rare,
                    ncol = k_max)
  NN_temp[nl_rare == "noise",] <-
    NN[nl_rare == "noise", -(k_max + 1)]
  NN_temp[nl_rare == "notNoise",] <- NN[nl_rare == "notNoise", -1]
  NN <- NN_temp

  k <- c()
  fl <- c()

  for (i in 1:n_rare) {
    cls <- z_notNoise[NN[i,]]

    if (all(cls == class_rare)) {
      k[i] <- k_max
    } else {
      k[i] <- which(cls == class_notRare)[1] - 1
      if (is.na(k[i])) {
        k[i] <- sum(cls == class_rare)
      }
    }

    if (k[i] == 0 & nl_rare[i] == "noise") {
      fl[i] <- "bad"
    }
    if (k[i] == 0 & nl_rare[i] == "notNoise") {
      fl[i] <- "lonely"
    }
    if (k[i] > 0) {
      fl[i] <- "good"
    }
  }

  n_syn <- round((perc_ov - 1) * n_rare)

  data_syn <- matrix(nrow = 0, ncol = p + 1)
  prop_sampling <- phi
  prop_sampling[fl == "bad"] <- 0

  repeat {
    i <-
      sample(1:n_rare,
             1,
             prob = pmax(prop_sampling, 1e-6))

    if (fl[i] == "lonely") {
      data_syn_step <- data_rare[i,]
      data_syn <- rbind(data_syn, data_syn_step)
    }
    if (fl[i] == "good") {
      NN_i <- NN[i, 1:k[i]]
      NN_i_rare <- NN_i[z_notNoise[NN_i] == class_rare]
      i_k <- sample(1:length(NN_i_rare), 1, prob = phi[NN_i_rare] + 1e-5)
      lambda <- runif(1)
      kk <- data_notNoise[NN_i_rare[i_k], , drop = FALSE]
      data_rare_i_temp <- data_rare[i, , drop = FALSE]
      data_syn_step <-
        data_rare_i_temp + (kk - data_rare_i_temp) * lambda
      data_syn <- rbind(data_syn, data_syn_step)
    }

    if (n_syn == nrow(data_syn)) {
      break
    }
  }

  return(data_syn)
}
