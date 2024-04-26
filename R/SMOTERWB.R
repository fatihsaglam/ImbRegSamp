#' @title SMOTEWB for Regression
#'
#' @description \code{SMOTERWB} applies undersampling and oversampling to
#' imbalanced regression data sets.
#'
#' @param x feature matrix or dataframe. Only numeric variables for now.
#' @param y target variable.
#' @param phi Relevance values for each sample. Default is NULL. If NULL,
#' inverse of kernel density estimation is used.
#' @param thresh_rel threshold to determine rare values. Default is 0.5.
#' @param perc_ov_lower percentage of oversampling of lower rare samples. If
#' NULL, it is automatically determined by algorithm. Cannot be lower than 1.
#' @param perc_ov_upper percentage of oversampling of upper rare samples. If
#' NULL, it is automatically determined by algorithm. Cannot be lower than 1.
#' @param perc_un percentage of undersampling of non-rare samples. If NULL, it
#' is automatically determined by algorithm. Cannot be higher than 1.
#' @param k_max to increase maximum number of neighbors. Determined automatically.
#' @param n_weak_classifier number of weak classifiers for boosting.
#' @param rel_method method for relevance function. Default is "PCHIP". Choices
#' are "PCHIP" and "density". Ignored if phi is given.
#' @param ... relevance function settings.
#'
#' @details
#' Despite its name, it can both undersample and oversample imbalanced data sets. It is
#' used to better estimate the rare values in regression. Algorithm is from
#' Sağlam & Cengiz (2022). There are three classes: lower rare, not rare and upper
#' rare. Lower rare class is samples that satisfy \eqn{\phi > treshold} and
#' \eqn{y < \tilde{y}}. Upper rare class is samples that satisfy
#' \eqn{\phi > treshold} and \eqn{y > \tilde{y}}. Other samples are considered
#' as not rare. \eqn{\tilde{y}} is median of \eqn{y}.
#'
#' @return an list object which includes:
#'  \item{x_new}{Balanced feature matrix}
#'  \item{y_new}{Balanced target variable}
#'  \item{groups_new}{categories of new data}
#'  \item{x_original}{original feature matrix}
#'  \item{y_original}{original target variable}
#'  \item{groups_original}{categories of original data}
#'  \item{x_syn}{Synthetic feature matrix}
#'  \item{y_syn}{Synthetic target variable}
#'  \item{phi}{relevance function values for y}
#'  \item{rel_model}{Details about relevance function. Can be used to calculate
#'  new relevance for test data.}
#'
#' @references
#' Sağlam, F., & Cengiz, M. A. (2022). A novel SMOTE-based resampling technique
#' trough noise detection and the boosting procedure. Expert Systems with
#' Applications, 200, 117023.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @examples
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' m_SMOTERWB <- SMOTERWB(x = x, y = y)
#'
#' plot(x, y)
#' plot(m_SMOTERWB$x_new, m_SMOTERWB$y_new)
#'
#' @importFrom ks kde
#' @importFrom stats bw.ucv
#' @importFrom stats median
#'
#' @rdname SMOTERWB
#' @export

SMOTERWB <-
  function(x,
           y,
           thresh_rel = 0.5,
           phi = NULL,
           type = 2,
           perc_ov_lower = NULL,
           perc_ov_upper = NULL,
           perc_un = NULL,
           k_max = NULL,
           n_weak_learner = 50,
           rel_method = "PCHIP",
           regressor = "lm",
           lr = 1,
           loss = "linear",
           ...) {
    data <- as.matrix(cbind(x, y))
    n <- nrow(data)
    p <- ncol(data) - 1

    if (is.null(phi)) {
      f_rel <- get(paste0("relevance_", rel_method))
      m_rel <- f_rel(y = y, ...)
      phi <- m_rel$rel
    } else {
      m_rel <- list()
      m_rel$rel_model <- list(
        method = "manual"
      )
      if (length(phi) != n) {
        stop("phi must be equal length to y")
      }
      if (!is.vector(phi)) {
        stop("phi must be vector of length y")
      }
      if (any(phi < 0)) {
        stop("phi cannot be negative")
      }
    }

    i_lowerRare <- (phi > thresh_rel & y < median(y))
    i_upperRare <- (phi > thresh_rel & y > median(y))
    i_notRare <- (phi <= thresh_rel)


    data_lowerRare <- data[i_lowerRare,]
    n_lowerRare <- nrow(data_lowerRare)
    data_notRare <- data[i_notRare, ]
    n_notRare <- nrow(data_notRare)
    data_upperRare <- data[i_upperRare,]
    n_upperRare <- nrow(data_upperRare)

    phi_lowerRare <- phi[i_lowerRare]
    phi_notRare <- phi[i_notRare]
    phi_upperRare <- phi[i_upperRare]

    n_effbump <- sum(n_notRare > 0,
                     n_lowerRare > 0,
                     n_upperRare > 0)

    if (is.null(perc_ov_lower)) {
      perc_ov_lower <- n / n_effbump / n_lowerRare
      if (perc_ov_lower < 1) {
        perc_ov_lower <- 1
      }
    } else {
      if (perc_ov_lower < 1) {
        perc_ov_lower <- 1
      }
    }
    if (is.null(perc_ov_upper)) {
      perc_ov_upper <- n / n_effbump / n_upperRare
      if (perc_ov_upper < 1) {
        perc_ov_upper <- 1
      }
    } else {
      if (perc_ov_upper < 1) {
        perc_ov_upper <- 1
      }
    }
    if (is.null(perc_un)) {
      perc_un <- n / n_effbump / n_notRare
      if (perc_un > 1) {
        perc_un <- 1
      }
    } else {
      if (perc_un > 1) {
        perc_un <- 1
      }
    }

    cat(
      "perc_ov_lower:",
      perc_ov_lower,
      "| perc_ov_upper:",
      perc_ov_upper,
      "| perc_un:",
      perc_un,
      "\n"
    )
    cat(
      "n_lowerRare:",
      n_lowerRare,
      "| n_upperRare:",
      n_upperRare,
      "| n_notRare:",
      n_notRare,
      "\n"
    )


    data_original <- rbind(
      data_notRare,
      data_lowerRare,
      data_upperRare
    )
    class_lowerRare <- "lowerRare"
    class_notRare <- "notRare"
    class_upperRare <- "upperRare"
    z <- c(rep(class_notRare, nrow(data_notRare)),
           rep(class_lowerRare, nrow(data_lowerRare)),
           rep(class_upperRare, nrow(data_upperRare)))
    z <- as.factor(z)
    phi <- c(phi_lowerRare, phi_notRare, phi_upperRare)

    # class_rare <- "rare"
    # class_notRare <- "notRare"
    # z <-
    #   as.factor(c(rep(class_rare, n_lowerRare), rep(class_notRare, n_notRare), rep(class_rare, n_upperRare)))
    x <- data_original[,1:p, drop = FALSE]
    y <- data_original[,p + 1]

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

    w <- w/(phi + 1e-5)
    w <- w/sum(w)

    w_lowerRare <- w[z == class_lowerRare]
    w_notRare <- w[z == class_notRare]
    w_upperRare <- w[z == class_upperRare]

    T_lowerRare <- (1/(3*n_lowerRare))
    T_notRare <- (1/(3*n_notRare))
    T_upperRare <- (1/(3*n_upperRare))

    # T_lowerRare <- T_lowerRare/(n*sum(T_lowerRare, T_notRare, T_upperRare))
    # T_notRare <- T_notRare/(n*sum(T_lowerRare, T_notRare, T_upperRare))
    # T_upperRare <- T_upperRare/(n*sum(T_lowerRare, T_notRare, T_upperRare))

    nl_lowerRare <- ifelse(w_lowerRare > T_lowerRare, "noise", "notNoise")
    nl_notRare <- ifelse(w_notRare > T_notRare, "noise", "notNoise")
    nl_upperRare <- ifelse(w_upperRare > T_upperRare, "noise", "notNoise")

    n_lowerRare_noise <- sum(nl_lowerRare == "noise")
    n_notRare_noise <- sum(nl_notRare == "noise")
    n_upperRare_noise <- sum(nl_upperRare == "noise")

    n_lowerRare_notNoise <- sum(nl_lowerRare == "notNoise")
    n_notRare_notNoise <- sum(nl_notRare == "notNoise")
    n_upperRare_notNoise <- sum(nl_upperRare == "notNoise")

    data_lowerRare_noise <- data_lowerRare[nl_lowerRare == "noise", , drop = FALSE]
    data_notRare_noise <-
      data_notRare[nl_notRare == "noise", , drop = FALSE]
    data_upperRare_noise <- data_upperRare[nl_upperRare == "noise", , drop = FALSE]

    data_lowerRare_notNoise <-
      data_lowerRare[nl_lowerRare == "notNoise", , drop = FALSE]
    data_notRare_notNoise <-
      data_notRare[nl_notRare == "notNoise", , drop = FALSE]
    data_upperRare_notNoise <-
      data_upperRare[nl_upperRare == "notNoise", , drop = FALSE]

    ### undersampling ####
    n_notRare_toBeRemoved <- round(n_notRare - n_notRare*perc_un)

    n_notRare_toBeRemoved_noise <- min(n_notRare_toBeRemoved, n_notRare_noise)
    n_notRare_toBeRemoved_notNoise <- n_notRare_toBeRemoved - n_notRare_toBeRemoved_noise

    i_notRare_toBeRemoved_noise <- sample(1:n_notRare_noise, n_notRare_toBeRemoved_noise)
    i_notRare_toBeRemoved_notNoise <- sample(1:n_notRare_notNoise, n_notRare_toBeRemoved_notNoise)

    if (sum(i_notRare_toBeRemoved_noise) == 0) {
      data_notRare_undersampled_noise <-
        data_notRare[nl_notRare == "noise", , drop = FALSE]
    } else {
      data_notRare_undersampled_noise <-
        data_notRare[nl_notRare == "noise", , drop = FALSE][-i_notRare_toBeRemoved_noise, , drop = FALSE]
    }


    if (sum(i_notRare_toBeRemoved_notNoise) == 0) {
      data_notRare_undersampled_notNoise <-
        data_notRare[nl_notRare == "notNoise", , drop = FALSE]
    }    else {
      data_notRare_undersampled_notNoise <-
        data_notRare[nl_notRare == "notNoise", , drop = FALSE][-i_notRare_toBeRemoved_notNoise, , drop = FALSE]
    }

    data_notRare_undersampled <- rbind(
      data_notRare_undersampled_noise,
      data_notRare_undersampled_notNoise
    )

    # i_notRare_undersampled <-
    #   sample(1:n_notRare, round(n_notRare * perc_un))
    # data_notRare_undersampled <-
    #   data_notRare[i_notRare_undersampled,]

    ### undersampling finished ###

    if (is.null(k_max)) {
      k_max_lower <- ceiling(n_notRare/n_lowerRare)
      k_max_upper <- ceiling(n_notRare/n_upperRare)
    } else {
      if (length(k_max) == 1) {
        k_max_lower <- k_max
        k_max_upper <- k_max
      } else if (length(k_max) == 2) {
        k_max_lower <- k_max[1]
        k_max_upper <- k_max[2]
      }
    }

    k_max_lower <- min(k_max_lower, n_lowerRare - 1)
    k_max_upper <- min(k_max_upper, n_upperRare - 1)

    data_syn_lower <-
      generator_SMOTERWB(
        data_rare = data_lowerRare,
        data_notRare = data_notRare,
        perc_ov = perc_ov_lower,
        k_max = k_max_lower,
        n_weak_learner = n_weak_learner,
        nl_rare = nl_lowerRare,
        nl_notRare = nl_notRare,
        class_rare = class_lowerRare,
        class_notRare = class_notRare,
        phi = phi_lowerRare
      )


    # data_rare = data_upperRare
    # data_notRare = data_notRare
    # perc_ov = perc_ov_upper
    # k_max = k_max_upper
    # n_weak_learner = n_weak_learner
    # nl_rare = nl_upperRare
    # nl_notRare = nl_notRare
    # class_rare = class_upperRare
    # class_notRare = class_notRare
    # phi = phi_upperRare

    data_syn_upper <-
      generator_SMOTERWB(
        data_rare = data_upperRare,
        data_notRare = data_notRare,
        perc_ov = perc_ov_upper,
        k_max = k_max_upper,
        n_weak_learner = n_weak_learner,
        nl_rare = nl_upperRare,
        nl_notRare = nl_notRare,
        class_rare = class_upperRare,
        class_notRare = class_notRare,
        phi = phi_upperRare
      )

    data_syn <- rbind(data_syn_lower,
                      data_syn_upper)
    data_new <- rbind(data_notRare_undersampled,
                      data_syn_lower,
                      data_syn_upper,
                      data_lowerRare,
                      data_upperRare)
    groups_new <- c(rep("notRare_undersampled", nrow(data_notRare_undersampled)),
                    rep("lowerSyn", nrow(data_syn_lower)),
                    rep("upperSyn", nrow(data_syn_upper)),
                    rep("lowerRare", nrow(data_lowerRare)),
                    rep("upperRare", nrow(data_upperRare)))
    groups_new <- as.factor(groups_new)

    nl_original = as.factor(c(nl_lowerRare, nl_notRare, nl_upperRare))

    results <- list(
      x_new = data_new[, 1:p],
      y_new = data_new[, p + 1],
      groups_new = groups_new,
      x_original = data_original[, 1:p],
      y_original = data_original[, p + 1],
      groups_original = z,
      x_syn = data_syn[, 1:p],
      y_syn = data_syn[, p + 1],
      nl_original = nl_original,
      w = w,
      noise_ratio = sum(nl_original == "noise")/n,
      phi = phi,
      rel_model = m_rel$rel_model
    )

    return(results)
  }

