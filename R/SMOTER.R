#' @title SMOTE for Regression
#'
#' @description \code{SMOTER} applies undersampling and oversampling to
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
#' @param k number of neighbors for links.
#' @param rel_method method for relevance function. Default is "PCHIP". Choices
#' are "PCHIP" and "density". Ignored if phi is given.
#' @param ... relevance function settings.
#'
#' @details
#' Despite its name, it can both undersample and oversample imbalanced data sets. It is
#' used to better estimate the rare values in regression. Algorithm is from
#' Torgo et al. (2013). There are three classes: lower rare, not rare and upper
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
#' Torgo, L., Ribeiro, R. P., Pfahringer, B., & Branco, P. (2013, September).
#' Smote for regression. In Portuguese conference on artificial intelligence
#' (pp. 378-389). Berlin, Heidelberg: Springer Berlin Heidelberg.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @examples
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' m_SMOTER <- SMOTER(x = x, y = y)
#'
#' plot(x, y)
#' plot(m_SMOTER$x_new, m_SMOTER$y_new)
#'
#' @importFrom ks kde
#' @importFrom stats bw.ucv
#' @importFrom stats median
#'
#' @rdname SMOTER
#' @export

SMOTER <-
  function(x,
           y,
           thresh_rel = 0.5,
           phi = NULL,
           perc_ov_lower = NULL,
           perc_ov_upper = NULL,
           perc_un = NULL,
           k = 5,
           rel_method = "PCHIP",
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
      m_rel$rel_model <- list(method = "manual")
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

    data_notRare <- data[i_notRare,]
    n_notRare <- nrow(data_notRare)
    data_lowerRare <- data[i_lowerRare, ]
    n_lowerRare <- nrow(data_lowerRare)
    data_upperRare <- data[i_upperRare, ]
    n_upperRare <- nrow(data_upperRare)

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

    ### undersampling ###
    i_notRare_undersampled <-
      sample(1:n_notRare, round(n_notRare * perc_un))
    data_notRare_undersampled <-
      data_notRare[i_notRare_undersampled, ]
    ### undersampling finished ###

    k_lower <- min(k, n_lowerRare - 1)
    k_upper <- min(k, n_upperRare - 1)

    data_lowerSyn <-
      generator_SMOTER(data_rare = data_lowerRare,
                       perc_ov = perc_ov_lower,
                       k = k_lower)
    data_upperSyn <-
      generator_SMOTER(data_rare = data_upperRare,
                       perc_ov = perc_ov_upper,
                       k = k_upper)

    data_syn <- rbind(data_lowerSyn,
                      data_upperSyn)
    data_new <- rbind(
      data_notRare_undersampled,
      data_lowerSyn,
      data_upperSyn,
      data_lowerRare,
      data_upperRare
    )
    groups_new <-
      c(
        rep("notRare_undersampled", nrow(data_notRare_undersampled)),
        rep("lowerSyn", nrow(data_lowerSyn)),
        rep("upperSyn", nrow(data_upperSyn)),
        rep("lowerRare", nrow(data_lowerRare)),
        rep("upperRare", nrow(data_upperRare))
      )
    groups_new <- as.factor(groups_new)

    data_original <- rbind(data_lowerRare,
                           data_notRare,
                           data_upperRare)

    groups_original <- c(rep("lowerRare", nrow(data_lowerRare)),
                         rep("notRare", nrow(data_notRare)),
                         rep("upperRare", nrow(data_upperRare)))
    groups_original <- as.factor(groups_original)

    results <- list(
      x_new = data_new[, 1:p],
      y_new = data_new[, p + 1],
      groups_new = groups_new,
      x_original = data_original[, 1:p],
      y_original = data_original[, p + 1],
      groups_original = groups_original,
      x_syn = data_syn[, 1:p],
      y_syn = data_syn[, p + 1],
      phi = phi,
      rel_model = m_rel$rel_model
    )

    return(results)
  }
