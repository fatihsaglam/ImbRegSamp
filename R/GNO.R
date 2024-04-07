#' @title Gaussian Noise Oversampling
#'
#' @description \code{GNO} applies undersampling and oversampling to
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
#' @param rel_method method for relevance function. Default is "PCHIP". Choices
#' are "PCHIP" and "density". Ignored if phi is given.
#' @param pert effects the variance of noises. Default is 0.02, same as the paper
#' (Branco et al., 2019).
#' @param ... relevance function settings.
#'
#' @details
#' GNO (Branco et al., 2019) generates synthetic data by adding Gaussian noise.
#' It compares the length of link to the median distance to other rare samples.
#' If the first is higher, SMOGN generates by adding noise. If the second is
#' higher, it generates sample on the link.
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
#' Branco, Paula, Luis Torgo, and Rita P. Ribeiro. 2019. "Pre-Processing
#' Approaches for Imbalanced Distributions in Regression." Neurocomputing
#' 343: 76–99. https://doi.org/https://doi.org/10.1016/j.neucom.2018.11.100.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @examples
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' m_GNO <- GNO(x = x, y = y)
#'
#' plot(x, y)
#' plot(m_GNO$x_new, m_GNO$y_new)
#'
#' @importFrom stats median
#'
#' @rdname GNO
#' @export

GNO <- function(x,
                  y,
                  phi = NULL,
                  thresh_rel = 0.5,
                  perc_ov_lower = NULL,
                  perc_ov_upper = NULL,
                  perc_un = NULL,
                  rel_method = "PCHIP",
                  pert = 0.02,
                  ...) {

  data <- as.matrix(cbind(x, y))
  n <- nrow(data)
  p <- ncol(data) - 1

  if (is.null(phi)) {
    f_rel <- get(paste0("relevance_", rel_method))
    m_rel <- f_rel(y = y, ...)
    phi <- m_rel$rel
  } else {
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

  i_rare_lower <- (phi > thresh_rel & y < median(y))
  i_rare_upper <- (phi > thresh_rel & y > median(y))
  i_notRare <- (phi <= thresh_rel)

  data_notRare <- data[i_notRare, ]
  n_notRare <- nrow(data_notRare)
  data_rare_lower <- data[i_rare_lower,]
  n_rare_lower <- nrow(data_rare_lower)
  data_rare_upper <- data[i_rare_upper,]
  n_rare_upper <- nrow(data_rare_upper)

  n_effbump <- sum(n_notRare > 0,
                   n_rare_lower > 0,
                   n_rare_upper > 0)


  if (is.null(perc_ov_lower)) {
    perc_ov_lower <- n / n_effbump / n_rare_lower
    if (perc_ov_lower < 1) {
      perc_ov_lower <- 1
    }
  } else {
    if (perc_ov_lower < 1) {
      perc_ov_lower <- 1
    }
  }
  if (is.null(perc_ov_upper)) {
    perc_ov_upper <- n / n_effbump / n_rare_upper
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
    "n_rare_lower:",
    n_rare_lower,
    "| n_rare_upper:",
    n_rare_upper,
    "| n_notRare:",
    n_notRare,
    "\n"
  )

  ### undersampling ###
  i_notRare_undersampled <-
    sample(1:n_notRare, round(n_notRare * perc_un))
  data_notRare_undersampled <-
    data_notRare[i_notRare_undersampled,]
  ### undersampling finished ###

  data_syn_lower <-
    generator_GNO(data_rare = data_rare_lower,
                    perc_ov = perc_ov_lower, pert = pert)
  data_syn_upper <-
    generator_GNO(data_rare = data_rare_upper,
                    perc_ov = perc_ov_upper, pert = pert)

  data_syn <- rbind(data_syn_lower,
                    data_syn_upper)
  data_new <- rbind(data_notRare_undersampled,
                    data_syn_lower,
                    data_syn_upper,
                    data_rare_lower,
                    data_rare_upper)
  groups_new <- c(rep("notRare_undersampled", nrow(data_notRare_undersampled)),
                  rep("syn_lower", nrow(data_syn_lower)),
                  rep("syn_upper", nrow(data_syn_upper)),
                  rep("rare_lower", nrow(data_rare_lower)),
                  rep("rare_upper", nrow(data_rare_upper)))
  groups_new <- as.factor(groups_new)

  data_original <- rbind(
    data_notRare,
    data_rare_lower,
    data_rare_upper
  )

  groups_original <- c(rep("notRare", nrow(data_notRare)),
                       rep("rare_lower", nrow(data_rare_lower)),
                       rep("rare_upper", nrow(data_rare_upper)))
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
