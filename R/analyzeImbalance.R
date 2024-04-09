#' @title Analyze Imbalance
#'
#' @description Gets information about imbalance ratios based on relevance function.
#'
#' @param x feature matrix or dataframe. Only numeric variables for now.
#' @param y target variable.
#' @param phi Relevance values for each sample. Default is NULL. If NULL,
#' inverse of kernel density estimation is used.
#' @param thresh_rel threshold to determine rare values. Default is 0.5.
#' @param rel_method method for relevance function. Default is "PCHIP". Choices
#' are "PCHIP" and "density". Ignored if phi is given.
#' @param ... relevance function settings.
#'
#' @details
#' \code{analyzeImbalance} determines rare instances based on specified relevance
#' function. Then creates a group variable to indicate the characteristics of
#' respected instance. Function reorders datasets so output "groups" variable is
#' based on output "y_original". Imbalance ratios are calculated based on number of
#' samples in both upper and lower rare area. It is useful to see the
#' consequences of specified relevance function.
#'
#' @return an list object which includes:
#'  \item{x_original}{original feature matrix}
#'  \item{y_original}{original target variable}
#'  \item{groups_original}{categories of original data}
#'  \item{phi}{relevance function values for y}
#'  \item{rel_model}{Details about relevance function. Can be used to calculate
#'  new relevance for test data.}
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @examples
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' analyzeImbalance(x = x, y = y)
#'
#'
#' @importFrom stats median
#'
#' @rdname analyzeImbalance
#' @export

analyzeImbalance <- function(
    x,
    y,
    phi = NULL,
    thresh_rel = 0.5,
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
  n_rare <- n_rare_lower + n_rare_upper

  imbRate <- if (n_rare > 0) n_notRare/n_rare else NA
  imbRate_lower <- if (n_rare_lower > 0) n_notRare/n_rare_lower else NA
  imbRate_upper <- if (n_rare_lower > 0) n_notRare/n_rare_upper else NA

  min_y_rare_upper <- if (n_rare_upper > 0) min(data_rare_upper[,p+1]) else NA
  max_y_rare_lower <- if (n_rare_lower > 0) max(data_rare_lower[,p+1]) else NA

  cat(
    " Number of samples:\n",
    "n:",
    n,
    "| n_rare:",
    n_rare,
    "| n_rare_lower:",
    n_rare_lower,
    "| n_rare_upper:",
    n_rare_upper,
    "| n_notRare:",
    n_notRare,
    "\n-----\n",
    "Imbalance ratios:\n",
    "imbRate:",
    if (is.na(imbRate)) NA else formatC(imbRate, digits = 3, format = "f"),
    "| imbRate lower",
    if (is.na(imbRate_lower)) NA else formatC(imbRate_lower, digits = 3, format = "f"),
    "| imbRate upper",
    if (is.na(imbRate_upper)) NA else formatC(imbRate_upper, digits = 3, format = "f"),
    "\n------\n",
    "Bumps:\n",
    "min y_rare_upper:",
    min_y_rare_upper,
    "| max y_rare_lower:",
    if (n_rare_lower > 0) max_y_rare_lower else "no rare",
    "\n"
  )

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
    x_original = data_original[, 1:p],
    y_original = data_original[, p + 1],
    groups_original = groups_original,
    phi = phi,
    rel_model = m_rel$rel_model,
    imbRate = imbRate,
    imbRate_lower = imbRate_lower,
    imbRate_upper = imbRate_upper,
    min_y_rare_upper = min_y_rare_upper,
    max_y_rare_lower = max_y_rare_lower,
    n_notRare = n_notRare,
    n_rare = n_rare,
    n_rare_upper = n_rare_upper,
    n_rare_lower = n_rare_lower
  )

  return(results)

}
