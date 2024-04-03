#' @title Weighted Sum of Squared Error
#'
#' @description Calculation of weighted sum of squared error for weighted
#' regression models
#'
#' @param truth true values.
#' @param pred predictions.
#' @param weight a numerical vector of length same as truth and pred. Default is
#' NULL. If NULL, it is equal to all samples are equally weighted.
#'
#' @details
#' Calculates weighted sum of squared error. It can be used as a performance or
#' model selection metric in imbalanced regression.
#'
#' @return a numerical value.
#'
#' @author Fatih Sağlam, saglamf89@gmail.com
#'
#' @examples
#'
#' x <- rnorm(100)
#' err <- rnorm(100)
#' y <- 2 + x^2 + err
#'
#' x_train = x[1:70]
#' y_train = y[1:70]
#'
#' x_test = x[71:100]
#' y_test = y[71:100]
#'
#' ### resampling with SMOTER to predict rare values better
#' m_SMOTER <- SMOTER(x = x_train, y = y_train)
#'
#' data_train <- data.frame(x = m_SMOTER$x_new, y = m_SMOTER$y_new)
#' m_lm <- lm(y~., data = data_train)
#'
#' pred_test <- predict(m_lm, newdata = data.frame(x = x_test))
#'
#' ### classical sum of squares
#' SSE_notweighted <- WSSE(truth = y_test, pred = pred_test)
#' SSE_notweighted
#'
#' ### weights for test data
#' m_rel <- relevance_density(y = y_train, y_new = y_test)
#'
#' ### weighted sum of squares
#' SSE_weighted <- WSSE(truth = y_test, pred = pred_test, weight = m_rel$rel)
#' SSE_weighted
#'
#' @rdname WSSE
#' @export

WSSE <- function(truth, pred, weight = NULL) {
  if (is.null(weight)) {
    weight <- rep(1, length(truth))
  }
  if (any(weight < 0)) {
    stop("weight cannot be negative")
  }
  if (sum(weight) != length(truth)) {
    weight <- weight/sum(weight)*length(truth)
  }

  sum((truth - pred)^2*weight)
}
