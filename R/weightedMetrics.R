#' @title Weighted Performance Metrics
#'
#' @description Calculation of weighted errors for weighted
#' regression models
#'
#' @param truth true values.
#' @param pred predictions.
#' @param weight a numerical vector of length same as truth and pred. Default is
#' NULL. If NULL, it is equal to all samples are equally weighted.
#'
#' @details
#' Calculates weighted errors. They can be used as a performance or
#' model selection metric in imbalanced regression when weights are determined
#' same as train data. An example is given below. Implemented methods are
#'
#' WSSE: Weighted sum of squared errors.
#'
#' WMSE: Weighted mean squared error.
#'
#' WMAD: Weighted mean absolute deviation.
#'
#' WMAPE: weighted mean absolute percentage error
#'
#' More will come.
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
#' ### weights for test data
#' m_rel <- relevance_PCHIP(y = y_train, y_new = y_test)
#'
#' ### weighted sum of squares
#' SSE_weighted <- WSSE(truth = y_test, pred = pred_test, weight = m_rel$rel)
#' SSE_weighted

#' ### weighted mean squared error
#' MSE_weighted <- WMSE(truth = y_test, pred = pred_test, weight = m_rel$rel)
#' MSE_weighted
#'
#' ### weighted mean absolute deviation
#' MAD_weighted <- WMAD(truth = y_test, pred = pred_test, weight = m_rel$rel)
#' MAD_weighted
#'
#' ### weighted mean absolute percentage error
#' MAPE_weighted <- WMAPE(truth = y_test, pred = pred_test, weight = m_rel$rel)
#' MAPE_weighted
#'
#' @rdname metrics
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

#'
#' @rdname metrics
#' @export

WMSE <- function(truth, pred, weight = NULL) {
  if (is.null(weight)) {
    weight <- rep(1, length(truth))
  }
  if (any(weight < 0)) {
    stop("weight cannot be negative")
  }
  if (sum(weight) != length(truth)) {
    weight <- weight/sum(weight)*length(truth)
  }

  mean((truth - pred)^2*weight)
}

#'
#' @rdname metrics
#' @export

WMAD <- function(truth, pred, weight = NULL) {
  if (is.null(weight)) {
    weight <- rep(1, length(truth))
  }
  if (any(weight < 0)) {
    stop("weight cannot be negative")
  }
  if (sum(weight) != length(truth)) {
    weight <- weight/sum(weight)*length(truth)
  }

  mean(abs(truth - pred)*weight)
}

#'
#' @rdname metrics
#' @export

WMAPE <- function(truth, pred, weight = NULL) {
  if (is.null(weight)) {
    weight <- rep(1, length(truth))
  }
  if (any(weight < 0)) {
    stop("weight cannot be negative")
  }
  if (sum(weight) != length(truth)) {
    weight <- weight/sum(weight)*length(truth)
  }

  mean(abs(truth - pred)/truth*weight)
}
