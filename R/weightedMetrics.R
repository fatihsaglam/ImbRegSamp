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
#' - WSSE: Weighted sum of squared errors.
#'
#' - WMSE: Weighted mean squared error.
#'
#' - WMAD: Weighted mean absolute deviation.
#'
#' - WMAPE: weighted mean absolute percentage error.
#'
#' - MU: Mean utility.
#'
#' - NMU: Normalized mean utility.
#'
#' - precision for regression
#'
#' - recall for regression
#'
#' - Fbeta score for regression. Default beta is 1.
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

  sum((truth - pred)^2*weight)/sum(weight)
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

  sum(abs(truth - pred)*weight)/sum(weight)
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

  sum(abs(truth - pred)/truth*weight)/sum(weight)
}

#'
#' @rdname metrics
#' @export
MU <- function(truth, pred, phi = NULL, y_train) {
  U <- utilityScores(y_test = truth, pred_test = pred, y_train = y_train)
  MU <- mean(U)
  return(MU)
}

#'
#' @rdname metrics
#' @export
NMU <- function(truth, pred, phi = NULL, y_train) {
  U <- utilityScores(y_test = truth, pred_test = pred, y_train = y_train)
  NMU <- mean(U)/2 + 1/2
  return(NMU)
}

#'
#' @rdname metrics
#' @export
precisionReg <- function(truth, pred, phi = NULL, thresh = 0.5, y_train, p = 0.5) {
  u <- utilityScores(y_test = truth, pred_test = pred, y_train = y_train, p = p)
  m_rel_truth <- ImbRegSamp::relevance_PCHIP(y = y_train, y_new = truth)
  m_rel_pred <- ImbRegSamp::relevance_PCHIP(y = y_train, y_new = pred)
  z_truth <- as.numeric(m_rel_truth$rel > thresh)
  z_pred <- as.numeric(m_rel_pred$rel > thresh)

  sum((1 + u)[z_truth == 1 & z_pred == 1])/
    (sum((1 + m_rel_truth$rel)[z_truth == 1 & z_pred == 1]) +
       sum((2 - p*(1 - m_rel_truth$rel))[z_truth == 0 & z_pred == 1]))
}

#'
#' @rdname metrics
#' @export
recallReg <- function(truth, pred, phi = NULL, thresh = 0.5, y_train, p = 0.5) {
  u <- utilityScores(y_test = truth, pred_test = pred, y_train = y_train, p = p)
  m_rel_truth <- ImbRegSamp::relevance_PCHIP(y = y_train, y_new = truth)
  m_rel_pred <- ImbRegSamp::relevance_PCHIP(y = y_train, y_new = pred)
  z_truth <- as.numeric(m_rel_truth$rel > thresh)
  z_pred <- as.numeric(m_rel_pred$rel > thresh)

  sum((1 + u)[z_truth == 1 & z_pred == 1])/
    sum((1 + m_rel_truth$rel)[z_truth == 1])
}

#'
#' @rdname metrics
#' @export
Fbeta <- function(truth, pred, phi = NULL, thresh = 0.5, y_train, p = 0.5, beta = 1) {
  prec <- f_reg_precision(truth = truth, pred = pred, phi = phi, thresh = thresh, y_train = y_train, p = p)
  rec <- f_reg_recall(truth = truth, pred = pred, phi = phi, thresh = thresh, y_train = y_train, p = p)

  if (any(c(prec, rec) == 0)) {
    return(0)
  } else {
    return((beta^2 + 1)*prec*rec/(beta^2*prec + rec))
  }
}
