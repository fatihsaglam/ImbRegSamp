#' @title Boosted Weights for SMOTE with Boosting for Regression (SMOTERWB)
#'
#' @description Calculation of Boosted Weights for SMOTE with Boosting for regression (SMOTERWB).
#'
#' @param x feature matrix.
#' @param y a factor class variable. Can work with more than two classes.
#' @param n_iter number of trees.
#' @param regressor regression method for weak learner. "lm" or "rpart". Default is "lm".
#' @param lr learning rate. Default is 0.1.
#' @param loss loss function for boosting update. Default is "linear" which is
#' no transformation. Other options are "square" and "exponential".
#'
#' @details
#' Calculation of Boosted Weights for SMOTERWB with Boosting.
#'
#' @return a vector of case weights.
#'  \item{w}{last weight vector of boosting process}
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#' @import rpart
#' @importFrom stats predict
#'
#'
#' @noRd

boosted_weights_regression <-
  function(x,
           y,
           n_iter = 20,
           regressor = "lm",
           lr = 0.1,
           loss = "linear") {
    dat <- data.frame(x, y)

    if (regressor == "lm") {
      weakLearner <- function(x, y, dat, w) {
        m <- lm.wfit(x = x, y = y, w = w)
        pred <- m$fitted.values
        return(pred)
      }
    } else if (regressor == "rpart") {
      weakLearner <- function(x, y, dat, w) {
        m <- rpart::rpart(
          y ~ .,
          data = dat,
          weights = w,
          control = rpart::rpart.control(
            minsplit = 3,
            cp = 0.01,
            maxdepth = 30
          )
        )
        pred <- predict(m, data = dat)
        return(pred)
      }
    }

    n <- length(y)
    w <- rep(1 / n, n)

    for (i in 1:n_iter) {
      pred <- weakLearner(x = x,
                          y = y,
                          dat = dat,
                          w = w)
      err <- abs(pred - y)
      i_mask <- w > 0
      w_masked <- w[i_mask]
      err_masked <- err[i_mask]
      err_max <- max(err_masked)

      if (err_max != 0) {
        err_masked <- err_masked / err_max
      }
      if (loss == "square") {
        err_masked <- err_masked ^ 2
      } else if (loss == "exponential") {
        err_masked <- 1 - exp(-err_masked)
      }

      err_estimator <- sum(w_masked * err_masked)

      if (err_estimator <= 0) {
        break
      } else if (err_estimator >= 0.5) {
        if (i > 1) {
          next
        }
      }
      beta <- err_estimator / (1 - err_estimator)
      alpha <- lr * log(1 / beta)
      w[i_mask] <- w[i_mask] * beta ^ (lr * (1 - err_masked))

      w <- w / sum(w)
    }

    return(w)
  }
