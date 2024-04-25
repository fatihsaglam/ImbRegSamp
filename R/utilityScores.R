#' @title utilityScores
#'
#' @description Calculation of utility scores.
#'
#' @param y_test ...
#' @param pred_test ...
#' @param y_train ...
#' @param p ...
#'
#' @details
#' Internal function.
#'
#' @return a vector of utility scores
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#' @references
#' Torgo, L., & Ribeiro, R. (2007). Utility-based regression. In Knowledge
#' Discovery in Databases: PKDD 2007: 11th European Conference on Principles
#' and Practice of Knowledge Discovery in Databases, Warsaw, Poland,
#' September 17-21, 2007. Proceedings 11 (pp. 597-604). Springer Berlin
#' Heidelberg.
#'
#' @noRd

utilityScores <- function(y_test, pred_test, y_train, p = 0.5, loss = "absolute") {
  m_rel <- ImbRegSamp::relevance_PCHIP(y = y_train)
  bumps <- m_rel$rel_model$pars$S_points

  m_rel_y_new <- ImbRegSamp::relevance_PCHIP(
    y = y_train,
    y_new = y_test,
    S_points = bumps,
    S_phi = m_rel$rel_model$pars$S_phi,
    S_der = m_rel$rel_model$pars$S_der
  )

  m_rel_pred <- ImbRegSamp::relevance_PCHIP(
    y = y_train,
    y_new = pred_test,
    S_points = bumps,
    S_phi = m_rel$rel_model$pars$S_phi,
    S_der = m_rel$rel_model$pars$S_der
  )
  phi_y_test <- m_rel_y_new$rel
  phi_pred_test <- m_rel_pred$rel

  if (loss == "absolute") {
    loss <- abs(y_test - pred_test)
  } else if (loss == "squared") {
    loss <- (y_test - pred_test)^2
  }

  L_B <- ifelse(pred_test < y_test, abs(y_test - bumps[2]), abs(y_test - Inf))
  L_B <- pmin(L_B, Inf)

  L_C <- ifelse(pred_test < y_test, abs(y_test - bumps[1]), abs(y_test - bumps[3]))
  L_C <- pmin(L_C, Inf)

  gamma_B <- ifelse(loss < L_B, loss / L_B, 1)
  gamma_C <- ifelse(loss < L_C, loss / L_C, 1)

  phi_joint_p <- (1 - p) * phi_pred_test + p * phi_y_test

  utility_phi_p <- phi_y_test*(1 - gamma_B) - phi_joint_p*gamma_C

  return(utility_phi_p)
}
