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

utilityScores <- function(y_test, pred_test, y_train, p = 0.5) {
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
  phi_y_new <- m_rel_y_new$rel
  phi_pred <- m_rel_pred$rel

  loss <- abs(y_test - pred_test)

  L <- ifelse(pred_test < y_test, abs(y_test - bumps[1]), abs(y_test - bumps[3]))
  gamma <- ifelse(loss < L, loss / L, 1)

  phi_joint_p <- (1 - p) * phi_pred + p * phi_y_new
  utility_phi_p <- phi_y_new*(1 - gamma) - phi_joint_p*gamma

  return(utility_phi_p)
}
