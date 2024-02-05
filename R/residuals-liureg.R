#' Liu Regression Residuals
#'
#' @param object An object of class \code{liureg}.
#' @param ... Not used in this implementation.
#'
#' @return The returned object is a vector or matrix whose columns
#' are Liu residuals for each \code{lambda}.
#' @author Murat Genç
#' @method residuals liureg
#' @export
#'
#' @seealso [liureg()], [pressliu()], [residuals()]
#'
#' @examples
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' residuals(liu.mod)
residuals.liureg <- function(object, ...){
  y <- drop(object$y)
  y <- y + object$ym
  y - predict.liureg(object)
}
