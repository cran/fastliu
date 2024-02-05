#' Extract coefficient estimates from a liureg object
#'
#' Prints coefficient estimates from a
#' fitted \code{liureg} object.
#'
#' @param object A \code{liureg} object.
#' @param ... Not used in this implementation.
#'
#' @method coef liureg
#'
#' @return The returned object is a data.frame containing the coefficients path.
#' @author Murat Genç
#' @export
#'
#' @seealso [liureg()], [predict()], [summary()], [pressliu()], [residuals()]
#' @examples
#' data("Hitters")
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' coef(liu.mod)
coef.liureg <- function(object, ...){
  betaorj <- coef_liureg(object)
  rownames(betaorj) <- c("Intercept", object$cnames)
  colnames(betaorj) <- object$lnames
  betaorj
}
