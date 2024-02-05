#' Print method for liureg objects
#'
#' Prints coefficients paths for Liu regression models over a grid of values for
#' the regularization (biasing) parameter lambda.
#'
#' @param x An object of class \code{liureg}.
#' @param digits Number of decimal places in the coefficients data.frame.
#' @param ... Not used in this implementation.
#'
#' @return The returned object is a data.frame showing the coefficients path.
#' @author Murat Genç
#' @method print liureg
#' @export
#' @seealso [liureg()], [summary()], [pressliu()], [residuals()]
#'
#' @examples
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' print(liu.mod)
print.liureg <- function(x, digits = max(3, getOption("digits") - 3),...) {

  cat("\nLiu Regression Coefficient Estimates:\n\n")
  coefs = data.frame(coef.liureg(x), check.names = FALSE)
  print(round(coefs,digits=digits))
  cat("\n")

}
