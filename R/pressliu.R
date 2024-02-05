#' Predicted Residual Sum of Squares (PRESS)
#'
#' \code{pressliu} computes the predicted residual sum of squares (PRESS) based on a
#' Liu regression model.
#'
#' @param obj A \code{liureg} object.
#' @param digits Decimal places in the columns of data frame of PRESS values.
#'    Can be an integer of vector of integers.
#' @param ... Not used in this implementation.
#'
#' @details The PRESS statistic is based on the predicted leave-one-out residual sum of squares.
#' The statistic is computed as \eqn{{\sum\limits _{i=1}^{n}\left(\frac{\hat{e}_{\lambda i}}{1-h_{1-ii}}-\frac{e_{i}\left(h_{1-ii}-\tilde{\mathbf{H}}_{\lambda-ii}\right)}{\left(1-h_{1-ii}\right)\left(1-h_{ii}\right)}\right)^{2}}}
#' where \eqn{h_{ii}} is the \eqn{i}th diagonal element of the hat matrix corresponding
#' to the least squares estimator, \eqn{h_{1-ii}} is the \eqn{i}th diagonal
#' element of the hat matrix of the Liu estimator and \eqn{e_{\lambda i}}
#' is the residual at the specific value of \eqn{\lambda}.
#'
#' @return The returned object is a vector of PRESS values computed for each \code{lambda.}.
#' @author Murat Genç, Ömer Özbilen
#' @references Liu, K. (1993). A new class of blased estimate in linear regression.
#' *Communications in Statistics-Theory and Methods*, **22**(2), 393-402.
#' \doi{10.1080/03610929308831027}.
#'
#' Ozkale, M. R. and Kaciranlar, S. (2007). A prediction-oriented
#' criterion for choosing the biasing parameter in Liu estimation.
#' *Communications in Statistics-Theory and Methods*, **36**(10), 1889-1903.
#' \doi{10.1080/03610920601126522}.
#' @export
#' @seealso [liureg()], [pressliu()], [residuals()]
#'
#' @examples
#' data("Hitters")
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' pressliu(liu.mod)
pressliu <- function(obj, digits = 5L, ...){
  press <- round(pressliucpp(obj), digits)
  rownames(press) <- obj$lnames
  colnames(press) <- "PRESS"
  press
}
