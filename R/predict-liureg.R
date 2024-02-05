#' Predict method for liureg objects
#'
#' @param object A \code{liureg} object.
#' @param newdata A data frame of new values for \code{X} at which predictions are to be made.
#'    Can be a \code{data.frame} or a \code{matrix}.
#' @param ... Not used in this implementation.
#'
#' @return Depending on whether the \code{lambda} is a scalar or a vector,
#' the \code{predict.liureg} function returns a vector or matrix of predictions, respectively.
#' @author Murat Genç
#' @method predict liureg
#' @export
#'
#' @seealso [liureg()], [predict()], [summary()], [pressliu()], [residuals()]
#'
#' @examples
#' data("Hitters")
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' # Predictions based on original X matrix.
#' predict(liu.mod)
#' # Predictions based on newdata. newdata can be a matrix or a data.frame.
#' predict(liu.mod, newdata=X[1:5, ])
predict.liureg <- function(object, newdata, ...) {

  beta <- coef.liureg(object)
  p_1 <- nrow(beta) - 1

  if (missing(newdata)){
    newdata <- t(t(object$X) * drop(object$Xscale) + drop(object$Xm))
    rownames(newdata) <- object$rnames
  }
  if(is.null(rownames(newdata))) rownames(newdata) <- 1:nrow(newdata)
  if(!inherits(newdata, "matrix")) newdata <- as.matrix(newdata)
  if (ncol(newdata) != p_1) stop(paste0("The number of variables in newdata must be ", p_1))

  newdatafit <- predict_liureg(object, newdata)
  rownames(newdatafit) <- rownames(newdata)
  colnames(newdatafit) <- object$lnames
  newdatafit
}








