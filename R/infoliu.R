#' Information Criteria for Liu Regression
#'
#' For each value of \code{lambda}, infoliu calculates
#' the values of the AIC and BIC model selection criteria.
#' Model selection criteria are based on the degrees of the freedom,
#' \eqn{\texttt{df}=\mathrm{trace}\left(\mathbf{H}_\lambda\right)}
#' of the Liu regression model where \eqn{\mathbf{H}} is the hat matrix of
#' Liu regression model.
#'
#' @param obj A \code{liureg} object
#'
#' @return \code{infoliu} returns the matrix of information criteria
#'  for each value of the regularization parameter \code{lambda}.
#' @author Murat Genç
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' *IEEE Transaction on Automatic Control*, **9**(6), 716-723.
#' \doi{10.1109/TAC.1974.1100705}.
#'
#' Liu, K. (1993). A new class of blased estimate in linear regression.
#' *Communications in Statistics-Theory and Methods*, **22**(2), 393-402.
#' \doi{10.1080/03610929308831027}.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model.
#' *Annals of Statistics*, **6**(2), 461-464.
#' \doi{10.1214/aos/1176344136}.
#' @export
#'
#' @seealso [predict()], [summary()]
#' @examples
#' data("Hitters")
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' infoliu(liu.mod)
infoliu <- function(obj){
  diagH <- diagHliu(obj)
  SSRES <- obj$SSRES
  traceH <- colSums(diagH)
  n <- nrow(diagH)
  info <- data.frame(AIC = drop(n*log(SSRES/n)+2*traceH),
                     BIC = drop(n*log(SSRES)+log(n)*traceH))
  rownames(info) <- obj$lnames
  info
}
