#' Computation of Liu Tuning Parameter
#'
#' \code{lamest} computes the Liu tuning parameters provided in the literature.
#'    The tuning parameter estimates are based on
#' \itemize{
#' \item Liu (1993) <doi:10.1080/03610929308831027>,
#' \item Ozkale and Kaciranlar (2007) <doi:10.1080/03610920601126522>,
#' \item Liu (2011) <doi:10.1016/j.jspi.2010.05.030>.
#' }
#'
#' @param obj An object of class \code{liureg}.
#' @param ... Not used in this implemetation.
#'
#' @details The \code{lamest} function computes the following tuning
#' parameter estimates available in the literature.
#' \tabular{ll}{
#' \code{lam.mm} (Liu, 1993) \tab \eqn{\displaystyle{1-\hat{\sigma}^{2}\left(\frac{\sum\limits _{j=1}^{p}\frac{1}{\lambda_{j}\left(1+\lambda_{j}\right)}}{\sum\limits _{j=1}^{p}\frac{\hat{\alpha}_{j}^{2}}{\left(1+\lambda_{j}\right)^{2}}}\right)}}\cr
#' \tab \cr
#' \code{lam.CL}  (Liu, 1993) \tab \eqn{\displaystyle{1-\hat{\sigma}^{2}\left(\frac{\sum\limits _{j=1}^{p}\frac{1}{\left(1+\lambda_{j}\right)}}{\sum\limits _{j=1}^{p}\frac{\lambda_{j}\hat{\alpha}_{j}^{2}}{\left(1+\lambda_{j}\right)^{2}}}\right)}}\cr
#' \tab \cr
#' \code{lam.opt} (Liu, 1993) \tab \eqn{\displaystyle{\frac{\sum\limits _{j=1}^{p}\left(\frac{\alpha_{j}^{2}-\sigma^{2}}{\left(1+\lambda_{j}\right)^{2}}\right)}{\sum\limits _{j=1}^{p}\left(\frac{\sigma^{2}+\lambda_{j}\alpha_{j}^{2}}{\lambda_{j}\left(1+\lambda_{j}\right)^{2}}\right)}}}\cr
#' \tab \cr
#' \code{lam.OK} (Ozkale and Kaciranlar, 2007; Liu, 2011) \tab \eqn{\displaystyle{\frac{\sum\limits _{i=1}^{n}\frac{\tilde{e}_{i}}{1-g_{ii}}\left(\frac{\tilde{e}_{i}}{1-h_{1-ii}}-\frac{\hat{e}_{i}}{1-h_{ii}}\right)}{\sum\limits _{i=1}^{n}\left(\frac{\tilde{e}_{i}}{1-g_{ii}}-\frac{\hat{e}_{i}}{1-h_{ii}}\right)^{2}}}} with \eqn{\hat{e}_{i}=y_{i}-\mathbf{x}_{i}^{T}\left(\mathbf{X}^{T}\mathbf{X}-\mathbf{x}_{i}\mathbf{x}_{i}^{T}\right)^{-1}\left(\mathbf{X}^{T}\mathbf{y}-\mathbf{x}_{i}y_{i}\right)} and \eqn{\tilde{e}_{i}=y_{i}-\mathbf{x}_{i}^{T}\left(\mathbf{X}^{T}\mathbf{X}+\mathbf{I}-\mathbf{x}_{i}\mathbf{x}_{i}^{T}\right)^{-1}\left(\mathbf{X}^{T}\mathbf{y}-\mathbf{x}_{i}y_{i}\right)} where \eqn{g_{ii}} and \eqn{h_{ii}} are the \eqn{i}th diagonal elements of \eqn{\mathbf{G}=\mathbf{X}\left(\mathbf{X}^{T}\mathbf{X}+\mathbf{I}\right)^{-1}\mathbf{X}^{T}} and \eqn{\mathbf{H=}\mathbf{X}\left(\mathbf{X}^{T}\mathbf{X}\right)^{-1}\mathbf{X}^{T}}, respectively.\cr
#' \tab \cr
#' \code{lam.GCV} \tab This is the \eqn{\lambda} value corresponding to the minimum of the generalized cross-validition (GCV) values. The GCV is computed by \eqn{\frac{\mathrm{SSR}_{\lambda}}{n-1-\mathrm{trace}\left(\mathbf{H}_{\lambda}\right)}} where \eqn{\mathrm{SSR}_{\lambda}} is the residual sum of squares and \eqn{\mathrm{trace}\left(\mathbf{H}_{\lambda}\right)} is the trace of the hat matrix at corresponding value of \eqn{\lambda} from Liu regression.\cr
#' \tab \cr
#' }
#'
#' @return The return object is the Liu tuning parameter
#' estimates based on the literature.
#' @author Murat Genç and Ömer Özbilen
#' @references Liu, K. (1993). A new class of blased estimate in linear regression.
#' *Communications in Statistics-Theory and Methods*, **22**(2), 393-402.
#' \doi{10.1080/03610929308831027}.
#'
#' Liu, X. Q. (2011). Improved Liu estimator in a linear regression model.
#' *Journal of Statistical Planning and Inference*, **141**(1), 189-196.
#' \doi{10.1016/j.jspi.2010.05.030}.
#'
#' Ozkale, M. R. and Kaciranlar, S. (2007). A prediction-oriented
#' criterion for choosing the biasing parameter in Liu estimation.
#' *Communications in Statistics-Theory and Methods*, **36**(10), 1889-1903.
#' \doi{10.1080/03610920601126522}.
#' Imdadullah, M., Aslam, M., and Altaf, S., (2017).
#' liureg: A Comprehensive R Package for the Liu Estimation of Linear Regression Model with
#' Collinear Regressors.  *The R Journal*, **9**(2), 232-247.
#' @export
#'
#' @seealso [liureg()], [predict()], [summary()], [pressliu()], [residuals()]
#' @examples
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' lamest(liu.mod)
lamest <- function(obj, ...) {
  lmest <- liuoptlamcpp(obj)
  colnames(lmest) <- "lambda"
  rownames(lmest) <- c("lam.mm", "lam.CL", "lam.opt", "lam.OK", "lam.GCV")
  lmest
}
