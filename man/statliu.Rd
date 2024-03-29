% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statliu.R
\name{statliu}
\alias{statliu}
\title{Liu Regression Statistics}
\usage{
statliu(obj)
}
\arguments{
\item{obj}{An object of class \code{liureg}.}
}
\value{
The return object is the statistics related to the Liu regression.
}
\description{
\code{statliu} computes the statistics related to the Liu regression.
}
\details{
\tabular{ll}{
\code{EDF} (Liu, 1993; Hastie et al., 2009) \tab Effective degrees of freedom, \eqn{n-\mathrm{trace}\left(2\mathbf{H}_\lambda\right)-\mathbf{H}_\lambda\mathbf{H}_\lambda^T} for each \eqn{\lambda} where \eqn{n} is the number of the observations in the design matrix and \eqn{\mathbf{H}_\lambda} is the hat matrix of Liu regression at \eqn{\lambda}.\cr
\tab \cr
\code{sigma2} \tab Computed \eqn{\hat{\sigma}^2} from the Liu regression for each \eqn{\lambda}.\cr
\tab \cr
\code{VAR} \tab Variance from the Liu regression for each \eqn{\lambda}.\cr
\tab \cr
\code{BIAS2} \tab Squared-bias from the Liu regression for each \eqn{\lambda}.\cr
\tab \cr
\code{MSE} \tab Mean squared error (MSE) from the Liu regression for each \eqn{\lambda}.\cr
\tab \cr
\code{FVal} \tab F-statistics value from the Liu regression for each \eqn{\lambda}.\cr
\tab \cr
\code{GCV} \tab Generalized cross-validation (GCV) from the Liu regression for each \eqn{\lambda}. The GCV is computed by \eqn{\frac{\mathrm{SSR}_{\lambda}}{n-1-\mathrm{trace}\left(\mathbf{H}_{\lambda}\right)}} where \eqn{\mathrm{SSR}_{\lambda}} is the residual sum of squares and \eqn{\mathrm{trace}\left(\mathbf{H}_{\lambda}\right)} is the trace of the hat matrix at corresponding value of \eqn{\lambda} from Liu regression.\cr
\tab \cr
\code{R2} \tab R-squared from the Liu regression for each \eqn{\lambda}.\cr
\tab \cr
\code{AdjR2} \tab Adjusted R-squared from the Liu regression for each \eqn{\lambda}.\cr
\tab \cr
}
}
\examples{
Hitters <- na.omit(Hitters)
X <- model.matrix(Salary ~ ., Hitters)[, -1]
y <- Hitters$Salary
lam <- seq(0, 1, 0.01)
liu.mod <- liureg(X, y, lam)
stats <- statliu(liu.mod)
print(stats)
}
\references{
Liu, K. (1993). A new class of blased estimate in linear regression.
\emph{Communications in Statistics-Theory and Methods}, \strong{22}(2), 393-402.
\doi{10.1080/03610929308831027}.

Hastie, T., Tibshirani, R., Friedman, J. H., Friedman, J. H. (2009).
The elements of statistical learning: data mining, inference,
and prediction (Vol. 2, pp. 1-758). \emph{New York: Springer}.
}
\seealso{
\code{\link[=liureg]{liureg()}}, \code{\link[=summary]{summary()}}, \code{\link[=pressliu]{pressliu()}}, \code{\link[=residuals]{residuals()}}
}
\author{
Murat Genç
}
