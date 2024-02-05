#' Summarizing Liu Regression Fit
#'
#' \code{summary} method for \code{liureg} objects.
#'
#' @param object An object of class \code{liureg}.
#' @param digits Number of decimal places in the data frame of Liu regression statistics.
#' @param ... Not used in this implemetation.
#'
#' @details \code{summary.liureg} produces an object with S3 class \code{summary.liureg}.
#' The function returns a list of summary statistics of the Liu regression fit for the grid
#' of regularization parameter \eqn{\lambda} values. Each element of the output list includes:
#' \tabular{ll}{
#' \code{coefficients} \tab A \eqn{p\times 5} matrix with columns coefficient estimates,
#' scaled coefficient estimates, scaled standard errors, scaled \eqn{t-}values with corresponding
#' \eqn{p-}value.\cr
#' \tab \cr
#' \code{Statistics} \tab Liu related statistics \eqn{R^2}, \eqn{\textrm{adjusted}-R^2},
#' \eqn{F-}statistics, AIC, BIC and MSE values.\cr
#' \tab \cr
#' }
#'
#' @return The returned object is a list whose elements are Liu
#' regression coefficient estimates and statistics related to Liu regression.
#' @author Murat Genç
#' @method summary liureg
#' @export
#'
#' @seealso [liureg()], [coef()], [predict()], [residuals()]
#' @examples
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' summary(liu.mod)
summary.liureg <- function(object, digits, ...) {
  digits <- if (missing(digits)) digits <- c(5,5,5,3,4) else rep(digits, length.out=5)
  lambda <- object$lambda
  y <- object$y
  ym <- object$ym
  Xm <- object$Xm
  n <- length(y)
  p <- length(Xm)
  l <- length(lambda)

  vy <- var(y)
  Xm2 <- Xm^2

  lstat <- statliu(object)
  linfo <- infoliu(object)


  liucoef <- coef(object)

  beta1 <- object$coefliu #Standardized beta1
  beta0 <- ym - colSums(beta1 * drop(Xm)) #Standardized beta0
  liucoefstd <- rbind(beta0,beta1)

  covliu <- covliu(object)
  SEmat <- matrix(NA,p,l)
  for(k in 1:l){
    SEmat[,k] <- sqrt(diag(covliu[[k]]))
  }
  rownames(SEmat) <- rownames(beta1)
  colnames(SEmat) <- colnames(beta1)

  tbeta1 <- beta1 / SEmat
  pbeta1 <- 2 * (1 - pnorm(abs(tbeta1)))

  sebeta0<-rep(0,l)
  tbeta0<-rep(0,l)
  pbeta0<-rep(0,l)

  coefmat <- matrix(list(), l, 1)
  rownames(coefmat) <- colnames(beta1)
  colnames(coefmat) <- "Coefficients"
  for (i in 1:l) {
    sebeta0[i]<-sqrt(vy/n + sum(Xm2 * SEmat[, i]^2))  #Xm2 %*% SEmat[, i]^2
    tbeta0[i]<-as.numeric(beta0[i]/sebeta0[i])
    pbeta0[i]<-2*(1-pnorm(abs(tbeta0[i])))

    coefmat[[i,1]] <- data.frame(Estimate=liucoef[,i],
                        `Estimate (St)`=liucoefstd[,i],
                        `StdErr (St)`=c(sebeta0[i], SEmat[,i]),
                        `t-val (St)`=c(tbeta0[i], tbeta1[,i]),
                        `Pr(>|t|)`=c(pbeta0[i], pbeta1[,i]), check.names = FALSE)
    pvals <- coefmat[[i,1]][,5]
    for(j in 1:5) coefmat[[i,1]][,j] <- round(coefmat[[i,1]][,j], digits[j])

    coefmat[[i,1]]$`Pr(>|t|)` <- paste(coefmat[[i,1]]$`Pr(>|t|)`,
                               symnum(coefmat[[i,1]]$`Pr(>|t|)`, corr = FALSE, na = FALSE,
                                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                      symbols = c("***", "** ", "*  ", ".  ", "   ")))

    cat(paste0("\nlambda = ", lambda[i],"\n"))
    cat("\nCoefficient Estimates:\n")
    print(coefmat[[i,1]], digits)
    cat("\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

    stats <- data.frame(R2=lstat$R2[i], AdjR2=lstat$AdjR2[i], Fval=lstat$FVal[i],
                        AIC=linfo[i,1], BIC=linfo[i,2], MSE=lstat$MSE[i])
    rownames(stats)<-""
    cat("\nStatistics:\n")
    print(stats, digits = 4)
    cat("----------------------------------------------------------------------\n")
  }
  invisible(coefmat)
}

