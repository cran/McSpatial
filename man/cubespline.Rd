\name{cubespline}
\alias{cubespline}
\title{
Smooth cubic spline estimation 
}
\description{
Estimates a smooth cubic spline for a model of the form y = f(z) + XB +u.  
The function divides the range of z into knots+1 equal intervals.  
The regression is 
\eqn{y=cons+\lambda_1 (z-z_0) + \lambda_2 (z-z_0)^2 + \lambda_3 (z-z_0)^3 + \sum_{k=1}^{K} \gamma_k(z-z_k)^3 D_k + X \beta + u}{y=cons+\lambda_1 (z-z_0) + \lambda_2 (x-x_0)^2 + \lambda_3 (x-x_0)^3 + \sum_{k=1}^{K} \gamma_k(x-x_k)^3 D_k + X \beta + u} 
where \eqn{z_0 = min(z), z_1...z_K}{z_0 = min(z), z_1...z_K} 
are the knots, and \eqn{D_k=1}{D_k=1} if \eqn{z \geq z_k}{z >= z_k}.  
Estimation can be carried out for a fixed value of K or for a range of K.  
In the latter case, the function indicates the value of K that produces the lowest value of one of the following criteria:  
the AIC, the Schwarz information criterion, or the gcv.  
}


\usage{
cubespline(form,knots=1,mink=1,maxk=1,crit="gcv",data=NULL)
}

\arguments{
  \item{form }{Model formula.  The spline is used with the first explanatory variable.}
  \item{knots }{If knots is specified, fits a cubic spline with K = knots.  
Default is knots=1, mink = 1, and maxk = 1, which implies a cubic spline with a single knot.}
  \item{mink }{The lower bound to search for the value of K that minimizes crit.  mink can take any value greater than zero.  The default is mink = 1}
  \item{maxk }{The upper bound to search for the value of K that minimizes crit.  maxk must be great than or equal to mink. The default is maxk = 1.}
  \item{crit }{The selection criterion.  Must be in quotes.  The default is the generalized cross-validation criterion, or "gcv".  
Options include the Akaike information criterion, "aic", and the Schwarz criterion, "sc".  Let nreg be the number of explanatory variables in the regression and sig2 the estimated variance.  The formulas for the available crit options are   \cr
	gcv = n*(n*sig2)/((n-nreg)^2) \cr
	aic = log(sig2) + 2* nreg /n \cr
	sc  = log(sig2) + log(n)*nreg /n}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}


\value{
  \item{yhat }{The predicted values of the dependent variable at the original data points}
  \item{rss }{The residual sum of squares}
  \item{sig2 }{The estimated error variance}
  \item{aic }{The value for AIC}
  \item{sc }{The value for sc}
  \item{gcv }{The value for gcv}
  \item{coef }{The estimated coefficient vector, B}
  \item{splinehat }{The predicted values for z alone, normalized to have the same mean as the dependent variable.  
If no X variables are included in the regression, splinehat = yhat.}
  \item{knots }{The vector of knots}

}
\references{
McMillen, Daniel P., "Testing for Monocentricity," in Richard J. Arnott and Daniel P. McMillen, eds., \emph{A Companion to Urban Economics}, Blackwell, Malden MA (2006), 128-140.

McMillen, Daniel P., "Issues in Spatial Data Analysis," \emph{Journal of Regional Science} 50 (2010), 119-141.  

Suits, Daniel B., Andrew Mason, and Louis Chan, "Spline Functions Fitted by Standard Regression Methods," \emph{Review of Economics and Statistics} 60 (1978), 132-139.

}

\seealso{
  \link{cparlwr}

  \link{fourier}

  \link{lwr}

  \link{lwrgrid}

  \link{semip}
}


\examples{
data(cookdata)
fardata <- cookdata[!is.na(cookdata$LNFAR),]
par(ask=TRUE)

# single variable
o <- order(fardata$DCBD)
fit1 <- cubespline(LNFAR~DCBD, mink=1, maxk=10,data=fardata)
c(fit1$rss, fit1$sig2, fit1$aic, fit1$sc, fit1$gcv, fit1$knots)
plot(fardata$DCBD[o], fardata$LNFAR[o], xlab="Distance from CBD", ylab="Log FAR")
lines(fardata$DCBD[o], fit1$splinehat[o], col="red")

# multiple explanatory variables
fit2 <- cubespline(fardata$LNFAR~fardata$DCBD+fardata$AGE,  mink=1, maxk=10)
c(fit2$rss, fit2$sig2, fit2$aic, fit2$sc, fit2$gcv, fit2$knots)
plot(fardata$DCBD[o], fardata$LNFAR[o], xlab="Distance from CBD", ylab="Log FAR")
lines(fardata$DCBD[o], fit2$splinehat[o], col="red")

# pre-specified number of knots
fit3 <- cubespline(LNFAR~DCBD+AGE,  knots=4, data=fardata)

}

\keyword{Parmametric Models}
\keyword{Series Expansions}

