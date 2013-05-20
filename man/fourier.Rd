\name{fourier}
\alias{fourier}
\title{Fourier expansion smoothing}
\description{
Estimates a model of the form y=f(z)+XB+u using a fourier expansion for z.  
The variable z is first transformed to \eqn{z = 2\pi\frac{(z-min(z))}{max(z)-min(z))}}{z = 2 pi (z-min(z))/(max(z)-min(z))}. 
The fourier model is \eqn{y = \alpha_1 z + \alpha_2 z^2 + \sum(\lambda_q sin(qz) + \delta_q cos(qz)) + X \beta +u}{y = \alpha_1 z + \alpha_2 z^2 + \sum_{q=1}^Q(\lambda_q sin(qz) + \delta_q cos(qz)) + X \beta +u}. 
Estimation can be carried out for a fixed value of Q or for a range of Q. In the latter case, 
the function indicates the value of Q that produces the lowest value of one of the following criteria: the AIC, the Schwarz information criterion, or the gcv.
}
\usage{
fourier(form,q=1,minq=0,maxq=0,crit="gcv",data=NULL)
}

\arguments{
  \item{form }{Model formula.  The expansion is applied to the first explanatory variable.}
  \item{q }{If q is specified and minq=maxq, fits a fourier expansion with Q set to q.  
Default is q=1, which implies a model with z, z^2, sin(z), cos(z) and X as explanatory variables.}
  \item{minq }{The lower bound to search for the value of Q that minimizes crit.  minq can take any value greater than zero.  Default:  not used.}
  \item{maxq }{The upper bound to search for the value of Q that minimizes crit.  maxq must be great than or equal to minq. Default:  not used.}
  \item{crit }{The selection criterion.  Must be in quotes.  The default is the generalized cross-validation criterion, or "gcv".  Options include the Akaike information criterion, "aic", and the Schwarz criterion, "sc".  Let nreg be the number of explanatory variables in the regression and sig2 the estimated variance.  The formulas for the available crit options are   \cr
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
  \item{fourierhat }{The predicted values for z alone, normalized to have the same mean as the dependent variable.  
If no X variables are included in the regression, fourierhat = yhat.}
  \item{q }{The value of Q used in the final estimated model.}

}
\references{
Gallant, Ronald, "On the Bias in Flexible functional Forms and an Essentially Unbiased Form:  The Fourier Flexible Form," \emph{Journal of Econometrics} 15 (1981), 211-245.

Gallant, Ronald, "Unbiased Determination of Production Technologies," \emph{Journal of Econometrics} 20 (1982), 285-323.

McMillen, Daniel P. and Jonathan Domborw, "A Flexible Fourier Approach to Repeat Sales Price Indexes," \emph{Real Estate Economics} 29 (2001), 207-225.

McMillen, Daniel P., "Neighborhood Price Indexes in Chicago:  A Fourier Repeat Sales Approach," \emph{Journal of Economic Geography} 3 (2003), 57-73.

McMillen, Daniel P., "Issues in Spatial Data Analysis," \emph{Journal of Regional Science} 50 (2010), 119-141.  

}

\seealso{

  \link{cparlwr}

  \link{cubespline}

  \link{lwr}

  \link{lwrgrid}

  \link{semip}
}



\examples{
set.seed(23849103)
n = 1000
x <- runif(n,0,2*pi)
x <- sort(x)
ybase <- x - .1*(x^2) + sin(x) - cos(x) -.5*sin(2*x) + .5*cos(2*x)
sig = sd(ybase)/2
y <- ybase + rnorm(n,0,sig)

par(ask=TRUE)
plot(x,y)
lines(x,ybase,col="red")

fit <- fourier(y~x,minq=1,maxq=10)
plot(x,ybase,type="l",xlab="x",ylab="y")
lines(x,fit$yhat,col="red")
legend("topright",c("Base","Fourier"),col=c("black","red"),lwd=1)


}

\keyword{Parmametric Models}
\keyword{Series Expansions}

