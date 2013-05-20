\name{lwrgrid}
\alias{lwrgrid}
\title{LWR Bandwidth or Window Selection}
\description{
Finds the value of a user-provided array of window or bandwidth values that provides the lowest \emph{cv} or \emph{gcv} for an LWR model.  
Calls \emph{lwr} and returns its full output for the chosen value of \emph{h}.
}

\usage{ 
lwrgrid(form,window=0,bandwidth=0,kern="tcub",method="gcv",print=TRUE,
  distance="Mahal",target=NULL,data=NULL)
}

\arguments{
  \item{form }{Model formula}
  \item{window }{Vector of possible window sizes.   Default:  none. }
  \item{bandwidth }{Vector of possible bandwidths.  Default:  none.}
  \item{method }{Specifies "gcv" or "cv" criterion function.  Default:  method="gcv".}
  \item{print }{If TRUE, prints \emph{gcv} or \emph{cv}  values for each value of the window or bandwidth.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{lwr} looks for the first two letters to determine which variable is latitude and which is longitude, 
so data set must be attached first or specified using the data option; options like data$latitude will not work.  Default:  Mahal. }
  \item{target}{If \emph{target = NULL}, uses the \emph{maketarget} command to form targets using \emph{alpha = window}.  
In this case, the target points are found using a nearest neighbor procedure even if \emph{bandwidth} > 0. 
If \emph{target="alldata"}, each observation is used as a target value for \emph{x}.  
A set of target values can be supplied directly.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}


\value{
  \item{target}{The target points for the original estimation of the function.}
  \item{ytarget}{The predicted values of \emph{y} at the original target points.}
  \item{dtarget1}{The estimated derivatives \emph{dy/dx1} at the target points.}
  \item{dtarget2}{The estimated derivatives \emph{dy/dx2} at the target points.  All zeros if the model has only one explanatory variable.}
  \item{ytarget.se}{Standard errors for the predicted values of \emph{y} at the target points.}
  \item{dtarget1.se}{Standard errors for the derivatives \emph{dy/dx1} at the target points. }
  \item{dtarget2.se}{Standard errors for the derivatives \emph{dy/dx2} at the target points.  All zeros if the model has only one explanatory variable.}
  \item{yhat}{The predicted values of \emph{y} for the full data set.}
  \item{dhat1}{The estimated derivatives \emph{dy/dx1} for the full data set.}
  \item{dhat2}{The estimated derivatives \emph{dy/dx2} for the full data set.  All zeros if the model has only one explanatory variable.}
  \item{yhat.se}{Standard errors for the predicted values of \emph{y} for the full data set.}
  \item{dhat1.se}{Standard errors for the estimated derivatives \emph{dy/dx1} for the full data set.}
  \item{dhat2.se}{Standard errors for the estimated derivatives \emph{dy/dx2} for the full data set.   All zeros if the model has only one explanatory variable.}
  \item{df1}{\emph{tr(L)}, a measure of the degrees of freedom used in estimation.}
  \item{df2}{\emph{tr(L'L)}, an alternative measure of the degrees of freedom used in estimation.  }
  \item{sig2}{Estimated residual variance, \emph{sig2 = rss/(n-2*df1+df2)}.}
  \item{cv}{Cross-validation measure.  \emph{cv = mean(((y-yhat)/(1-infl))^2)}, 
where \emph{yhat} is vector of predicted values for \emph{y} and \emph{infl} is the vector of diagonal terms for \emph{L}.}
  \item{gcv}{\emph{gcv = n*(n*sig2)/((n-nreg)^2)}, where \emph{sig2} is the estimated residual variance and \emph{nreg = 2*df1 - df2.}}
  \item{infl}{A vector containing the diagonal elements of \emph{L}.}
  \item{minh}{Value of window or bandwidth the minimizes the criterion function.}
}

\references{
Cleveland, William S. and Susan J. Devlin, "Locally Weighted Regression:  An Approach to Regression Analysis by Local Fitting," \emph{Journal of the American Statistical Association} 83 (1988), 596-610. 

Loader, Clive.  \emph{Local Regression and Likelihood}.  New York:  Springer, 1999. 

McMillen, Daniel P., "Issues in Spatial Data Analysis," \emph{Journal of Regional Science} 50 (2010), 119-141.  

McMillen, Daniel P., "Employment Densities, Spatial Autocorrelation, and Subcenters in Large Metropolitan Areas," \emph{Journal of Regional Science} 44 (2004), 225-243. 

McMillen, Daniel P. and John F. McDonald, "A Nonparametric Analysis of Employment Density in a Polycentric City," \emph{Journal of Regional Science} 37 (1997), 591-612. 

McMillen, Daniel P. and Christian Redfearn, "Estimation and Hypothesis Testing for Nonparametric Hedonic House Price Functions," \emph{Journal of Regional Science} 50 (2010), 712-733.

Pagan, Adrian and Aman Ullah.  \emph{Nonparametric Econometrics}.  New York:  Cambridge University Press, 1999. 

Silverman, A. W., \emph{Density Estimation for Statistics and Data Analysis}, Chapman and Hall, New York (1986). 

}

\seealso{
  \link{lwr}

  \link{maketarget}
}

\examples{
par(ask=TRUE)
n = 1000
z1 <- runif(n,0,2*pi)
z1 <- sort(z1)
z2 <- runif(n,0,2*pi)
o1 <- order(z1)
o2 <- order(z2)
ybase1 <-  z1 - .1*(z1^2) + sin(z1) - cos(z1) - .5*sin(2*z1) + .5*cos(2*z1) 
ybase2 <- -z2 + .1*(z2^2) - sin(z2) + cos(z2) + .5*sin(2*z2) - .5*cos(2*z2)
ybase <- ybase1+ybase2
sig = sd(ybase)/2
y <- ybase + rnorm(n,0,sig)
summary(lm(y~ybase))

# Single variable estimation
fit1 <- lwrgrid(y~z1,window=seq(.10,.30,.05))
c(fit1$df1,fit1$df2,2*fit1$df1-fit1$df2)
plot(z1[o1],ybase1[o1],type="l",ylim=c(min(ybase1,fit1$yhat),max(ybase1,fit1$yhat)),
  xlab="z1",ylab="y")
# Make predicted and actual values have the same means
fit1$yhat <- fit1$yhat - mean(fit1$yhat) + mean(ybase1)
lines(z1[o1],fit1$yhat[o1], col="red")
legend("topright", c("Base", "LWR"), col=c("black","red"),lwd=1)
fit2 <- lwrgrid(y~z2,window=seq(.10,.90,.10))
fit2$yhat <- fit2$yhat - mean(fit2$yhat) + mean(ybase2)
c(fit2$df1,fit2$df2,2*fit2$df1-fit2$df2)
plot(z2[o2],ybase2[o2],type="l",ylim=c(min(ybase2,fit2$yhat),max(ybase2,fit2$yhat)),
    xlab="z1",ylab="y")
lines(z2[o2],fit2$yhat[o2], col="red")
legend("topright", c("Base", "LWR"), col=c("black","red"),lwd=1)

#both variables
fit3 <- lwrgrid(y~z1+z2,window=seq(.03,.09,.02))
yhat1 <- fit3$yhat - mean(fit3$yhat) + mean(ybase1)
plot(z1[o1],yhat1[o1], xlab="z1",ylab="y")
lines(z1[o1],ybase1[o1],col="red")
yhat2 <- fit3$yhat - mean(fit3$yhat) + mean(ybase2)
plot(z2[o2],yhat2[o2], xlab="z2",ylab="y")
lines(z2[o2],ybase2[o2],col="red")


}


\keyword{Locally Weighted Regression}
\keyword{Nonparametric}

