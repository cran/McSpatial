\name{cparlwrgrid}
\alias{cparlwrgrid}
\title{
Conditionally parametric LWR regression bandwidth or window selection
}
\description{
Finds the value of a user-provided array of window or bandwidth values that provides the lowest \emph{cv} or \emph{gcv} for a \emph{CPAR} model.  
Calls \emph{cparlwr} and returns its full output for the chosen value of \emph{h}.
}
\usage{ 
cparlwrgrid(form,nonpar,window=0,bandwidth=0,kern="tcub",method="gcv",
  print=TRUE,distance="Mahal",alldata=FALSE,data=NULL)
}

\arguments{
  \item{form }{Model formula}
  \item{nonpar }{List of either one or two variables for \emph{z}.  
Formats: \emph{cparlwr(y~xlist, nonpar=~z1, ...)} or \emph{cparlwr(y~xlist, nonpar=~z1+z2, ...)}.  
Important:  note the "~" before the first \emph{z} variable.
}
  \item{window }{Window size.   Default:  not used. }
  \item{bandwidth }{Bandwidth.   Default:  not used.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{method }{Specifies "gcv" or "cv" criterion function.  Default:  method="gcv".}
  \item{print }{If TRUE, prints \emph{gcv} or \emph{cv}  values for each value of the window or bandwidth.}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{cparlwr} looks for the first two letters to determine which variable is latitude and which is longitude, 
so the data set must be attached first or specified using the data option; options like data$latitude will not work.  Default:  Mahal. }
  \item{alldata }{If \emph{alldata=T}, each observation is used as a target value for \emph{z}.  
When \emph{alldata=F}, the function is estimated at a set of points chosen by the \emph{locfit} program using an adaptive decision tree approach, 
and the \emph{akima} program is used to interpolate to the full set of observations.  Specifying \emph{alldata=T} can lead to long estimation times.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}

\value{
  \item{target}{The target points for the original estimation of the function.}
  \item{ytarget}{The predicted values of \emph{y} at the target values \emph{z}.}
  \item{xcoef.target}{Estimated coefficients, \emph{B(z)}, at the target values of \emph{z}.}
  \item{ytarget.se}{Standard errors for the predicted \emph{y} at the target values of \emph{z}.}
  \item{xcoef.target.se}{Standard errors for \emph{B(z)} at the target values of \emph{z}.}
  \item{yhat}{Predicted values of \emph{y} at the original data points.}
  \item{xcoef}{Estimated coefficients, \emph{B(z)}, at the original data points.}
  \item{yhat.se}{Standard errors for the predicted values of \emph{y}, full data set.}
  \item{xcoef.se}{Standard errors for \emph{B(z)} with \emph{z} evaluated at all points in the data set.}
  \item{df1}{\emph{tr(L)}, a measure of the degrees of freedom used in estimation.}
  \item{df2}{\emph{tr(L'L)}, an alternative measure of the degrees of freedom used in estimation.  }
  \item{sig2}{Estimated residual variance, \emph{sig2 = rss/(n-2*df1+df2)}.}
  \item{cv}{Cross-validation measure.  \emph{cv = mean(((y-yhat)/(1-infl))^2) }, 
where \emph{yhat} is the vector of predicted values for \emph{y} and \emph{infl} is the vector of diagonal terms for \emph{L}. }
  \item{gcv}{gcv = n*(n*sig2)/((n-nreg)^2), where \emph{sig2} is the estimated residual variance and \emph{nreg = 2*df1 - df2}.}
  \item{infl}{A vector containing the diagonal elements of \emph{L}.}
}

\references{
Cleveland, William S. and Susan J. Devlin, "Locally Weighted Regression:  An Approach to Regression Analysis by Local Fitting," 
\emph{Journal of the American Statistical Association} 83 (1988), 596-610. 

Loader, Clive. \emph{ Local Regression and Likelihood.}  New York:  Springer, 1999. 

McMillen, Daniel P., "One Hundred Fifty Years of Land Values in Chicago:  A Nonparametric Approach," \emph{Journal of Urban Economics} 40 (1996), 100-124. 

McMillen, Daniel P., "Issues in Spatial Data Analysis," \emph{Journal of Regional Science} 50 (2010), 119-141.  

McMillen, Daniel P., "Employment Densities, Spatial Autocorrelation, and Subcenters in Large Metropolitan Areas," \emph{Journal of Regional Science} 44 (2004), 225-243. 

McMillen, Daniel P. and John F. McDonald, "A Nonparametric Analysis of Employment Density in a Polycentric City," \emph{Journal of Regional Science} 37 (1997), 591-612. 

McMillen, Daniel P. and Christian Redfearn, ``Estimation and Hypothesis Testing for Nonparametric Hedonic House Price Functions,'' 
\emph{Journal of Regional Science} 50 (2010), 712-733.

Pagan, Adrian and Aman Ullah.  \emph{Nonparametric Econometrics.}  New York:  Cambridge University Press, 1999.

}

\seealso{\link{cparlwr}}

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
fit1 <- cparlwrgrid(y~z1,nonpar=~z1,window=seq(.10,.40,.10))
c(fit1$df1,fit1$df2,2*fit1$df1-fit1$df2)
plot(z1[o1],ybase1[o1],type="l",ylim=c(min(ybase1,fit1$yhat),max(ybase1,fit1$yhat)),
  xlab="z1",ylab="y")
# Make predicted and actual values have the same means
fit1$yhat <- fit1$yhat - mean(fit1$yhat) + mean(ybase1)
lines(z1[o1],fit1$yhat[o1], col="red")
legend("topright", c("Base", "LWR"), col=c("black","red"),lwd=1)
fit2 <- cparlwrgrid(y~z2,nonpar=~z2,window=seq(.10,.40,.10))
fit2$yhat <- fit2$yhat - mean(fit2$yhat) + mean(ybase2)
c(fit2$df1,fit2$df2,2*fit2$df1-fit2$df2)
plot(z2[o2],ybase2[o2],type="l",ylim=c(min(ybase2,fit2$yhat),max(ybase2,fit2$yhat)),
    xlab="z1",ylab="y")
lines(z2[o2],fit2$yhat[o2], col="red")
legend("topright", c("Base", "LWR"), col=c("black","red"),lwd=1)

#both variables
fit3 <- cparlwrgrid(y~z1+z2,nonpar=~z1+z2,window=seq(.10,.20,.05))
yhat1 <- fit3$yhat - mean(fit3$yhat) + mean(ybase1)
plot(z1[o1],yhat1[o1], xlab="z1",ylab="y")
lines(z1[o1],ybase1[o1],col="red")
yhat2 <- fit3$yhat - mean(fit3$yhat) + mean(ybase2)
plot(z2[o2],yhat2[o2], xlab="z2",ylab="y")
lines(z2[o2],ybase2[o2],col="red")
}


\keyword{Conditionally Parametric}
\keyword{Nonparametric}
\keyword{Locally Weighted Regression}

