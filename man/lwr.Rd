\name{lwr}
\alias{lwr}
\title{
Locally Weighted Regression
}
\description{
Estimates a model of the form \emph{y = f(x)} using locally weighted regression.  
\emph{x} can include either one or two variables.  Returns estimated values, derivatives, and standard errors for both \emph{f(x)} and \emph{df(x)/dx}.
}

\usage{ 
lwr(form,window=.25,bandwidth=0,kern="tcub",distance="Mahal",
  target=NULL,data=NULL)
}

\arguments{
  \item{form}{Model formula}
  \item{window}{Window size.   Default:  0.25.}
  \item{bandwidth}{Bandwidth.   Default: not used.}
  \item{kern}{Kernel weighting function.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance}{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{lwr} looks for the first two letters to determine which variable is latitude and which is longitude, 
so the data set must be attached first or specified using the data option; options like data$latitude will not work.  
Default:  Mahal. }
  \item{target}{If \emph{target = NULL}, uses the \emph{maketarget} command to form targets using the values specified for \emph{window}, 
\emph{bandwidth}, and \emph{kern}.  If \emph{target="alldata"}, each observation is used as a target value for \emph{x}.  
A set of target values can be supplied directly.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
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
}
\details{
The estimated value of \emph{y} at a target value \eqn{x_0}{x0} is the predicted value from a weighted least squares regression of \emph{y} on \eqn{x-x_0}{x-x0} 
with weights given by \eqn{K(\psi/h)}, where \eqn{\psi} is a measure of the distance between \emph{x} and \eqn{x_0}{x0} and \emph{h} is the bandwidth or window. 


When \emph{x} includes a single variable, \eqn{\psi = x-x_0.}{\psi = x-x0.}
When \emph{x} includes two variables, the method for specifying \eqn{\psi} depends on the \emph{distance} option. 
If \emph{distance="Mahal"} or \emph{distance="Euclid"}, the \emph{i}th row of the matrix \emph{X} = (x1, x2) is transformed such  
that \eqn{x_i} = sqrt(\emph{x_i * V * t(x_i)}). Under the "Mahal" option, \emph{V} is the inverse of cov(\emph{X}).  
Under the \emph{"Euclid"} option, \emph{V} is the inverse of diag(cov(\emph{X})).  
By reducing x from two dimensions to one, this transformation leads again to the simple 
kernel weighting function \eqn{K((x- x_0 )/(sd(x)*h))}{K((x- x0 )/(sd(x)*h))}.  


The great circle formula is used to define \emph{K} when \emph{distance = "Latlong"}; in this case, the explanatory variable list must be specified as 
\emph{~latitude+longitude} (or \emph{~lo+la} or \emph{~lat+long}, etc), with the longitude and latitude variables expressed in degrees 
(e.g.,  -87.627800 and 41.881998 for one observation of longitude and latitude, respectively).   
The order in which latitude and longitude are listed does not matter and the function only looks for the 
first two letters to determine which variable is latitude and which is the longitude.  
It is important to note that the great circle distance measure is left in miles rather than being standardized.  
Thus, the window option should be specified when \emph{distance = "Latlong"} or the bandwidth should be adjusted to account for the scale.  
The kernel weighting function becomes \emph{K(distance/h)} under the \emph{"Latlong"} option. 


\emph{h} is specified by the bandwidth or window options.  
The intercept, \eqn{\alpha}, provides an estimate of \emph{y} at \eqn{x_0}{x0} and \eqn{\beta} provides an estimate of the slope, \emph{dy/dx} at \eqn{x_0}{x0}.  
When \emph{target="alldata"}, each data point in turn is used as a target point, \eqn{x_0}{x0}.  


Since each estimate is a linear function of all \emph{n} values for \emph{y}, 
the full set of estimates takes the form \eqn{\hat{y} = LY}{yhat = LY}, where \emph{L} is an \emph{nxn} matrix.
Loader (1999) suggests two measures of the number of degrees of freedom used in estimation, 
\emph{df1 = tr(L)} and \emph{df2 = tr(L'L)}, both of which are stored by \emph{lwr}.  
The diagonal elements of \emph{tr(L)} are stored in the array \emph{infl}.  
Again following Loader (1999), the degrees of freedom correction used to estimate the error variance, 
\eqn{\hat{\sigma^2}}{sig2}, is \emph{df = 2*df1 -  df2}.  Let \emph{e} represent the vector of residuals, \eqn{e = y-\hat{y}}{e = y-yhat}.  
The estimated variance is \eqn{\hat{\sigma^2} = \sum_ie_i^2/(n-df)}{sig2 = sum(e^2)/(n-df)}. The covariance matrix is 
\deqn{\hat{\sigma}^2(\sum_{i=1}^n Z_i K(\psi_i/h) Z_i^\top)^{-1}(\sum_{i=1}^n Z_i (K(\psi_i/h))^2 Z_i^\top )(\sum_{i=1}^n Z_i K(\psi_i/h) Z_i^\top)^{-1}.}{\sigma^2(\sum Z_i K(\psi_i/h) Z_i')^{-1}(\sum Z_i (K(\psi_i/h))^2 Z_i' )(\sum Z_i K(\psi_i/h) Z_i')^{-1}.}
where \eqn{Z = (1 \; x-x_0).}{Z = (1 x-x0).}  


Estimation can be very slow when \emph{target = "alldata"}.  The \emph{maketarget} command can be used to identify target points. 
The \emph{smooth12} command is then used to interpolate the coefficient estimates, the standard errors, 
and the values used to form \emph{df1} and \emph{df2}.  


\emph{h} can be specified to be either a fixed bandwidth or a window size set to a percentage of the sample size.  
Optionally, the \emph{lwrgrid} command can be used to specify a vector of values for \emph{h} with \emph{lwr} picking the one that minimizes a criterion function.  
In general, the \emph{window} option will be preferable because it provides more accurate estimates in regions where \emph{x} is relatively sparse. 


Available kernel weighting functions include the following: 

\tabular{lll}{
Kernel \tab Call abbreviation \tab Kernel function K(z) \cr
Rectangular \tab ``rect'' \tab \eqn{\frac{1}{2} I(|z| <1)}{1/2 * I(|z|<1)} \cr
Triangular \tab ``tria'' \tab \eqn{(1-|z|)I(|z|<1)}{(1-|z|) * I(|z|<1)}\cr
Epanechnikov \tab ``epan'' \tab \eqn{\frac{3}{4} (1-z^2) * I(|z| <1)}{3/4 * (1-z^2)*I(|z| < 1)} \cr
Bi-Square \tab ``bisq'' \tab \eqn{\frac{15}{16} (1-z^2)^2 * I(|z| <1)}{15/16 * (1-z^2)^2 * I(|z| < 1)} \cr
Tri-Cube \tab ``tcub'' \tab \eqn{\frac{70}{81} (1-|z|^3)^3 * I(|z| <1)}{70/81 * (1-|z|^3)^3 * I(|z| < 1)} \cr
Tri-Weight \tab ``trwt'' \tab \eqn{\frac{35}{32} (1-z^2)^3 * I(|z| <1)}{35/32 * (1-z^2)^3 * I(|z| < 1)} \cr
Gaussian \tab ``gauss'' \tab \eqn{(2\pi)^{-.5} e^{-z^2/2}}{2pi^{-.5} exp(-z^2/2)} \cr }
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
  \link{cparlwr}

  \link{cubespline}

  \link{fourier}

  \link{lwrgrid}

  \link{maketarget}

  \link{semip}
}

\examples{
# 1. Monte Carlo data
n = 1000
x <- runif(n,0,2*pi)
x <- sort(x)
ybase <- x - .1*(x^2) + sin(x) - cos(x) -.5*sin(2*x) + .5*cos(2*x)
sig = sd(ybase)/2
y <- ybase + rnorm(n,0,sig)
par(ask=TRUE)
plot(x,y)
lines(x,ybase,col="red")
fit <- lwr(y~x, window=.15)
# plot 95% confidence intervals for predicted y 
predse <- sqrt(fit$sig2 + fit$yhat.se^2)
lower <- fit$yhat + qnorm(.025)*predse
upper <- fit$yhat + qnorm(.975)*predse
plot(x, ybase, type="l", ylim=c(min(lower), max(upper)), 
  main="Estimated Function", xlab="x", ylab="y")
lines(x, fit$yhat, col="red")
lines(x, lower, lty="dashed", col="red")
lines(x, upper, lty="dashed", col="red")
legend("topleft", c("Base", "Predicted", "95 Percent CI"), 
 col=c("black", "red", "red"), lty=c("solid", "solid", "dashed"), lwd=1)

# plot 95%  confidence intervals for slopes
dxbase <- 1 - .2*x + cos(x) + sin(x) - cos(2*x) - sin(2*x)
lower <- fit$dhat1 + qnorm(.025)*fit$dhat1.se
upper <- fit$dhat1 + qnorm(.975)*fit$dhat1.se
plot(x, dxbase, type="l", ylim=c(min(lower), max(upper)), 
  main="Estimated Slopes", xlab="x", ylab="y")
lines(x, fit$dhat1, col="red")
lines(x, lower, lty="dashed", col="red")
lines(x, upper, lty="dashed", col="red")
legend("topright", c("Base", "Predicted", "95 Percent CI"), 
 col=c("black", "red", "red"),lty=c("solid", "solid", "dashed"), lwd=1)

# Derivative estimates with larger window size
fit <- lwr(y~x,window=.20)
lower <- fit$dhat1 + qnorm(.025)*fit$dhat1.se
upper <- fit$dhat1 + qnorm(.975)*fit$dhat1.se
plot(x, dxbase, type="l", ylim=c(min(lower), max(upper)), 
  main="Estimated Slopes", xlab="x", ylab="y")
lines(x, fit$dhat1, col="red")
lines(x, lower, lty="dashed", col="red")
lines(x, upper, lty="dashed", col="red")
legend("topright", c("Base", "Predicted", "95 Percent CI"), 
 col=c("black", "red", "red"), lty=c("solid", "solid", "dashed"), lwd=1)

\dontrun{
#2. Population density data
library(RColorBrewer)

cook <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cook$obs <- seq(1:nrow(cook))
# measure distance to Chicago city center
lmat <- coordinates(cook)
cook$LONGITUDE <- lmat[,1]
cook$LATITUDE  <- lmat[,2]
cook$DCBD <- geodistance(longvar=cook$LONGITUDE,latvar=cook$LATITUDE,
  lotarget=-87.627800,latarget=41.881998,dcoor=FALSE)$dist
# population density = population/acres,  acres = square mile x 640
cook$LNDENS <- log(cook$POPULATION/(cook$AREA*640))
densdata <- data.frame(cook[cook$POPULATION>0,])
par(ask=TRUE)

# lndens = f(longitude, latitude), weights are function of straight-line distance
fit <- lwr(LNDENS~LONGITUDE+LATITUDE,  window=.10, 
   distance="Latlong",data=densdata)
c(fit$df1, fit$df2, 2*fit$df1-fit$df2)
cook$lwrhat[densdata$obs] <- fit$yhat
brks <- seq(min(cook$lwrhat,na.rm=TRUE),max(cook$lwrhat,na.rm=TRUE),length=9)
spplot(cook,"lwrhat",at=brks,col.regions=rev(brewer.pal(9,"RdBu")),
   main="Log Density LWR Estimates")

}
}

\keyword{Locally Weighted Regression}
\keyword{Nonparametric}


