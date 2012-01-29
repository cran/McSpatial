\name{qreglwr}
\alias{qreglwr}
\title{
Locally Weighted Quantile Regression
}
\description{
Estimates a model of the form \emph{y = f(x)} using locally weighted quantile regression for a set of user-provided quantiles.  
\emph{x} can include either one or two variables.  Returns estimated values, derivatives, and standard errors for both \emph{f(x)} and \emph{df(x)/dx}.
}

\usage{ 
qreglwr(form,taumat=c(.10,.25,.50,.75,.90), window=.25,bandwidth=0,
  kern="tcub", distance="Mahal",alldata=FALSE,data=NULL)
}

\arguments{
  \item{form }{Model formula}
  \item{taumat }{Vector of target quantiles.  Default:  \emph{taumat}=c(.10,.25,.50,.75,.90) }.
  \item{window }{Window size.   Default:  0.25. }
  \item{bandwidth }{
Bandwidth.   Default: not used. }
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{qreglwr} looks for the first two letters to determine which variable is latitude and which is longitude, 
so the data set must be attached first or specified using the data option; options like data$latitude will not work.  
Default:  Mahal. }
  \item{alldata }{If \emph{alldata=T}, each observation is used as a target value for \emph{x}.  
When \emph{alldata=F}, the function is estimated at a set of points chosen by the \emph{locfit} program using an adaptive decision tree approach, 
and the \emph{akima} program is used to interpolate to the full set of observations.  
Specifying \emph{alldata=T} can lead to long estimation times.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}


\value{
  \item{target}{The target points for the original estimation of the function.}
  \item{ytarget}{The matrix of predicted values of \emph{y} at the target points, by quantile.  Rows represent targets; columns are quantiles.}
  \item{dtarget1}{The matrix of estimated derivatives \emph{dy/dx1} at the target points, by quantile.  Rows represent targets; columns are quantiles.}
  \item{dtarget2}{The matrix of estimated derivatives \emph{dy/dx2} at the target points, by quantile.  Rows represent targets; columns are quantiles.  
All zeros if the model has only one explanatory variable.}
  \item{ytarget.se}{The matrix of standard errors for the predicted values of \emph{y} at the target points, by quantile.  
Rows represent targets; columns are quantiles.}
  \item{dtarget1.se}{The matrix of standard errors for the derivatives \emph{dy/dx1} at the target points, by quantile.  
Rows represent targets; columns are quantiles. }
  \item{dtarget2.se}{The matrix of standard errors for the derivatives \emph{dy/dx2} at the target points, by quantile.  
Rows represent targets; columns are quantiles.  All zeros if the model has only one explanatory variable.}
  \item{yhat}{The matrix of predicted values of \emph{y} for the full data set, by quantile.  
Dimension = n x length(taumat).}
  \item{dhat1}{The matrix of estimated derivatives \emph{dy/dx1} for the full data set, by quantile.  
Dimension = n x length(taumat).}
  \item{dhat2}{The matrix of estimated derivatives \emph{dy/dx2} for the full data set, by quantile.  
Dimension = n x length(taumat).  All zeros if the model has only one explanatory variable.}
  \item{yhat.se}{The matrix of standard errors for the predicted values of \emph{y} for the full data set, by quantile.  
Dimension = n x length(taumat).}
  \item{dhat1.se}{The matrix of standard errors for the estimated derivatives \emph{dy/dx1} for the full data set, by quantile.  
Dimension = n x length(taumat).}
  \item{dhat2.se}{The matrix of standard errors for the estimated derivatives \emph{dy/dx2} for the full data set, by quantile.  
Dimension = n x length(taumat).   All zeros if the model has only one explanatory variable.}
}
\details{
Serves as an interface to the \emph{quantreg} package.  
Uses a kernel weight function in \emph{quantreg's} "weight" option to estimate quantile regressions at a series of target values of \emph{x}.  
\emph{x} may include either one or two variables.
The target values are found using \emph{locfit's} adaptive decision tree approach.   
The predictions are then interpolated to the full set of \emph{x} values using the \emph{akima} package.  
If \emph{alldata=T}, the procedure is applied to every value of \emph{x} rather than a set of target points.  

The weights at a target value \eqn{x_0}{x0} are given by \eqn{K(\psi/h)}, 
where \eqn{\psi} is a measure of the distance between \emph{x} and \eqn{x_0}{x0} and \emph{h} is the bandwidth or window. 
When \emph{x} includes a single variable, \eqn{\psi = x-x_0.}{\psi = x-x0.}
When \emph{x} includes two variables, the method for specifying \eqn{\psi} depends on the \emph{distance} option. 
If \emph{distance="Mahal"} or \emph{distance="Euclid"}, the \emph{i}th row of the matrix \emph{X} = (x1, x2) is transformed such  
that \eqn{x_i = sqrt(x_i * V * t(x_i))}. Under the "Mahal" option, \emph{V} is the inverse of cov(\emph{X}).  
Under the \emph{"Euclid"} option, \emph{V} is the inverse of diag(cov(\emph{X})).  
By reducing x from two dimensions to one, this transformation leads again to the simple 
kernel weighting function \eqn{K((x- x_0 )/(sd(x)*h))}{K((x- x0 )/(sd(x)*h))}.  
\emph{h} is specified by the bandwidth or window options.  

The great circle formula is used to define \emph{K} when \emph{distance = "Latlong"}; in this case, the explanatory variable list must be specified as 
\emph{~latitude+longitude} (or \emph{~lo+la} or \emph{~lat+long}, etc), with the longitude and latitude variables expressed in degrees 
(e.g.,  -87.627800 and 41.881998 for one observation of longitude and latitude, respectively).   
The order in which latitude and longitude are listed does not matter and the function only looks for the 
first two letters to determine which variable is latitude and which is longitude.  
It is important to note that the great circle distance measure is left in miles rather than being standardized.  
Thus, the window option should be specified when \emph{distance = "Latlong"} or the bandwidth should be adjusted to account for the scale.  
The kernel weighting function becomes \emph{K(distance/h)} under the \emph{"Latlong"} option. 

Since \emph{qreglwr} estimates weighted quantile regressions of the dependent variable, \emph{y}, on \eqn{x-x_0}{x-x0}, 
the intercept provides an estimate of \emph{y} at \eqn{x_0}{x0} and \eqn{\beta} provides an estimate of the slope of the quantile line, 
\emph{dy/dx}, at \eqn{x_0}{x0}.  \emph{quantreg's} standard error for the intercept is stored in \emph{ytarget.se} (target points) 
and \emph{yhat.se} (all observations).  The standard errors for the slopes are stored in \emph{dtarget1.se, dtarget2.se, dhat1.se,} and \emph{dhat2.se}.

When \emph{alldata=T}, each data point in turn is used as a target point, \eqn{x_0}{x0}.  
Fixed bandwidths may prove too small if there are regions where \emph{x} is sparse.  A nearest neighbor approach is generally preferable (e.g, window=.50).  
Estimation can be very slow when \emph{alldata=T}.  
When \emph{alldata=F}, the package \emph{locfit} is used to find a good set of target points at which to evaluate the function.  
See Loader (1999, section 12.2) for a description of the algorithm used to determine the target points.  
The \emph{akima} package is then used to interpolate the coefficient estimates and standard errors.  


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

Koenker, Roger.  \emph{Quantile Regression.}  New York:  Cambridge University Press, 2005.
Chapter 7 and Appendix A.9.

Loader, Clive.  \emph{Local Regression and Likelihood}.  New York:  Springer, 1999. 

}

\seealso{
  \code{\link{lwr}}
}


\examples{
data(cookdata)
cookdata <- cookdata[cookdata$CHICAGO==1,]
cookdata$obs <- seq(1:nrow(cookdata))
cookdata <- cookdata[cookdata$CHICAGO==1&cookdata$POPULATION>0,]
par(ask=TRUE)

# lndens = f(dcbd)
fit <- lwr(LNDENS~DCBD,window=.20,data=cookdata)
fit1 <- qreglwr(LNDENS~DCBD,taumat=c(.10,.50,.90),window=.30,kern="rect",data=cookdata)
o <- order(cookdata$DCBD)
ymin = min(fit1$yhat)
ymax = max(fit1$yhat)
plot(cookdata$DCBD[o], fit$yhat[o], type="l", ylim=c(ymin,ymax),
  xlab="Distance to CBD", ylab="Log of Population Density")
lines(cookdata$DCBD[o], fit1$yhat[o,1], col="red", lty="dashed")
lines(cookdata$DCBD[o], fit1$yhat[o,2], col="red")
lines(cookdata$DCBD[o], fit1$yhat[o,3], col="red", lty="dashed")
legend("topright", c("LWR", "tau = 50", "tau = 10, 90"), col=c("black","red", "red"), 
  lwd=1, lty=c("solid","solid","dashed"))

\dontrun{
library(maptools)
library(RColorBrewer)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1,]
# lndens = f(longitude, latitude), weights are function of straight-line distance
fit <- qreglwr(LNDENS~LONGITUDE+LATITUDE,taumat=c(.10,.50,.90),window=.20,data=cookdata)
cmap$lwr10[cookdata$obs] <- fit$yhat[,1]
cmap$lwr50[cookdata$obs] <- fit$yhat[,2]
cmap$lwr90[cookdata$obs] <- fit$yhat[,3]
cmap$lwr1090[cookdata$obs] <- fit$yhat[,3] - fit$yhat[,1]
brks <- seq(min(cmap$lwr10,na.rm=TRUE),max(cmap$lwr10,na.rm=TRUE),length=9)
spplot(cmap,"lwr10",at=brks,col.regions=rev(brewer.pal(9,"RdBu")),
   main="Log Density Estimates, tau = .10")
brks <- seq(min(cmap$lwr50,na.rm=TRUE),max(cmap$lwr50,na.rm=TRUE),length=9)
spplot(cmap,"lwr50",at=brks,col.regions=rev(brewer.pal(9,"RdBu")),
   main="Log Density Estimates, tau = .50")
brks <- seq(min(cmap$lwr90,na.rm=TRUE),max(cmap$lwr90,na.rm=TRUE),length=9)
spplot(cmap,"lwr90",at=brks,col.regions=rev(brewer.pal(9,"RdBu")),
   main="Log Density Estimates, tau = .90")
brks <- seq(min(cmap$lwr1090,na.rm=TRUE),max(cmap$lwr1090,na.rm=TRUE),length=9)
spplot(cmap,"lwr1090",at=brks,col.regions=rev(brewer.pal(9,"RdBu")),
   main="Difference in Log Density, tau = .90 - .10")
}
}

\keyword{Quantile Regression}
\keyword{Nonparametric}

