\name{qregcpar}
\alias{qregcpar}
\title{
Conditionally Parametric LWR Quantile Estimation 
}
\description{
Estimates a model of the form \emph{y = XB(z) + u} using locally weighted quantile regression for a set of user-provided quantiles.
\emph{z} can include one or two variables.  

}

\usage{ 
qregcpar(form,nonpar,taumat=c(.10,.25,.50,.75,.90),
  window=.25,bandwidth=0,kern="tcub",distance="Mahal",
  target=NULL,data=NULL) 
}

\arguments{
  \item{form }{Model formula}
  \item{nonpar }{List of either one or two variables for \emph{z}.  
Formats: \emph{qregcpar(y~xlist, nonpar=~z1, ...)} or \emph{qregcpar(y~xlist, nonpar=~z1+z2, ...).}
Important:  note the "~" before the first \emph{z} variable.
}
  \item{taumat }{Vector of target quantiles.  Default:  \emph{taumat}=c(.10,.25,.50,.75,.90) }.
  \item{window }{Window size.   Default:  0.25. }
  \item{bandwidth }{Bandwidth.   Default:  not used.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{qregcpar} looks for the first two letters to determine which variable is latitude and which is longitude, 
so the data set must be attached first or specified using the data option; options like data$latitude will not work.  Default:  Mahal. }
  \item{target}{If \emph{target = NULL}, uses the \emph{maketarget} command to form targets using the values specified for \emph{window}, 
\emph{bandwidth}, and \emph{kern}.  If \emph{target="alldata"}, each observation is used as a target value for \emph{x}.  
A set of target values can be supplied directly. }

  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}

\value{
  \item{target}{The target points for the original estimation of the function.}
  \item{xcoef.target}{The matrix of estimated coefficients, \emph{B(z)}, at the target values of \emph{z}.  Dimension = \emph{ntarget x ntau x k}, 
where \emph{ntarget} = number of target poitns, \emph{ntau} = number of quantiles, and \emph{k} = number of explanatory variables including the intercept.}
  \item{xcoef.target.se}{The matrix of standard errors for \emph{B(z)} at the target values of \emph{z}.  Dimension = \emph{ntarget x ntau x k.}}
  \item{xcoef}{The matrix of estimated coefficients, \emph{B(z)}, at the original data points.  Dimension = \emph{n x ntau x k.}}
  \item{xcoef.se}{The matrix of standard errors for \emph{B(z)} with \emph{z} evaluated at all points in the data set.  Dimension = \emph{n x ntau x k.}}
  \item{yhat}{The matrix of predicted values of \emph{y} at the original data points.  Dimension = \emph{n x ntau.}}
}

\details{
The list of explanatory variables is specified in the base model formula while \emph{Z} is specified using \emph{nonpar}.  
\emph{X} can include any number of explanatory variables, but \emph{Z} must have at most two.  

The estimated value of \emph{y} at a target value \eqn{z_0} and a quantile \eqn{\tau} is the predicted value from a weighted quantile regression 
of \emph{y} on \emph{X} with weights given by \emph{K}.
When \emph{Z} includes a single variable, \emph{K} is a simple kernel weighting function:  \eqn{K((z - z_0 )/(sd(z)*h)) }.  
When \emph{Z} includes two variables (e.g, nonpar=~z1+z2), the method for specifying \emph{K} depends on the \emph{distance} option. 
Under either option, the \emph{i}th row of the matrix \emph{Z} = (z1, z2) is transformed such  
that \eqn{z_i = \sqrt{z_i * V * t(z_i)}}{z_i = sqrt(z_i * V * t(z_i))}. Under the "Mahal" option, \emph{V} is the inverse of cov(\emph{Z}).  
Under the \emph{"Euclid"} option, \emph{V} is the inverse of diag(cov(\emph{Z})).  
After this transformation, the weights again reduce to the simple kernel weighting function \eqn{K((z- z_0 )/(sd(z)*h))}.  

The great circle formula is used to define \emph{K} when \emph{distance = "Latlong"}; in this case, the variable list for \emph{nonpar} must be listed as 
\emph{nonpar = ~latitude+longitude} (or \emph{~lo+la} or \emph{~lat+long}, etc), with the longitude and latitude variables expressed in degrees 
(e.g.,  -87.627800 and 41.881998 for one observation of longitude and latitude, respectively).   
The order in which latitude and longitude are listed does not matter and the function only looks for the 
first two letters to determine which variable is latitude and which is the longitude.  
It is important to note that the great circle distance measure is left in miles rather than being standardized.  
Thus, the window option should be specified when \emph{distance = "Latlong"} or the bandwidth should be adjusted to account for the scale.  
The kernel weighting function becomes \emph{K(distance/h)} under the \emph{"Latlong"} option. 
\emph{h} is specified by the \emph{bandwidth} or \emph{window} option.  

For each quantile, the estimated coefficient matrix, \emph{xcoef}, includes an intercept (the first column in \emph{k} of \emph{xcoef}) 
and the coefficients for the explanatory variables.  The dimension of \emph{xcoef} is \emph{n x ntau x k}.  

Estimation can be very slow when \emph{target = "alldata"}.  The \emph{maketarget} command can be used to identify target points. 

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
Cleveland, William S. and Susan J. Devlin, "Locally Weighted Regression:  An Approach to Regression Analysis by Local Fitting," 
\emph{Journal of the American Statistical Association} 83 (1988), 596-610. 

Loader, Clive. \emph{ Local Regression and Likelihood.}  New York:  Springer, 1999. 

Koenker, Roger.  \emph{Quantile Regression.}  New York:  Cambridge University Press, 2005.
Chapter 7 and Appendix A.9.

McMillen, Daniel P., "Issues in Spatial Data Analysis," \emph{Journal of Regional Science} 50 (2010), 119-141.  

McMillen, Daniel P. and Christian Redfearn, ``Estimation and Hypothesis Testing for Nonparametric Hedonic House Price Functions,'' 
\emph{Journal of Regional Science} 50 (2010), 712-733.

Pagan, Adrian and Aman Ullah.  \emph{Nonparametric Econometrics.}  New York:  Cambridge University Press, 1999.

}

\seealso{
  \link{qreglwr}
}

\examples{
data(cookdata)
cookdata$obs <- seq(1,nrow(cookdata))
cookdata <- cookdata[!is.na(cookdata$FAR),]
par(ask=TRUE)

# 1. CPAR LWR estimates, y = a(DCBD) + b(dcbd)*DCBD + u
fit <- qregcpar(LNFAR~DCBD,nonpar=~DCBD, taumat=c(.10,.50,.90), 
  kern="bisq", window=.30, data=cookdata)
o <- order(cookdata$DCBD)
plot(cookdata$DCBD[o], fit$yhat[o,1],type="l", main="Log Floor Area Ratio",
  xlab="Distance from CBD",ylab="Log FAR")
lines(cookdata$DCBD[o], fit$yhat[o,2])
lines(cookdata$DCBD[o], fit$yhat[o,3])

\dontrun{
# 2. CPAR estimates, y = a(lat,long) + b(lat,long)xDCBD + u
fit <- qregcpar(LNFAR~DCBD, nonpar=~LATITUDE+LONGITUDE, taumat=c(.10,.90), 
  kern="bisq", window=.30, distance="LATLONG", data=cookdata)
plot(cookdata$DCBD, cookdata$LNFAR,main="Log Floor Area Ratio",
  xlab="Distance from CBD",ylab="Log FAR")
points(cookdata$DCBD, fit$yhat[,1], col="red")
plot(cookdata$DCBD, cookdata$LNFAR,main="Log Floor Area Ratio",
  xlab="Distance from CBD",ylab="Log FAR")
points(cookdata$DCBD, fit$yhat[,2], col="red")

library(RColorBrewer)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap$yhat10[cookdata$obs] <- fit$yhat[,1]
cmap$yhat90[cookdata$obs] <- fit$yhat[,2]
cmap$yhat1090 <- cmap$yhat90 - cmap$yhat10
brks <- seq(min(cmap$yhat1090,na.rm=TRUE),max(cmap$yhat1090,na.rm=TRUE),length=9)
spplot(cmap,"yhat1090",at=brks,col.regions=rev(brewer.pal(9,"RdBu")),
  main="Difference between .10 and.90 Quantiles")
}
}



\keyword{Quantile Regression}
\keyword{Nonparametric}
\keyword{Conditionally Parametric}
