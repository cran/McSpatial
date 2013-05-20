\name{cparlwr}
\alias{cparlwr}
\title{
Conditionally Parametric LWR Estimation
}
\description{
Estimates a model of the form \emph{y = XB(z) + u}, where \emph{z} can include one or two variables.  
"Geographically weighted regression" is a special case in which \emph{z = (latitude, longitude)} or some other measure of location.}

\usage{ 
cparlwr(form,nonpar,window=.25,bandwidth=0,kern="tcub",
  distance="Mahal",targetobs=NULL,data=NULL)
}

\arguments{
  \item{form }{Model formula}
  \item{nonpar }{List of either one or two variables for \emph{z}.  
Formats: \emph{cparlwr(y~xlist, nonpar=~z1, ...)} or \emph{cparlwr(y~xlist, nonpar=~z1+z2, ...)}.  
Important:  note the "~" before the first \emph{z} variable.
}
  \item{window }{Window size.   Default:  0.25. }
  \item{bandwidth }{Bandwidth.   Default:  not used.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{cparlwr} looks for the first two letters to determine which variable is latitude and which is longitude, 
so the data set must be attached first or specified using the data option; options like data$latitude will not work.  Default:  Mahal. }
  \item{targetobs}{If \emph{targetobs = NULL}, uses the \emph{maketarget} command to form targets.  
If \emph{target="alldata"}, each observation is used as a target value for \emph{x}.  
A set of target can also be supplied directly by listing the observation numbers of the target data points.  The observation numbers can be
identified using the \emph{obs} variable produced by the \emph{maketarget} command.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}

\value{
  \item{target}{The target points for the original estimation of the function.}
  \item{ytarget}{The predicted values of \emph{y} at the target values \emph{z}.}
  \item{xcoef.target}{Estimated coefficients, \emph{B(z)}, at the target values of \emph{z}.}
  \item{xcoef.target.se}{Standard errors for \emph{B(z)} at the target values of \emph{z}.}
  \item{yhat}{Predicted values of \emph{y} at the original data points.}
  \item{xcoef}{Estimated coefficients, \emph{B(z)}, at the original data points.}
  \item{xcoef.se}{Standard errors for \emph{B(z)} with \emph{z} evaluated at all points in the data set.}
  \item{df1}{\emph{tr(L)}, a measure of the degrees of freedom used in estimation.}
  \item{df2}{\emph{tr(L'L)}, an alternative measure of the degrees of freedom used in estimation.  }
  \item{sig2}{Estimated residual variance, \emph{sig2 = rss/(n-2*df1+df2)}.}
  \item{cv}{Cross-validation measure.  \emph{cv = mean(((y-yhat)/(1-infl))^2) }, 
where \emph{yhat} is the vector of predicted values for \emph{y} and \emph{infl} is the vector of diagonal terms for \emph{L}. }
  \item{gcv}{gcv = n*(n*sig2)/((n-nreg)^2), where \emph{sig2} is the estimated residual variance and \emph{nreg = 2*df1 - df2}.}
  \item{infl}{A vector containing the diagonal elements of \emph{L}.}
}

\details{
The list of explanatory variables is specified in the base model formula while \emph{Z} is specified using \emph{nonpar}.  
The model formula does not have to include an intercept, making it suitable for repeat sales estimation as well as other models.
\emph{X} can include any number of explanatory variables, but \emph{Z} must have at most two.  
\emph{cparlwr} is equivalent to the \emph{lwr} command when \emph{Z = X} and the formula includes an intercept, with one exception:  
the explanatory variables are not centered on the target points so the intercept does not provide a direct estimate of y. 
This affects the intercept and its standard errors but not the coefficients on the explanatory variables.  It also means that
\eqn{\hat{y} = \hat{\alpha} + X \hat{\beta} } rather than just \eqn{\hat{\alpha}. }  The estimated coefficient matrix, 
\emph{xcoef}, provides estimates of the slopes at \eqn{z_0} , i.e., \eqn{B(z_0)}

The estimated value of \emph{y} at a target value \eqn{z_0} is the predicted value from a weighted least squares regression 
of \emph{y} on \emph{X} with weights given by \emph{K}.
When \emph{Z} includes a single variable, \emph{K} is a simple kernel weighting function:  \eqn{K((z - z_0 )/(sd(z)*h)) }.  
When \emph{Z} includes two variables (e.g., nonpar=~z1+z2), the method for specifying \emph{K} depends on the \emph{distance} option. 
Under either option, the \emph{i}th row of the matrix \emph{Z} = (z1, z2) is transformed such  
that \eqn{z_i = \sqrt{z_i * V * t(z_i)}.}{z_i = sqrt(z_i * V * t(z_i)).} Under the "Mahal" option, \emph{V} is the inverse of cov(\emph{Z}).  
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
The function \emph{cparlwrgrid} can be used to search for the value of \emph{h} that minimizes the \emph{cv} or \emph{gcv} criterion.  

Since each estimate is a linear function of all \emph{n} values for \emph{y}, the full set of estimates takes the form \emph{yhat = LY}, where \emph{L} is an \emph{nxn} matrix. 
Loader (1999) suggests two measures of the number of degrees of freedom used in estimation:  \emph{df1 = tr(L)} and \emph{df2 = tr(L'L)}.  
The diagonal elements of \emph{tr(L)} are stored in the array \emph{infl}.  
Since the degrees of freedom measures can differ substantially when \emph{target="alldata"} rather than using a set of target points, 
it is a good idea to report final estimates using \emph{target="alldata"} when possible.

Again following Loader (1999), the degrees of freedom correction used to estimate the error variance, \emph{sig2}, is \emph{df = 2*df1 -  df2}.  
Let \emph{e} represent the vector of residuals, \emph{y - yhat}.  The estimated variance is \eqn{sig2 = \sum_i e_i^2 / (n-df)}{sig2 = \sum e^2/(n-df)}. 
The covariance matrix for \eqn{B(z_0)} 
is \deqn{\hat{\sigma}^2(\sum_{i=1}^n X_i K(\phi_i) X_i^\top)^{-1}(\sum_{i=1}^n X_i (K(\phi_i))^2 X_i^\top )(\sum_{i=1}^n X_i K(\phi_i) X_i^\top)^{-1}.}{sig2 (\sum XK(\phi)X')^(-1)(\sum XK(\phi)^2X')(\sum XK(\phi)X')^(-1) }

Estimation can be very slow when \emph{targetobs = "alldata"}.  The \emph{maketarget} command can be used to identify target points. 

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

McMillen, Daniel P., "One Hundred Fifty Years of Land Values in Chicago:  A Nonparametric Approach," \emph{Journal of Urban Economics} 40 (1996), 100-124. 

McMillen, Daniel P., "Issues in Spatial Data Analysis," \emph{Journal of Regional Science} 50 (2010), 119-141.  

McMillen, Daniel P., "Employment Densities, Spatial Autocorrelation, and Subcenters in Large Metropolitan Areas," \emph{Journal of Regional Science} 44 (2004), 225-243. 

McMillen, Daniel P. and John F. McDonald, "A Nonparametric Analysis of Employment Density in a Polycentric City," \emph{Journal of Regional Science} 37 (1997), 591-612. 

McMillen, Daniel P. and Christian Redfearn, ``Estimation and Hypothesis Testing for Nonparametric Hedonic House Price Functions,'' 
\emph{Journal of Regional Science} 50 (2010), 712-733.

Pagan, Adrian and Aman Ullah.  \emph{Nonparametric Econometrics.}  New York:  Cambridge University Press, 1999.

}

\examples{

data(cookdata)
par(ask=TRUE)
cookdata <- cookdata[cookdata$CHICAGO==1&!is.na(cookdata$LNFAR),]
fit1 <- cparlwr(LNFAR~DCBD,nonpar=~DCBD, window=.10, 
  data=cookdata)
fit2 <- cparlwr(LNFAR~DCBD,nonpar=~LONGITUDE+LATITUDE,window=.10,
  distance="LATLONG",data=cookdata)
cookdata$yhat1 <- fit1$yhat
cookdata$yhat2 <- fit2$yhat
o <- order(cookdata$DCBD)
plot(cookdata$DCBD[o], cookdata$LNFAR[o],main="Log Floor Area Ratio",
  xlab="Distance from CBD",ylab="Log FAR")
lines(cookdata$DCBD[o], cookdata$yhat1[o], col="red")
plot(cookdata$DCBD[o], cookdata$LNFAR[o],main="Log Floor Area Ratio",
  xlab="Distance from CBD",ylab="Log FAR")
points(cookdata$DCBD[o], cookdata$yhat2[o], col="red")

}


\seealso{
  \link{cparlwrgrid}

  \link{cubespline}

  \link{fourier}

  \link{lwr}

  \link{lwrgrid}

  \link{semip}
}

\keyword{Conditionally Parametric}
\keyword{Nonparametric}
\keyword{Locally Weighted Regression}


