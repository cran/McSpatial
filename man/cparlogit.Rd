\name{cparlogit}
\alias{cparlogit}
\title{
Conditionally Parametric logit for two or more choices
}
\description{
  Estimates a logit model with two choices by maximizing a locally weighted likelihood function -- the logit equivalent of cparlwr
}

\usage{ 
cparlogit(form,nonpar,window=.25,bandwidth=0,kern="tcub",
distance="Mahal",target=NULL,data=NULL,minp=NULL)  
}

\arguments{
  \item{form }{Model formula}
  \item{nonpar }{List of either one or two variables for \emph{z}.  
Formats: \emph{cparlogit(y~xlist, nonpar=~z1, ...)} or \emph{cparlogit(y~xlist, nonpar=~z1+z2, ...)}.  
Important:  note the "~" before the first \emph{z} variable.
}
  \item{window }{Window size.   Default:  0.25. }
  \item{bandwidth }{Bandwidth.   Default:  not used.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{cparlogit} looks for the first two letters to determine which variable is latitude and which is longitude, 
so the data set must be attached first or specified using the data option; options like data$latitude will not work.  Default:  Mahal. }

  \item{target}{If \emph{target = NULL}, uses the \emph{maketarget} command to form targets using the values specified for \emph{window}, 
\emph{bandwidth}, and \emph{kern}.  If \emph{target="alldata"}, each observation is used as a target value for \emph{x}.  
A set of target values can be supplied directly.}

  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
  \item{minp}{Specifies a limit for the estimated probability.  Any estimated probability lower than \emph{minp} will be set to \emph{minp} and
any probability higher than 1-\emph{minp} will be set to 1-\emph{minp}.  By default, the estimated probabilities are bounded by 0 and 1.}
}

\value{
  \item{target}{The target points for the original estimation of the function.}
  \item{xcoef.target}{Estimated coefficients, \emph{B(z)}, at the target values of \emph{z}.}
  \item{xcoef.target.se}{Standard errors for \emph{B(z)} at the target values of \emph{z}.}
  \item{xcoef}{Estimated coefficients, \emph{B(z)}, at the original data points.}
  \item{xcoef.se}{Standard errors for \emph{B(z)} with \emph{z} evaluated at all points in the data set.}
  \item{p}{The estimated probabilities.}
  \item{lnl}{The log-likelihood value.}
}

\details{
The list of explanatory variables is specified in the base model formula while \emph{Z} is specified using \emph{nonpar}.  
\emph{X} can include any number of explanatory variables, but \emph{Z} must have at most two. 

The model is estimated by maximizing the following weighted log-likelihood function at each target point:

\deqn{ \sum_{i=1}^n w_i \{ y_i log(P_i) + (1-y_i) log(1-P_i ) \} }{\sum w_i { y_i log(P_i) + (1-y_i) log(1-P_i ) }  }


where y is the discrete dependent variable, X is the set of explanatory variables, and
\eqn{P_i = \frac{exp(X_i \beta)}{1 + exp(X_i \beta)}.}{P_i = exp(X_i \beta) / (1+ exp(X_i \beta)).} 

When \emph{Z} includes a single variable, \eqn{w_i} is a simple kernel weighting function:  \eqn{ w_i = K((z_i - z_0 )/(sd(z)*h)) }.  
When \emph{Z} includes two variables (e.g., nonpar=~z1+z2), the method for specifying \emph{w} depends on the \emph{distance} option. 
Under either option, the \emph{i}th row of the matrix \emph{Z} = (z1, z2) is transformed such  
that \eqn{z_i = \sqrt{z_i * V * t(z_i)}.}{z_i = sqrt(z_i * V * t(z_i)).} Under the "Mahal" option, \emph{V} is the inverse of cov(\emph{Z}).  
Under the \emph{"Euclid"} option, \emph{V} is the inverse of diag(cov(\emph{Z})).  
After this transformation, the weights again reduce to the simple kernel weighting function \eqn{K((z_i - z_0 )/(sd(z)*h))}.   
\emph{h} is specified by the \emph{bandwidth} or \emph{window} option.  

The great circle formula is used to construct the distances used to form the weights when \emph{distance = "Latlong"}; 
in this case, the variable list for \emph{nonpar} must be listed as 
\emph{nonpar = ~latitude+longitude} (or \emph{~lo+la} or \emph{~lat+long}, etc), with the longitude and latitude variables expressed in degrees 
(e.g.,  -87.627800 and 41.881998 for one observation of longitude and latitude, respectively).   
The order in which latitude and longitude are listed does not matter and the function only looks for the 
first two letters to determine which variable is latitude and which is longitude.  
It is important to note that the great circle distance measure is left in miles rather than being standardized.  
Thus, the window option should be specified when \emph{distance = "Latlong"} or the bandwidth should be adjusted to account for the scale.  
The kernel weighting function becomes \emph{K(distance/h)} under the \emph{"Latlong"} option. 

Following White (1982), the covariance matrix for a quasi-maximum likelihood model is \eqn{A^{-1}BA^{-1} }, where

\deqn{A = \sum_{i=1}^n w_i \frac{\partial^2 LnL_i}{\partial \beta \partial \beta ^\prime} }{A = \sum w_i d^2LnL_i/d\beta d\beta' }
\deqn{B = \sum_{i=1}^n w_i^2 \frac{\partial LnL_i}{\partial \beta}\frac{\partial LnL_i}{\partial \beta ^\prime} }{B = \sum w_i^2 (dLnL_i/d\beta)(dLnL_i/d\beta')  }
For the logit model, 
\deqn{ A = \sum_{i=1}^n w_i P_i(1 - P_i) X_i X_i ^\prime }{ A = \sum w_i P_i(1 - P_i) X_i X_i' }
\deqn{ B = \sum_{i=1}^n w_i^2 (y_i - P_i)^2 X_i X_i ^\prime }{ B = \sum w_i^2 (y_i - P_i)^2  X_i X_i' }

The covariance matrix is calculated at all target points and the implied standard errors are then interpolated to each data point.

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

Fan, Jianqing, Nancy E. Heckman, and M.P. Wand, "Local Polynomial Kernel Regression for Generalized Linear Models and Quasi-Likelihood Functions,"
\emph{Journal of the American Statistical Association} 90 (1995), 141-150.

Loader, Clive. \emph{Local Regression and Likelihood.}  New York:  Springer, 1999. 

McMillen, Daniel P. and John F. McDonald, "Locally Weighted Maximum Likelihood Estimation:  Monte Carlo Evidence and an Application,"
in Luc Anselin, Raymond J.G.M. Florax, and Sergio J. Rey, eds., \emph{Advances in Spatial Econometrics}, Springer-Verlag, New York (2004), 225-239.

Tibshirani, Robert and Trevor Hastie, "Local Likelihood Estimation," \emph{Journal of the American Statistical Association} 82 (1987), 559-568.

}


\seealso{

  \link{cparprobit}

  \link{cparmlogit}

  \link{gmmlogit}

  \link{gmmprobit}

  \link{splogit}

  \link{spprobit}

  \link{spprobitml}

}

\examples{
set.seed(5647)
data(cookdata)
cookdata <- cookdata[!is.na(cookdata$AGE),]
n = nrow(cookdata)
cookdata$ystar <- cookdata$DCBD - .5*cookdata$AGE
cookdata$y <- cookdata$ystar - mean(cookdata$ystar) + rnorm(n,sd=4) > 0

fit  <- cparlogit(y~DCBD+AGE,~LONGITUDE+LATITUDE,window=.5,
  distance="Latlong",data=cookdata,minp=0.001)
}


\keyword{Discrete Choice Models}
\keyword{Logit}
\keyword{Conditionally Parametric}
\keyword{Nonparametric}

