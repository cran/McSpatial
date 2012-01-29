\name{cparprobit}
\alias{cparprobit}
\title{
Conditionally Parametric Probit
}
\description{
  Estimates a probit model by maximizing a locally weighted likelihood function -- the probit equivalent of cparlwr
}

\usage{ 
  cparprobit(form,nonpar,window=.25,bandwidth=0,kern="tcub",
    distance="Mahal",alldata=FALSE,data=NULL) 
}

\arguments{
  \item{form }{Model formula}
  \item{nonpar }{List of either one or two variables for \emph{z}.  
Formats: \emph{cparprobit(y~xlist, nonpar=~z1, ...)} or \emph{cparprobit(y~xlist, nonpar=~z1+z2, ...)}.  
Important:  note the "~" before the first \emph{z} variable.
}
  \item{window }{Window size.   Default:  0.25. }
  \item{bandwidth }{Bandwidth.   Default:  not used.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{cparprobit} looks for the first two letters to determine which variable is latitude and which is longitude, 
so the data set must be attached first or specified using the data option; options like data$latitude will not work.  Default:  Mahal. }
  \item{alldata }{If \emph{alldata=T}, each observation is used as a target value for \emph{z}.  
When \emph{alldata=F}, the function is estimated at a set of points chosen by the \emph{locfit} program using an adaptive decision tree approach, 
and the \emph{akima} program is used to interpolate to the full set of observations.  Specifying \emph{alldata=T} can lead to long estimation times.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
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

\deqn{ \sum_{i=1}^n w_i \{ y_i log(\Phi (X_i \beta)) + (1-y_i) log(1-\Phi (X_i \beta) ) \} }{\sum w_i { y_i log(\Phi (X_i \beta)) + (1-y_i) log(1-\Phi_i (X \beta) ) }  }

where y is the discrete dependent variable and X is the set of explanatory variables.
When \emph{Z} includes a single variable, \eqn{w_i} is a simple kernel weighting function:  \eqn{ w_i = K((z_i - z_0 )/(sd(z)*h)) }.  
When \emph{Z} includes two variables (e.g., nonpar=~z1+z2), the method for specifying \emph{w} depends on the \emph{distance} option. 
Under either option, the \emph{i}th row of the matrix \emph{Z} = (z1, z2) is transformed such  
that \eqn{z_i = \sqrt{z_i * V * t(z_i)}.}{z_i = sqrt(z_i * V * t(z_i)).} Under the "Mahal" option, \emph{V} is the inverse of cov(\emph{Z}).  
Under the \emph{"Euclid"} option, \emph{V} is the inverse of diag(cov(\emph{Z})).  
After this transformation, the weights again reduce to the simple kernel weighting function \eqn{K((z_i - z_0 )/(sd(z)*h))}.   
\emph{h} is specified by the \emph{bandwidth} or \emph{window} option.  

The great circle formula is used to constuct the distances used to form the weights when \emph{distance = "Latlong"}; 
in this case, the variable list for \emph{nonpar} must be listed as 
\emph{nonpar = ~latitude+longitude} (or \emph{~lo+la} or \emph{~lat+long}, etc), with the longitude and latitude variables expressed in degrees 
(e.g.,  -87.627800 and 41.881998 for one observation of longitude and latitude, respectively).   
The order in which latitude and longitude are listed does not matter and the function only looks for the 
first two letters to determine which variable is latitude and which is the longitude.  
It is important to note that the great circle distance measure is left in miles rather than being standardized.  
Thus, the window option should be specified when \emph{distance = "Latlong"} or the bandwidth should be adjusted to account for the scale.  
The kernel weighting function becomes \emph{K(distance/h)} under the \emph{"Latlong"} option. 

Following White (1982), the covariance matrix for a quasi-maximum likelihood model is \eqn{A^{-1}BA^{-1} }, where

\deqn{A = \sum_{i=1}^n w_i \frac{\partial^2 LnL_i}{\partial \beta \partial \beta ^\prime} }{A = \sum w_i d^2LnL_i/d\beta d\beta' }
\deqn{B = \sum_{i=1}^n w_i^2 \frac{\partial LnL_i}{\partial \beta}\frac{\partial LnL_i}{\partial \beta ^\prime} }{B = \sum w_i^2 (dLnL_i/d\beta)(dLnL_i/d\beta')  }
For the probit model, 
\deqn{ A = \sum_{i=1}^n w_i \frac{ \phi_i^2}{\Phi_i(1-\Phi_i)} X_i X_i ^\prime }{ A = \sum w_i \phi_i^2/(\Phi_i(1-\Phi_i)) X_i X_i' }
\deqn{ B = \sum_{i=1}^n w_i^2 \frac{ (y_i - X_i \beta)^2 \phi_i^2}{(\Phi_i(1-\Phi_i))^2} X_i X_i ^\prime }{ B = \sum w_i^2 (y_i - X_i \beta)^2 \phi_i^2/((\Phi_i(1-\Phi_i))^2) X_i X_i' }

The covariance matrix is calculated at each target point and the implied standard errors are then interpolated to each data point.

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

Fan, Jianqing, Nancy E. Heckman, and M.P. Wand, "Local Polynomial Kernel Regression for Generalized Linear Models and Quasi-Likelihood Functions,"
\emph{Journal of the American Statistical Association} 90 (1995), 141-150.

Loader, Clive. \emph{Local Regression and Likelihood.}  New York:  Springer, 1999. 

McMillen, Daniel P. and John F. McDonald, "Locally Weighted Maximum Likelihood Estimation:  Monte Carlo Evidence and an Application,"
in Luc Anselin, Raymond J.G.M. Florax, and Sergio J. Rey, eds., \emph{Advances in Spatial Econometrics}, Springer-Verlag, New York (2004), 225-239.

Tibshirani, Robert and Trevor Hastie, "Local Likelihood Estimation," \emph{Journal of the American Statistical Association} 82 (1987), 559-568.

}

\seealso{
  \link{cparlogit}

  \link{cparmlogit}

  \link{gmmlogit}

  \link{gmmprobit}

  \link{splogit}

  \link{spprobit}

  \link{spprobitml}

}

\examples{

library(spdep)
data(cookdata)
cookdata <- cookdata[cookdata$CHICAGO==1,]
cookdata$DCBD2 <- cookdata$DCBD^2

# Monte Carlo
set.seed(484909)
n = length(cookdata$DCBD)
o <- order(cookdata$DCBD)
cookdata$ybase <- 2*cookdata$DCBD -.15*cookdata$DCBD2 
cookdata$yvar <- cookdata$ybase + 2*rnorm(n)
fit <- lm(yvar~DCBD+DCBD2,data=cookdata)
summary(fit)
cookdata$yvar <- cookdata$yvar>5.6
table(cookdata$yvar)
fit1 <- glm(yvar~DCBD+DCBD2,family=binomial(link="probit"),data=cookdata)
summary(fit1)
fit2 <- glm(yvar~DCBD,family=binomial(link="probit"),data=cookdata)
summary(fit2)
nlfit1 <- cparprobit(yvar~DCBD,nonpar=~DCBD, window=.5, kern="rect", data=cookdata) 
nlfit1$lnl
nlfit2 <- cparprobit(yvar~DCBD+DCBD2,nonpar=~DCBD, window=.5, kern="rect", data=cookdata) 
nlfit2$lnl
plot(cookdata$DCBD[o],fitted(fit1)[o],xlab="Distance from CBD",ylab="p",type="l")
lines(cookdata$DCBD[o],fitted(fit2)[o],col="red")
lines(cookdata$DCBD[o],nlfit1$p[o],col="blue")
lines(cookdata$DCBD[o],nlfit2$p[o],col="green")
legend("topright",c("True Model","Linear","CPAR-linear","CPAR-True"),
  col=c("black","red","blue","green"),lwd=1)


}

\keyword{Conditionally Parametric}
\keyword{Nonparametric}
\keyword{Discrete Choice Models}

