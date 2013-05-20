\name{cparmlogit}
\alias{cparmlogit}
\title{
Conditionally parametric logit for two or more choices
}
\description{
  Estimates a multinomial logit model with two or more choices by maximizing a locally weighted likelihood function -- the logit equivalent of cparlwr
}

\usage{ 
cparmlogit(form,nonpar,window=.25,bandwidth=0,kern="tcub",
distance="Mahal",target=NULL,data=NULL)  
}

\arguments{
  \item{form }{Model formula}
  \item{nonpar }{List of either one or two variables for \emph{z}.  
Formats: \emph{cparmlogit(y~xlist, nonpar=~z1, ...)} or \emph{cparmlogit(y~xlist, nonpar=~z1+z2, ...)}.  
Important:  note the "~" before the first \emph{z} variable.
}
  \item{window }{Window size.   Default:  0.25. }
  \item{bandwidth }{Bandwidth.   Default:  not used.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{cparmlogit} looks for the first two letters to determine which variable is latitude and which is longitude, 
so the data set must be attached first or specified using the data option; options like data$latitude will not work.  Default:  Mahal. }
  \item{target}{If \emph{target = NULL}, uses the \emph{maketarget} command to form targets using the values specified for \emph{window}, 
\emph{bandwidth}, and \emph{kern}.  If \emph{target="alldata"}, each observation is used as a target value for \emph{x}.  
A set of target values can be supplied directly.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}

\value{
  \item{target}{The target points for the original estimation of the function.}
  \item{xcoef.target}{Estimated coefficients, \emph{B(z)}, at the target values of \emph{z}.}
  \item{xcoef.target.se}{Standard errors for \emph{B(z)} at the target values of \emph{z}.}
  \item{xcoef}{Estimated coefficients, \emph{B(z)}, at the original data points.}
  \item{xcoef.se}{Standard errors for \emph{B(z)} with \emph{z} evaluated at all points in the data set.}
  \item{pmat}{The n x K+1 matrix of estimated probabilities.}
  \item{lnl}{The log-likelihood value.}
}

\details{
The list of explanatory variables is specified in the base model formula while \emph{Z} is specified using \emph{nonpar}.  
\emph{X} can include any number of explanatory variables, but \emph{Z} must have at most two. 

The model is estimated by maximizing the following weighted log-likelihood function at each target point:

\deqn{ \sum_{i=1}^n \sum_{j=1}^K w_i I(y_i=j) log(P(X_i \beta_j)) }{ \sum_i \sum_j w_i I(y_i=j) log(P(X_i \beta_j)) }

where y is the discrete dependent variable with K+1 choices, X is the set of explanatory variables, and
\eqn{P(X_i \beta_j) = \frac{exp(X_i \beta_j)}{\sum_j exp(X_i \beta_j)}.}{P(X_i \beta_j) = exp(X_i \beta_j) / \sum_j exp(X_i \beta_j).} 
For the base value, y=0, the coefficients are normalized to \eqn{\beta_0 = 0}.

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

The covariance matrix is calculated at each target point and the implied standard errors are then interpolated to each data point.
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

Fan, Jianqing, Nancy E. Heckman, and M.P. Wand, "Local Polynomial Kernel Regression for Generalized Linear Models and Quasi-Likelihood Functions,"
\emph{Journal of the American Statistical Association} 90 (1995), 141-150.

Loader, Clive. \emph{Local Regression and Likelihood.}  New York:  Springer, 1999. 

McMillen, Daniel P. and John F. McDonald, "Locally Weighted Maximum Likelihood Estimation:  Monte Carlo Evidence and an Application,"
in Luc Anselin, Raymond J.G.M. Florax, and Sergio J. Rey, eds., \emph{Advances in Spatial Econometrics}, Springer-Verlag, New York (2004), 225-239.

Tibshirani, Robert and Trevor Hastie, "Local Likelihood Estimation," \emph{Journal of the American Statistical Association} 82 (1987), 559-568.

}


\seealso{
  \link{cparlogit}

  \link{cparprobit}

  \link{gmmlogit}

  \link{gmmprobit}

  \link{splogit}

  \link{spprobit}

  \link{spprobitml}

}

\examples{
library(mlogit)
set.seed(5647)
n = 1000
x <- runif(n,0,pi*sqrt(12))
o <- order(x)
x <- x[o]
form <- yvar~x
nonpar <- ~x

# 2 choices
ybase <- x + rlogis(n)
yvar <- ybase>.5*pi*sqrt(12)
table(yvar)
fit <- glm(yvar~x,family=binomial(link="logit"))
summary(fit)
p <- fitted(fit)
fit1 <- cparmlogit(yvar~x,nonpar=~x,window=.5,kern="tcub")
fit1$lnl
colMeans(fit1$xcoef)
colMeans(fit1$xcoef.se)
cor(p,fit1$pmat)
plot(x,p,xlab="x",ylab="Prob(y=1)",type="l")
lines(x,fit1$pmat[,2],col="red")
legend("topleft",c("Standard Logit","CPAR"),col=c("black","red"),lwd=1)

\dontrun{
par(ask=TRUE)
# 3 choices
ybase1 <- -.5*pi*sqrt(12) + x + rlogis(n)
ybase2 <-  -.5*pi*sqrt(12)/2 + x/2 + rlogis(n)
yvar <- ifelse(ybase1>ybase2,1,2)
yvar <- ifelse(ybase1<0&ybase2<0,0,yvar)
table(yvar)
mdata <- data.frame(yvar,x)
fit <- mlogit(yvar~0 | x, data=mdata, shape="wide")
summary(fit)
fit1 <- cparmlogit(yvar~x,nonpar=~x,window=.5,kern="tcub")
fit1$lnl
colMeans(fit1$xcoef)
colMeans(fit1$xcoef.se)
cor(fit$probabilities,fit1$pmat)
plot(x,fit$probabilities[,1],xlab="x",ylab="Prob(y=1)",type="l",main="Prob(y=0)")
lines(x,fit1$pmat[,1],col="red")
legend("topright",c("Standard Logit","CPAR"),col=c("black","red"),lwd=1)
plot(x,fit$probabilities[,2],xlab="x",ylab="Prob(y=1)",type="l",main="Prob(y=1)")
lines(x,fit1$pmat[,2],col="red")
legend("topleft",c("Standard Logit","CPAR"),col=c("black","red"),lwd=1)
plot(x,fit$probabilities[,3],xlab="x",ylab="Prob(y=1)",type="l",main="Prob(y=2)")
lines(x,fit1$pmat[,3],col="red")
legend("topleft",c("Standard Logit","CPAR"),col=c("black","red"),lwd=1)

# 2 choices, quadratic
x2 <- x^2
ybase <- x - .1*(x^2) + rlogis(n)
yvar <- ybase>median(ybase)
table(yvar)
fit <- glm(yvar~x+x2,family=binomial(link="logit"))
summary(fit)
p <- fitted(fit)
fit1 <- cparmlogit(yvar~x,nonpar=~x,window=.25,kern="tcub")
fit1$lnl
colMeans(fit1$xcoef)
colMeans(fit1$xcoef.se)
cor(p,fit1$pmat)
plot(x,p,xlab="x",ylab="Prob(y=1)",type="l")
lines(x,fit1$pmat[,2],col="red")
legend("topleft",c("Standard Logit","CPAR"),col=c("black","red"),lwd=1)
}
}

\keyword{Conditionally Parametric}
\keyword{Nonparametric}
\keyword{Discrete Choice Models}

