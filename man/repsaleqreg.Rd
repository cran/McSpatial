\name{repsaleqreg}
\alias{repsaleqreg}
\title{
  Quantile Repeat Sales Estimation
}
\description{
  Median-Based Repeat Sales Estimation
}

\usage{
repsaleqreg(price0,time0,price1,time1,mergefirst=1,
  graph=TRUE,graph.conf=TRUE,conf=.95,print=TRUE) 
}

\arguments{
  \item{price0}{Earlier price in repeat sales pair }
  \item{time0}{Earlier time in repeat sales pair }
  \item{price1}{Later price in repeat sales pair }
  \item{time1}{Later time in repeat sales pair }
  \item{mergefirst}{Number of initial periods with coefficients constrained to zero. Default:  \emph{mergefirst=1}}
  \item{graph}{If TRUE, graph results.  Default:  \emph{graph=T}}
  \item{graph.conf}{If TRUE, add confidence intervals to graph.  Default:  \emph{graph.conf=T} }
  \item{conf}{Confidence level for intervals.  Default:  .95}
  \item{print}{If \emph{print=T}, prints the regression results. Default:  \emph{print=T}.}
}

\value{
  \item{fit}{Full quantile regression model.}
  \item{pindex}{The estimated price index.}
  \item{lo}{The lower bounds for the price index confidence intervals.}
  \item{hi}{The upper bounds for the price index confidence intervals.}
}


\examples{
set.seed(189)
n = 2000
# sale dates range from 0-10
# drawn uniformly from all possible time0, time1 combinations with time0<time1
tmat <- expand.grid(seq(0,10), seq(0,10))
tmat <- tmat[tmat[,1]<tmat[,2], ]
tobs <- sample(seq(1:nrow(tmat)),n,replace=TRUE)
time0 <- tmat[tobs,1]
time1 <- tmat[tobs,2]
timesale <- time1-time0
table(timesale)

# constant variance; index ranges from 0 at time 0 to 1 at time 10
y0 <- time0/10 + rnorm(n,0,.2)
y1 <- time1/10 + rnorm(n,0,.2)
fit <- repsaleqreg(price0=y0, price1=y1, time0=time0, time1=time1)

# variance rises with timesale
# var(u0) = .2^2; var(u1) = (.2 + timesale/10)^2
# var(u1-u0) = var(u0) + var(u1) = 2*(.2^2) + .4*timesale/10 + (timesale^2)/100
y0 <- time0/10 + rnorm(n,0,.2)
y1 <- time1/10 + rnorm(n,0,.2+timesale/10)
par(ask=TRUE)
fit <- repsaleqreg(price0=y0, price1=y1, time0=time0, time1=time1)
summary(fit$pindex)

}


\details{
The repeat sales model is
\deqn{ y_t - y_s = \delta_t - \delta_s + u_t - u_s}{ y(t) - y(s) = \delta(t) - \delta(s) + u(t) - u(s)}

where \emph{y} is the log of sales price, \emph{s} denotes the earlier sale in a repeat sales pair, and \emph{t} denotes the later sale.  
Each entry of the data set should represent a repeat sales pair, with \eqn{price0 = y_s}{price0 = y(s)}, \eqn{price1 = y_t}{price1 = y(t)}, 
\eqn{time0 = s}, and \eqn{time1 = t}.  The function \emph{repsaledata} can help transfer a standard hedonic data set to a set of 
repeat sales pairs.

Repeat sales estimates are sometimes very sensitive to sales from the first few time periods, particularly when the sample size is small.
The option \emph{mergefirst} indicates the number of time periods for which the price index is constrained to equal zero.  The default is 
\emph{mergefirst = 1}, meaning that the price index equals zero for just the first time period.  

The \emph{repsaleqreg} function uses the \emph{quantreg} package to estimate a quantile regression for the .50 quantile, i.e., the median.  
A median-based estimator is less sensitive to outliers than linear regression.  McMillen and Thorsnes (2006) show that the quantile approach is less 
sensitive to the inclusion of properties that have undergone renovations between sales.  \emph{repsaleqreg} first fits a standard quantile model, 
including the interecpt.  The coefficient vector is then rotated to have a zero intercept using the formula for transforming unrestricted linear 
regression estimates to the restricted (zero intercept) values: 

  fit <- rq(dy~x) \cr
  b <- fit$coef \cr
  fit1 <- summary(fit,covariance=TRUE)  \cr 
  vmat <- fit1$cov \cr
  k = length(b1) \cr
  rmat <- diag(k) \cr
  rmat[,1] <- rmat[,1] - vmat[1,]/vmat[1,1] \cr
  bmat <- rmat%*%b1 
}

\references{
Case, Karl and Robert Shiller, "Prices of Single-Family Homes since 1970:  New Indexes for Four Cities," \emph{New England Economic Review} (1987), 45-56.

McMillen, Daniel P. and Paul Thorsnes, "Housing Renovations and the Quantile Repeat Sales Price Index," \emph{Real Estate Economics} 34 (2006), 567-587.

}

\seealso{
  \link{repsale}

  \link{repsaledata}

  \link{repsalefourier}
}

\keyword{Repeat Sales}
\keyword{Quantile Regression}
