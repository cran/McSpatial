\name{repsalefourier}
\alias{repsalefourier}
\title{
  Repeat Sales Estimation using Fourier Expansions
}
\description{
  Standard and Weighted Least Squares Repeat Sales Estimation using Fourier Expansions
}

\usage{
repsalefourier(price0,time0,price1,time1,mergefirst=1,q=1, graph=TRUE,
  graph.conf=TRUE,conf=.95,stage3=FALSE,stage3_xlist=~timesale,
  print=TRUE) 
}

\arguments{
  \item{price0}{Earlier price in repeat sales pair }
  \item{time0}{Earlier time in repeat sales pair }
  \item{price1}{Later price in repeat sales pair }
  \item{time1}{Later time in repeat sales pair }
  \item{mergefirst}{Number of initial periods with coefficients constrained to zero. Default:  \emph{mergefirst=1}}
  \item{q}{Sets \emph{Q} for the fourier expansion.  Default: \emph{q}=1.}
  \item{graph}{If TRUE, graph results.  Default:  \emph{graph=T}}
  \item{graph.conf}{If TRUE, add confidence intervals to graph.  Default:  \emph{graph.conf=T} }
  \item{conf}{Confidence level for intervals.  Default:  .95}
  \item{stage3}{If \emph{stage3 = NULL}, no corrections for heteroskedasticty.  If \emph{stage3="abs"}, uses the absolute value of the 
    first-stage residuals as the dependent variable in the second-stage regression.  If \emph{stage3="square"}, uses the square of the
    first-stage residuals as the dependent variable.  Default:  \emph{stage3=NULL}.}
  \item{stage3_xlist}{List of explanatory variables for heteroskedasticity.  By default, the single variable \emph{timesale = time1-time0} is
  constructed and used as the explanatory variable when \emph{stage3="abs"} or \emph{stage3="square"}.  Alternatively, a formula can be provided for
  a user-specified list of explanatory variables, e.g., \emph{stage3_xlist=~x1+x2}.  \emph{Important}:  note the "~" before the variable list.}
  \item{print}{If \emph{print=T}, prints the regression results.  Prints one stage only -- the first stage when \emph{stage=NULL} and 
the final stage when \emph{stage3="square"} or \emph{stage3="abs"}.  Default:  \emph{print=T}.}
}

\value{
  \item{fit}{Full regression model.}
  \item{pindex}{The estimated price index.}
  \item{lo}{The lower bounds for the price index confidence intervals.}
  \item{hi}{The upper bounds for the price index confidence intervals.}
  \item{dy}{The dependent variable for the repeat sales regression, \emph{dy = price1-price0.}}
  \item{xmat}{The matrix of explanatory variables for the repeat sales regressions.  \eqn{dim(xmat) = 2 + 2Q}. }
}

\examples{
set.seed(189)
n = 2000
# sale dates range from 0-50
# drawn uniformly from all possible time0, time1 combinations with time0<time1
tmat <- expand.grid(seq(0,50), seq(0,50))
tmat <- tmat[tmat[,1]<tmat[,2], ]
tobs <- sample(seq(1:nrow(tmat)),n,replace=TRUE)
time0 <- tmat[tobs,1]
time1 <- tmat[tobs,2]
timesale <- time1-time0
timesale2 <- timesale^2

par(ask=TRUE)
z0 <- 2*pi*time0/50
z0sq <- z0^2
sin0 <- sin(z0)
cos0 <- cos(z0)
z1 <- 2*pi*time1/50
z1sq <- z1^2
sin1 <- sin(z1)
cos1 <- cos(z1)
ybase0 <- z0 + .05*z0sq -.5*sin0 - .5*cos0
miny <- min(ybase0)
ybase0 <- ybase0-miny
ybase1 <- z1 + .05*z1sq -.5*sin1 - .5*cos1 - miny
maxy <- max(ybase1)
ybase0 <- ybase0/maxy
ybase1 <- ybase1/maxy
summary(data.frame(ybase0,ybase1))
sig1 = sd(c(ybase0,ybase1))/2
y0 <- ybase0 + rnorm(n,0,sig1)
y1 <- ybase1 + rnorm(n,0,sig1)
fit <- lm(y0~z0+z0sq+sin0+cos0)
summary(fit)
plot(time0,fitted(fit))
fit <- lm(y1~z1+z1sq+sin1+cos1)
summary(fit)
plot(time1,fitted(fit))

fit1 <- repsale(price1=y1,price0=y0,time1=time1,time0=time0,graph=FALSE,
  mergefirst=5)
fit2 <- repsalefourier(price1=y1,price0=y0,time1=time1,time0=time0,q=1,
  graph=FALSE,mergefirst=5)
timevar <- seq(0,50)
plot(timevar,fit1$pindex,type="l",xlab="Time",ylab="Index",
  ylim=c(min(fit1$pindex),max(fit2$pindex)))
lines(timevar,fit2$pindex)



# variance rises with timesale
# var(u0) = sig1^2; var(u1) = (sig1 + timesale/50)^2
# var(u1-u0) = var(u0) + var(u1) = 2*(sig1^2) + 2*sig1*timesale/10 + (timesale^2)/2500
y0 <- ybase0 + rnorm(n,0,sig1)
y1 <- ybase1 + rnorm(n,0,sig1+timesale/50)
par(ask=TRUE)
fit1 <- repsalefourier(price0=y0, price1=y1, time0=time0, time1=time1,
  graph=FALSE)
fit2 <- repsalefourier(price0=y0, price1=y1, time0=time0, time1=time1,
  graph=FALSE,stage3="abs",stage3_xlist=~timesale+timesale2)
plot(timevar,fit1$lo,type="l",xlab="Time",ylab="Index",
  ylim=c(min(fit1$lo,fit2$lo),max(fit1$hi,fit2$hi)))
lines(timevar,fit1$hi)
lines(timevar,fit2$lo,col="red")
lines(timevar,fit2$hi,col="red")

}


\details{
The repeat sales model is
\deqn{ y_t - y_s = \delta_t - \delta_s + u_t - u_s}{ y(t) - y(s) = \delta(t) - \delta(s) + u(t) - u(s)}

where \emph{y} is the log of sale price, \emph{s} denotes the earlier sale in a repeat sales pair, and \emph{t} denotes the later sale.  
Each entry of the data set should represent a repeat sales pair, with \eqn{price0 = y_s}{price0 = y(s)}, \eqn{price1 = y_t}{price1 = y(t)}, 
\eqn{time0 = s}, and \eqn{time1 = t}.  The function \emph{repsaledata} can help transfer a standard hedonic data set to a set of 
repeat sales pairs.

The repeat sales model can be derived from a hedonic price function with the form 
\eqn{y_{i,t} = \delta_t + X_i \beta + u_{i,t}} where \eqn{X_i} is a vector of variables that are assumed constant over time.
\emph{repsalefourier} replaces \eqn{\delta_t} with a smooth continuous function, \eqn{g(T_i)} where \eqn{T_i} denotes the time of sale for observation \emph{i}.
Letting \eqn{g(T_i) = \alpha_0 + \alpha_1 z_i + \alpha_2 z_i^2 + \sum_{i=1}^Q \{\lambda_q sin(qz_i) + \gamma_q cos(qz_i) \} },
where \eqn{z_i = 2 \pi (T_i - min(T_i))/(max(T_i) - min(T_i)) }, the repeat sales model becomes \eqn{y_{i,t} - y_{i,s} = g(T_i) - g(T_i^s)} =
\deqn{\alpha_1 (z_i - z_i^s) + \alpha_2 (z_i^2 - z_i^{s2}) + \sum_{q=1}^Q \{ \lambda_q (sin(qz_i) - sin(qz_i^s)) + \gamma_q (cos(qz_i) - cos(z_i^s)) \} + u_{i,t} - u_{i,t-s} }
After imposing the constraint that the price index in the base time period equals zero, the index is constructed from the estimated regression using the following expression:
\deqn{ g(T_i) = \alpha_1 z_i + \alpha_2 z_i^2 + \sum_{q=1}^Q \{ \lambda_q  sin(qz_i) + \gamma_q (cos(qz_i) - 1) \} }
More details can be found in McMillen and Dombrow (2001).

Repeat sales estimates are sometimes very sensitive to sales from the first few time periods, particularly when the sample size is small.
The option \emph{mergefirst} indicates the number of time periods for which the price index is constrained to equal zero.  The default is 
\emph{mergefirst = 1}, meaning that the price index equals zero for just the first time period.  The \emph{repsalefourier} command does not have an 
option for including an intercept in the model.  

Following Case and Shiller (1987), many authors use a three-stage procedure to construct repeat sales price indexes that are adjusted for
heteroskedasticity related to the length of time between sales.  Common specifications for the second-stage function are 
\eqn{e^2 = \alpha_0 + \alpha_1 (t-s)}{e^2 = \alpha0 + \alpha1 (t-s)} or \eqn{|e| = \alpha0 + \alpha1 (t-s)}, where \emph{e} represents the first-stage residuals.  
The first equation implies an error variance of \eqn{\hat{\sigma^2} = \hat{e^2}}{\sigma^2 = e^2} and the second equation leads to 
\eqn{\hat{\sigma^2} = \hat{|e|}^2.}{\sigma^2 = |e|^2.}  The \emph{repsalefourier} function uses a standard \emph{F} test to 
determine whether the slope cofficients are significant in the second-stage regression.  The results are reported if \emph{print=T}.

The third-stage equation is 

\deqn{ \frac{y_t - y_s}{\hat{\sigma}} = \frac{g(T_i) - g(T_i^s)}{\hat{\sigma}} + \frac{u_t - u_s}{\hat{\sigma}}}{ (y(t) - y(s))/\sigma = (g(T) - g(T_s))/\sigma + (u(t) - u(s))/\sigma}

This equation is estimated by regressing \eqn{y_t - y_s}{y(t) - y(s)} on \eqn{z, z^2, sin(z)...sin(Qz), cos(z)...cos(Qz)} 
using the \emph{weights} option in \emph{lm} with \emph{weights} = \eqn{1/\hat{\sigma^2}}{1/sigma^2}

}

\references{
Case, Karl and Robert Shiller, "Prices of Single-Family Homes since 1970:  New Indexes for Four Cities," \emph{New England Economic Review} (1987), 45-56.

McMillen, Daniel P. and Jonathan Dombrow, "A Flexible Fourier Approach to Repeat Sales Price Indexes," \emph{Real Estate Economics} 29 (2001), 207-225.

}

\seealso{
  \link{repsale}

  \link{repsaledata}

  \link{repsaleqreg}
}

\keyword{Repeat Sales}
\keyword{Series Expansions}

