\name{repsale}
\alias{repsale}
\title{
  Repeat Sales Estimation
}
\description{
  Standard and Weighted Least Squares Repeat Sales Estimation
}

\usage{
repsale(price0,time0,price1,time1,mergefirst=1,
 graph=TRUE,graph.conf=TRUE,conf=.95,
  stage3=FALSE,stage3_xlist=~timesale,print=TRUE) 
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
  \item{xmat}{The matrix of explanatory variables for the repeat sales regressions.  \emph{dim(xmat) = nt - mergefirst,} where
\emph{nt} = the number of time periods and \emph{mergefirst} is specified in the call to \emph{repsale.} }
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
fit <- repsale(price0=y0, price1=y1, time0=time0, time1=time1)

# variance rises with timesale
# var(u0) = .2^2; var(u1) = (.2 + timesale/10)^2
# var(u1-u0) = var(u0) + var(u1) = 2*(.2^2) + .4*timesale/10 + (timesale^2)/100
y0 <- time0/10 + rnorm(n,0,.2)
y1 <- time1/10 + rnorm(n,0,.2+timesale/10)
par(ask=TRUE)
fit <- repsale(price0=y0, price1=y1, time0=time0, time1=time1)
summary(fit$pindex)
fit <- repsale(price0=y0, price1=y1, time0=time0, time1=time1, stage3="abs")
summary(fit$pindex)
timesale2 <- timesale^2
fit <- repsale(price0=y0, price1=y1, time0=time0, time1=time1, stage3="square", 
  stage3_xlist=~timesale+timesale2)
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
\emph{mergefirst = 1}, meaning that the price index equals zero for just the first time period.  The \emph{repsale} command does not have an 
option for including an intercept in the model.  

Following Case and Shiller (1987), many authors use a three-stage procedure to construct repeat sales price indexes that are adjusted for
heteroskedasticity related to the length of time between sales.  Common specifications for the second-stage function are 
\eqn{e^2 = \alpha_0 + \alpha_1 (t-s)}{e^2 = \alpha0 + \alpha1 (t-s)} or \eqn{|e| = \alpha0 + \alpha1 (t-s)}, where \emph{e} represents the first-stage residuals.  
The first equation implies an error variance of \eqn{\hat{\sigma^2} = \hat{e^2}}{\sigma^2 = e^2} and the second equation leads to 
\eqn{\hat{\sigma^2} = \hat{|e|}^2.}{\sigma^2 = |e|^2.}  The \emph{repsale} function uses a standard \emph{F} test to 
determine whether the slope cofficients are significant in the second-stage regression.  The results are reported if \emph{print=T}.

The third-stage equation is 

\deqn{ \frac{y_t - y_s}{\hat{\sigma}} = \frac{\delta_t - \delta_s}{\hat{\sigma}} + \frac{u_t - u_s}{\hat{\sigma}}}{ (y(t) - y(s))/\sigma = (\delta(t) - \delta(s))/\sigma + (u(t) - u(s))/\sigma}

This equation is estimated by regressing \eqn{y_t - y_s}{y(t) - y(s)} on the series of indicator variables implied by \eqn{\delta_t - \delta_s}{\delta(t) - \delta(s)} 
using the \emph{weights} option in \emph{lm} with \emph{weights} = \eqn{1/\hat{\sigma^2}}{1/sigma^2}
}

\references{
Case, Karl and Robert Shiller, "Prices of Single-Family Homes since 1970:  New Indexes for Four Cities," \emph{New England Economic Review} (1987), 45-56.
}

\seealso{
  \link{repsaledata}

  \link{repsalefourier}

  \link{repsaleqreg}
}

\keyword{Repeat Sales}

