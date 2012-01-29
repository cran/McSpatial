\name{matchqreg}
\alias{matchqreg}
\title{Sample quantiles and means over time for a matched sample data set}
\description{Calculates and graphs sample means and quantiles over time.  
Intended for but not limited to a data set constructed with \emph{matchprop} or \emph{matchmahal} }

\usage{
matchqreg(form,taumat=c(.10,.25,.50,.75,.90), qreglwr.smooth=TRUE, 
  window=.50,bandwidth=0,kern="tcub", alldata=FALSE,
  graph.yhat=TRUE,graph.mean=TRUE,data) 

}

\arguments{
  \item{form}{ A formula of the type \emph{y~x}, where \emph{x} represents time. }
  \item{taumat}{Vector of quantiles.   Default:  \emph{taumat}=c(.10, .25, .50, .75, .90). }
  \item{qreglwr.smooth}{If \emph{qreglwr.smooth=T}, uses \emph{qreglwr} to smooth the quantile series.  
     If \emph{qreglwr.smooth=F}, calculates period by period quantiles.}
  \item{window}{Window size to be passed to \emph{qreglwr} if \emph{qreglwr.smooth=T.}  Default:  0.50.}
  \item{bandwidth}{Bandwidth to be passed to \emph{qreglwr} if \emph{qreglwr.smooth=T.}  Default:  0, i.e., not used.}
  \item{kern}{Kernel weighting function to be passed to \emph{qreglwr} if \emph{qreglwr.smooth=T.}   Default is the tri-cube. 
Options include "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{alldata}{ Indicates how the \emph{alldata} option should be treated for \emph{qreglwr} if \emph{qreglwr.smooth=T}.  Default:  \emph{alldata=F} } 
  \item{graph.yhat}{If \emph{graph.yhat=T}, graphs the series of quantile lines. Default:  \emph{graph.yhat=T.}}
  \item{graph.mean}{If \emph{graph.mean=T}, graphs the means over time.  Default:  \emph{graph.yhat=T.}}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}


\value{
\item{yhat}{Matrix of quantiles for y; actual quantiles if \emph{qreglwr.smooth=F} and smoothed values if \emph{qreglowr.smooth=T.} 
Rows represent time periods and columns represent quantiles.}
\item{ymean}{Average value of y for each time period.}
\item{timevect}{Vector of target quantile values.}
}


\examples{
set.seed(189)
n = 500
# sale dates range from 0-10
# mean and variance of x increase over time, from 1 to 2
# price index for y increases from 0 to 1
timesale <- array(0,dim=n)
x <- rnorm(n,0,1)
for (j in seq(1,10)) {
  timesale <- c(timesale, array(j, dim=n))
  x <- c(x, rnorm(n,j/10,1+j/10))
}
n = length(x)
y <- x*1 + timesale/10 + rnorm(n, 0, sd(x)/2)
fit <- lm(y~x+factor(timesale))
summary(fit)
heddata <- data.frame(y,x,timesale)
summary(heddata)

par(ask=TRUE)
matchdata <- matchprop(timesale~x,data=heddata,ytreat=0,
  distance="logit",discard="both")
table(matchdata$timesale)
fit <- matchqreg(y~timesale,qreglwr.smooth=FALSE, 
  graph.yhat=TRUE,graph.mean=TRUE,data=matchdata) 

}


\details{

Calculates means and quantiles of \emph{y} for each time period present in the variable on the right hand side of the model formula.
The quantiles can be varied with the \emph{taumat} option.  If \emph{qreglwr.smooth=T}, matchqreg uses the qreglwr command to smooth
the quantile lines and stores the results in the matrix \emph{yhat}.  The unsmoothed, actual quantile values are stored in \emph{yhat} if
\emph{qreglwr.smooth=F}.  The \emph{window, bandwidth, kern,} and \emph{alldata} options are passed on to \emph{qreglwr} if \emph{qreglwr.smooth=T}.  

Although \emph{matchqreg} is meant to follow the \emph{matchprop} or \emph{matchmahal} command, it can be applied to any data set. 

}

\references{
Deng, Yongheng, Sing Tien Foo, and Daniel P. McMillen, "Private Residential Price Indices in Singapore," \emph{Regional Science and Urban Economics}, (forthcoming).

Ho, D., Imai, K., King, G, Stuart, E., "Matching as Nonparametric Preprocessing for Reducing Model Dependence in Parametric Causal Inference," \emph{Political Analysis} 15 (2007), 199-236.

Ho, D., Imai, K., King, G, Stuart, E., "MatchIt:  Nonparametric preprocessing for parametric causal inference," \emph{Journal of Statistical Software} 42 (2011), 1-28..

McMillen, Daniel P., "Repeat Sales as a Matching Estimator," \emph{Real Estate Economics} (forthcoming).


}


\seealso{
  \code{\link{matchmahal}}

  \code{\link{matchprop}}

  \code{\link{qreglwr}}

}


\keyword{Matching}
\keyword{Quantile Regression}
