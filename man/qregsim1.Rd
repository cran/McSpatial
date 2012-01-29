\name{qregsim1 }
\alias{qregsim1 }
\title{Changes in Distributions Implied by Quantile Regression Estimates}
\description{Uses quantile regression results to simulate the effects of an explanatory variable on the distribution of the dependent variable. }
\usage{
qregsim1(formall, formx, bmat, xmin=NULL, xmax=NULL, 
  graphx=TRUE, graphb=TRUE, graphy=TRUE, graphxb=TRUE, graphols=FALSE,
  histogram=FALSE, histfreq=FALSE, yname=NULL, xname=NULL, nsim=20000, 
  bwadjust=1, data=NULL)  
}

\arguments{
  \item{formall}{Formula with the dependent variable and \emph{all} explanatory variables, as provided to the \emph{qregbmat} command.}
  \item{formx}{The explanatory variable to be analyzed.  Specified as a formula, e.g., formx=~x.}
  \item{bmat}{Matrix of coefficient estimates from the \emph{qregbmat} command.}
  \item{xmin}{Lower value of the explanatory variable for discrete changes.}
  \item{xmax}{Upper value of the explanatory variable for discrete changes.}
  \item{graphx}{If \emph{graphx=T}, presents the kernel density function for the explanatory variable.}
  \item{graphb}{If \emph{graphb=T}, presents graph of the estimated quantile coefficients.}
  \item{graphy}{If \emph{graphy=T}, presents density functions for the actual and predicted values of the dependent variable.}
  \item{graphxb}{If \emph{graphxb=T}, presents graphs of the counterfactual predictions.}
  \item{graphols}{If \emph{graphols=T}, OLS density estimates are added to all graphs.}
  \item{histogram}{If \emph{histogram=T} and \emph{graphx=T}, the density function for the explanatory variable is presented as a histogram. 
Not relevant if \emph{graphx=F}.}
  \item{histfreq}{If \emph{histogram=T} and \emph{graphx=T}, the histogram is presented using frequencies rather than densities.
Not relevant if \emph{graphx=F} or \emph{histogram=F}.}
  \item{yname}{A label used for the dependent variable in the density graphs, e.g., \emph{yname = "Log of Sale Price"}.}
  \item{xname}{A label for graphs involving the explanatory variable, e.g., \emph{xname = "x1"}.}
  \item{nsim }{Number of simulations for quantile distributions.  Default:  \emph{nsim=20000}. }
  \item{bwadjust}{Factor used to adjust bandwidths for kernel densities.  Smoother functions are produced when \emph{bwadjust}>1.  
Passed directly to the \emph{density} function's \emph{adjust} option.  Default:  \emph{bwadjust=1}.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}


\value{
  \item{xtarget}{The values for the horizontal axis of the quantile counterfactual density graphs. }
  \item{densx1}{Counterfactual density function. When \emph{xmin=NULL} and \emph{xmax=NULL}, \emph{densx1} is the density function for the predicted
values, \eqn{X \beta}, when X is drawn from the actual values and \eqn{\beta} is drawn \emph{bmat}.  
Otherwise, \emph{densx1} is the density function when the x specified by \emph{formx} is set to \emph{xmin} while \eqn{\beta} values are drawn from \emph{bmat}
and all other X variables are drawn from their actual values.}
  \item{densx2}{Counterfactual density function. When \emph{xmin=NULL} and \emph{xmax=NULL}, \emph{densx2} is the density function for the predicted
values, \eqn{X \beta}, when the x specified by \emph{formx} is set to its mean value while \eqn{\beta} values are drawn from \emph{bmat}
and all other X variables are drawn from their actual values.
Otherwise, \emph{densx2} is the density function when when the x specified by \emph{formx} is set to \emph{xmax}.}
  \item{ytarget}{The values for the horizontal axis for the density functions for the actual and predicted values of the dependent variable.  
\emph{ytarget=NULL} if \emph{graphy=F}.}
  \item{densy1}{The values of the density function for the actual values of the dependent variable.  \emph{densy1=NULL} if \emph{graphy=F}.}
  \item{densy2}{The values of the density function for the quantile regression predictions of the dependent variable.  \emph{densy2=NULL} if \emph{graphy=F}.}
}
\details{
The \emph{qregsim1} function uses a simulation procedure based on the Machado-Mata (2005) approach to simulate the effect of an explanatory variable on the
distribution of the dependent variable.  The command is intended to follow the \emph{qregbmat} command. The following results are produced:

1. The estimated density function for the explanatory variable (if \emph{graphx=T}). Presented as a histogram if \emph{histogram=T}.

2. The estimated quantile regression coefficients for the variable specified by \emph{formx}.  This graph is produced if \emph{graphb=T}. 

3. Graphs of the density functions for the actual values of the dependent variable and the predicted values from the OLS (if \emph{olsgraph=T}) and quantile regressions.
Using a procedure based on Machado-Mata (2005), the X matrix used for the quantile predictions is based on \emph{nsim} random draws with replacement from the
original matrix of explanatory variables.  Similarly, the coefficients are based on \emph{nsim} draws with replacement from the matrix of estimated quantile 
regression coefficients, which are supplied by \emph{bmat}.  The OLS estimates are based on the actual values of X rather than on simulated values.

4. The counterfactual density functions.  To illustrate how these functions are calculated, consider a simple case in which the full model estimated by
\emph{qregmat} has two explanatory variables, x1 and x2, and suppose that \emph{formx=~x1} for the \emph{qregsim1} command. 
The density functions differ depending on whether values are specified for \emph{xmin} and \emph{xmax}:


Case 1:  \emph{xmin=NULL} and \emph{xmax=NULL}. In this case, the base set of predicted values of y1 = x1*b1 + x2*b2 is calculated using \emph{nsim} draws 
with replacement from the original set of values for x1 and x2 and from the quantile regression estimates of b1 and b2. 
The alternative (i.e., counterfactual) set of predicted values is
calculated using the mean value of x1 in place of the random draws, i.e., y2 = mean(x1)*b1 + x2*b2.  The coefficients, b1, are evaluated at the same
set of random draws as used for the base set of predictions.  A similar procedure is used for the OLS predictions if \emph{olsgraph=T}, but the 
explanatory variables are evaluated at actual values rather than randomly drawn values.
The density function for the base set of quantile predicted values is stored in \emph{densx1} and the density function for the quantile predictions
with x1 set to mean(x1) is stored in \emph{densx2}.  

Case 2:  Values specified for \emph{xmin} and \emph{xmax}.  In this case, the base set of predicted values is y1 = xmin*b1 + x2*b2, with x2, b1, and b2 drawn
with replacement from the original values.  The alternative set of predicted values is y2 = xmax*b1 + x2*b2. The density function for y1 is stored in
densx1 and the density function for y2 is stored in densx2.  A similar procedure is used for the OLS predictions if \emph{olsgraph=T}, but x2 is evaluated 
at actual values rather than randomly drawn values.  This procedure is comparable to calculating a discrete derivative.


Case 3:  Value specified for one of \emph{xmin} or \emph{xmax} while the other equals \emph{NULL}.  This case is equivalent to case 2 with \emph{xmin} set to x1 while
\eqn{xmax = x1 + \delta}, where \eqn{\delta} is the value provided by the non-\emph{NULL} value.  For example, if \emph{xmin}=1 and \emph{xmax = NULL}, 
the density functions are evaluated at y1 = x1*b1 + x2*b2 and y2 = (x1+1)*b1 + x2*b2.  The same results would be produced by \emph{xmin = NULL} and \emph{xmax}=1.
As in case 2, quantile estimates are evaluated at simulated values for x1 and x2 while OLS estimates are presented at actual values.
}


\references{
Koenker, Roger.  \emph{Quantile Regression}.  New York:  Cambridge University Press, 2005.

Machado, J.A.F. and Mata, J., "Counterfactual Decomposition of Changes in Wage Distributions using Quantile Regression," \emph{Journal of Applied Econometrics} 20 (2005), 445-465. 

McMillen, Daniel P., "Changes in the Distribution of House Prices over Time:  Structural Characteristics, Neighborhood or Coefficients?" \emph{Journal of Urban Economics} 64 (2008), 573-589. 

 }

\examples{

par(ask=TRUE)
data(matchdata)
matchdata$age <- matchdata$year - matchdata$yrbuilt 
bmat <- qregbmat(lnprice~lnland+lnbldg+age+factor(year), data=matchdata,
  graph.factor=TRUE)

fit <- qregsim1(lnprice~lnland+lnbldg+age+factor(year),~lnland, bmat, 
  data=matchdata)
lo = signif(quantile(matchdata$lnland,.25),4)
hi = signif(quantile(matchdata$lnland,.75),4)
#Simulation for 1995
smalldata <- matchdata[matchdata$year==1995,]
bmat95 <- bmat[,1:4]
fit <- qregsim1(lnprice~lnland+lnbldg+age,~lnland, bmat95, 
  xmin=lo, xmax=hi, data=smalldata)
#Simulation for 2005
smalldata <- matchdata[matchdata$year==2005,]
bmat05 <- bmat[,1:4]
bmat05[,1] <- bmat05[,1] + bmat[,5]
fit <- qregsim1(lnprice~lnland+lnbldg+age,~lnland, bmat05, 
  xmin=lo, xmax=hi, data=smalldata)

}




\seealso{
  \code{\link{qregbmat}}

  \code{\link{qregsim2}}

  \code{\link{qregcpar}}

  \code{\link{qreglwr}}

}


\keyword{Quantile Regression}
\keyword{Machado-Mata}

