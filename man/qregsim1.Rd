\name{qregsim1 }
\alias{qregsim1 }
\title{Changes in Distributions Implied by Quantile Regression Estimates}
\description{Uses quantile regression results to simulate the effects of an explanatory variable on the distribution of the dependent variable. }
\usage{
qregsim1(formall, formx, bmat, taumat, xvalues=NULL, ytarget=NULL, 
  xcolors=NULL,graphx=TRUE, graphy=TRUE, graphsim=TRUE, histogram=FALSE, 
  histfreq=FALSE, yname=NULL, xname=NULL, nsim=0, bwadjust=1, 
  legloc="topright",data=NULL) 
}

\arguments{
  \item{formall}{Formula with the dependent variable and \emph{all} explanatory variables, as provided to the \emph{qregbmat} or \emph{qregcpar} command.}
  \item{formx}{The explanatory variable to be analyzed.  Specified as a formula, e.g., formx=~x.}
  \item{bmat}{Matrix of coefficient estimates from the \emph{qregbmat} or \emph{qregcpar} command.}
  \item{taumat}{The vector of quantile values represented in \emph{bmat}.}
  \item{xvalues}{Vector of explanatory variable values for discrete changes.  If \emph{xvalues = NULL}, the values are set to 
\emph{xvalues <- quantile(x,c(.25,.75))}. }
  \item{ytarget}{Vector of target values for the density functions involving y.  The values in \emph{ytarget} are used as the target points for the 
predicted and actual values of y.  The same values of \emph{ytarget} are also used as the target values for the simulations.}

  \item{xcolors}{Vector of colors for the density function graphs.  Default for two values of \emph{xvalues} is c("black","red").  For more than two values,
the default is drawn from the \emph{RColorBrewer} package:  \emph{xcolors = brewer.pal(nx,"Blues")}. }
  \item{graphx}{If \emph{graphx=T}, presents the kernel density function for the explanatory variable.}
  \item{graphy}{If \emph{graphy=T}, presents density functions for the actual and predicted values of the dependent variable.}
  \item{graphsim}{If \emph{graphsim=T}, presents graphs of the density functions for the predicted values of y at the values specified in \emph{xvalues}.}
  \item{histogram}{If \emph{histogram=T} and \emph{graphx=T}, the density function for the explanatory variable is presented as a histogram. 
Not relevant if \emph{graphx=F}.}
  \item{histfreq}{If \emph{histogram=T} and \emph{graphx=T}, the histogram is presented using frequencies rather than densities.
Not relevant if \emph{graphx=F} or \emph{histogram=F}.}
  \item{yname}{A label used for the dependent variable in the density graphs, e.g., \emph{yname = "Log of Sale Price"}.}
  \item{xname}{A label for graphs involving the explanatory variable, e.g., \emph{xname = "x1"}.}
  \item{nsim }{Number of simulations for quantile distributions.  Default:  \emph{nsim=0}. }
  \item{bwadjust}{Factor used to adjust bandwidths for kernel densities.  Smoother functions are produced when \emph{bwadjust}>1.  
Passed directly to the \emph{density} function's \emph{adjust} option.  Default:  \emph{bwadjust=1}.}
  \item{legloc}{The legend location.}
  \item{data}{A data frame containing the data. Default:  use data in the current working directory.}
}


\value{
  \item{ytarget}{The values for the horizontal axis of the quantile predictions at \emph{xvalues}. }
  \item{densyhat}{Matrix of predictions for the predicted values of y at the values specified in \emph{xvalues}.  The dimension of the matrix
is \emph{n x length(xvalues)}.}
  \item{densy1}{The values of the density function for the actual values of the dependent variable.  \emph{densy1=NULL} if \emph{graphy=F}.}
  \item{densy2}{The values of the density function for the quantile regression predictions of the dependent variable.  \emph{densy2=NULL} if \emph{graphy=F}.}
}
\details{
The conditional quantile function is \eqn{y(\tau) = \alpha(\tau) + \beta(\tau)*x + \lambda(\tau)*z}.  The complete model specification is listed in
\emph{formall}, while x is specified in \emph{formx}, e.g., \emph{formall <- y~x+z} and \emph{formx <- ~x}. When \emph{nsim = 0}, the \emph{qregsim1} command simply
calculates predicted values of y at each value of x listed in \emph{xvalues} and at each \eqn{\tau} list in \emph{taumat}.  
Thus, the first column of \emph{densyhat} holds the estimated density function for 
\eqn{\hat{y}(\tau) = \hat{\alpha}(\tau) + \hat{\beta}(\tau)*xvalues[1] + \hat{\lambda}(\tau)*z}, 
the second column holds the predictions at \emph{xvalues[2]}, and so on.  The estimates are evaluated at each value of \eqn{\tau}, 
which leads to an \emph{n x length(taumat)} set of predictions for each value of \emph{xvalues}.  Kernel density estimates are then calculated for these predictions.  


The \emph{qregsim1} command can follow either \emph{qregbmat} or \emph{qregcpar}. All that differs is the dimension of \emph{bmat}:  using \emph{qregcpar},
the coefficients vary by observation.


If \emph{nsim > 0}, the \emph{qregsim1} function uses a simulation procedure based on the Machado-Mata (2005) approach to simulate the effect of an 
explanatory variable on the distribution of the dependent variable.  The function begins by drawing (with replacement) \emph{nsim} values from the 
rows of the explanatory variable matrix and \emph{nsim} values of \eqn{\tau}.  With \emph{nsim} values of both the explanatory variables and coefficient vectors,
the predictions are simply \eqn{\hat{y}_j(\tau) = \hat{\alpha}(\tau) + \hat{\beta}(\tau)*xvalues[j] + \hat{\lambda}(\tau)*z} for
\eqn{j = 1,...,length(xvalues)}.  


The \emph{n x length(xvalues) matrix \emph{densyhat}} holds the full set of predictions.



In addition to the predicted values of \eqn{y} at \emph{xvalues}, the \emph{qregsim1} command can produce the following:


1. The estimated density function for the explanatory variable (if \emph{graphx=T}). Presented as a histogram if \emph{histogram=T}.

2. The estimated quantile regression coefficients for the variable specified by \emph{formx}.  This graph is produced if \emph{graphb=T}. 

}


\references{
Koenker, Roger.  \emph{Quantile Regression}.  New York:  Cambridge University Press, 2005.

Machado, J.A.F. and Mata, J., "Counterfactual Decomposition of Changes in Wage Distributions using Quantile Regression," 
\emph{Journal of Applied Econometrics} 20 (2005), 445-465. 

McMillen, Daniel P., "Changes in the Distribution of House Prices over Time:  Structural Characteristics, Neighborhood or Coefficients?" 
\emph{Journal of Urban Economics} 64 (2008), 573-589. 

}

\examples{

par(ask=TRUE)
data(matchdata)
matchdata$age <- matchdata$year - matchdata$yrbuilt 
tvect <- seq(.10,.90,.10)
bmat <- qregbmat(lnprice~lnland+lnbldg+age+factor(year), data=matchdata,
  graph.factor=TRUE,taumat=tvect,graphb=FALSE)

xvect <- signif(quantile(matchdata$lnland,c(.25,.75)),4)
fit <- qregsim1(lnprice~lnland+lnbldg+age+factor(year),~lnland, bmat, 
  tvect,xvalues=xvect,data=matchdata)
#Simulation for 1995
smalldata <- matchdata[matchdata$year==1995,]
bmat95 <- bmat[,1:4]
fit <- qregsim1(lnprice~lnland+lnbldg+age,~lnland,bmat95,tvect, 
  xvalues=xvect, data=smalldata)
#Simulation for 2005
smalldata <- matchdata[matchdata$year==2005,]
bmat05 <- bmat[,1:4]
bmat05[,1] <- bmat05[,1] + bmat[,5]
fit <- qregsim1(lnprice~lnland+lnbldg+age,~lnland,bmat05,tvect, 
  xvalues=xvect, data=smalldata)

}




\seealso{
  \code{\link{qregbmat}}

  \code{\link{qregsim2}}

  \code{\link{qregcpar}}

}


\keyword{Quantile Regression}
\keyword{Machado-Mata}

