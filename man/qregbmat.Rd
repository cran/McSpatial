\name{qregbmat }
\alias{qregbmat }
\title{Quantile Regression for Multiple Quantiles}
\description{Returns estimated coefficients from a series of quantile regressions. }
\usage{ qregbmat(form,taumat=seq(.10,.90,.10),graphb=TRUE,graph.factor=FALSE,
  data=NULL) }

\arguments{
  \item{form }{Model formula}
  \item{taumat }{Vector of target quantiles.  Default:  \emph{taumat=seq(.10,.90,.10)} }
  \item{graphb}{If \emph{graphb=TRUE}, prints graphs of the coefficient estimates.  Default:  \emph{graphb=TRUE}. }
  \item{graph.factor}{If \emph{graph.factor=TRUE} and \emph{graphb=TRUE}, prints graphs of the coefficient estimates for any factor variables.  
    Default:  \emph{graph.factor=TRUE}.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}


\value{
  Returns the length(taumat) x k matrix of estimated coefficients, where k is the number of explanatory variables.
}
\details{
Estimates a series of quantile regressions using the \emph{quantreg} packages.  The quantiles are listed in \emph{taumat}.
The \emph{qregbmat} command is intended primarily as a first stage before the \emph{qregsim1} or \emph{qregsim2} commands.
 }

\references{
Koenker, Roger.  \emph{Quantile Regression}.  New York:  Cambridge University Press, 2005.
 }

\examples{

par(ask=TRUE)
data(matchdata)
matchdata$age <- matchdata$year - matchdata$yrbuilt 
bmat <- qregbmat(lnprice~lnland+lnbldg+age+factor(year), data=matchdata,
  graph.factor=TRUE)
summary(bmat)
}


\seealso{
  \code{\link{qregsim1}}

  \code{\link{qregsim2}}

  \code{\link{qregcpar}}

  \code{\link{qreglwr}}

}

\keyword{Quantile Regression}
\keyword{Machado-Mata}

