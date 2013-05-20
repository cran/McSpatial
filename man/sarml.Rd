\name{sarml}
\alias{sarml}
\title{
Spatial AR Maximum-Likelihood Estimation
}
\description{
Estimates the model \eqn{Y = \rho WY + X \beta + u} by maximizing the log-likelihood function.
}
\usage{
sarml(form,wmat=NULL,shpfile=NULL,wy=NULL,eigvar=NULL,startrho=NULL,
  print=TRUE,data=NULL)
}

\arguments{
  \item{form }{Model formula}
  \item{wmat }{The \emph{W} matrix.  If not specified, \emph{W} will be calculated from the shape file.  Default:  \emph{W} = NULL.}
  \item{shpfile }{Shape file.  Needed unless (1) \emph{wy} and \emph{eigvar} are both provided, or (2) \emph{wmat} and \emph{eigvar} are provided}
  \item{wy }{The \emph{WY} variable.  Default:  not specified; program attempts to calculate \emph{WY} using \emph{wmat} or \emph{shpfile}.}
  \item{eigvar }{The vector of eigenvalues for \emph{W}.  Default:  not provided.  
\emph{shpfile} must be specified to calculate the eigenvalues within the \emph{sarml} command.}
  \item{startrho }{A starting value for \eqn{\rho}.  Default:  \emph{startrho=0.}  Estimation will generally be faster if \emph{startrho = 0.}}
  \item{print }{If \emph{print=F}, no results are printed.  Default:  \emph{print=T.}}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}
\details{
The primary motivation for the \emph{sarml} command is to provide a convenient way to estimate multiple spatial AR models 
without having to calculate the eigenvalues of \emph{W} each time.  Under the assumption that the errors, \emph{u}, are 
independently and identically distributed normal, the log-likelihood function for the spatial AR model is  
\deqn{lnl = -\frac{n}{2}log(\pi) - \frac{n}{2}log(\sigma^2) - \frac{1}{2 \sigma^2}\sum_i u_i^2 + \sum_ilog(1 - \rho*eigvar_i)}{lnl = -n*log(pi)/2 - n*log(s2)/2 - \sum(u^2)/(2*s2) + \sum log(1 - rho*eigvar)} 
where \emph{eigvar} is the vector of eigenvalues of \emph{W}.   
Though \emph{spdep} provides a convenient and fast method for calculating the eigenvalues from a shape file, 
the calculation can nonetheless take a very long time for large data sets.  
The \emph{sarml} command allows the user to input the vector of eigenvalues directly, 
which saves time when several models are estimated using the same \emph{W} matrix.   
Unless a vector of eigenvalues is provided using the \emph{eigvar} option, 
the eigenvalues are calculated from the shape file (provided using the \emph{shpfile} option) using the \emph{spdep} package.  




Conditional on the value of \eqn{\rho}, the maximum likelihood estimate of \eqn{\beta} is simply the vector of coefficients from a regression of \eqn{Y - \rho WY} on X.  
The estimate of the error variance also has a closed form solution:  \eqn{\hat{\sigma^2} = \sum_iu_i^2/n}{\sigma^2 = \sum u^2/n}. 
Substituting these estimates into the log-likelihood functions leads to the following concentrated log-likelihood function:

\deqn{lc = -\frac{n}{2}(log(\pi)+1) - \frac{n}{2}log(\sum_iu_i^2) + \sum_i log(1-\rho*eigvar_i)}{lc = -n*(log(pi)+1)/2 - n*log(\sum(u^2))/2 + \sum log(1-rho*eigvar)}


Working with the concentrated likelihood function reduces the optimization problem to a one-dimensional search for the value of \eqn{\rho} that maximizes \emph{lc}.  
Unless a value is provided for \emph{startrho}, the \emph{sarml} procedure begins by using the \emph{optimize} command to find the value of \eqn{\rho} 
that maximizes \emph{lc}.  
This estimate of \eqn{\rho} (or the value provided by the \emph{startrho} option) is then used to calculate the implied values of \eqn{\beta} and \eqn{\sigma^2}, 
and these values are used as starting values to maximize \emph{lnl} using the \emph{nlm} command. \cr



The covariance matrix for the estimates of \eqn{\beta} and \eqn{\rho}, \emph{vmat}, is the inverse of \eqn{(1/\sigma^2)V}.
\emph{V} has partitions \eqn{V_{11} = X'X}{V11 = X'X},   \eqn{V_{12} = X'WY}{V12 = X'WY},    \eqn{V_{21} = Y'W'X}{V21 = Y'W'X},   
and \eqn{V_{22} = Y'W'WY + \sigma^2 \sum_i \frac{eigvar_i}{(1 - \rho*eigvar_i)^2 }}{V22 = Y'W'WY + s2*\sum (eigvar/(1-rho*eigvar))^2}.

}
\value{
  \item{beta }{The estimated vector of coefficients, \eqn{\beta}.}
  \item{rho }{The estimated value of \eqn{\rho}.}
  \item{sig2 }{The estimated error variance, \eqn{\sigma^2}.}
  \item{vmat }{The covariance matrix for \eqn{(\beta, \rho^2)}.}
  \item{eigvar }{The vector of eigenvalues}
}


\examples{
library(spdep)

cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
samppop <- cmap$POPULATION>0&cmap$AREA>0
cmap <- cmap[samppop,]
cmap$lndens <- log(cmap$POPULATION/cmap$AREA)
lmat <- coordinates(cmap)
cmap$LONGITUDE <- lmat[,1]
cmap$LATITUDE  <- lmat[,2]
cmap$dcbd <- geodistance(longvar=cmap$LONGITUDE,latvar=cmap$LATITUDE,
  lotarget=-87.627800,latarget=41.881998)$dist

fit <- makew(shpfile=cmap,eigenvalues=TRUE)
wmat <- fit$wmat
eigvar <- fit$eigvar

# input w, calculate eigvar within sarml
fit <- sarml(lndens~dcbd,wmat=wmat,eigvar=eigvar,data=cmap)

}
\seealso{
  \code{\link{makew}}

  \code{\link{qregspiv}}

  \code{\link{qregbmat}}

  \code{\link{qregsim1}}

  \code{\link{qregsim2}}

  \code{\link{qregcpar}}

  \code{\link{qreglwr}}

}

\keyword{Spatial AR Model}
\keyword{Parametric Models}


