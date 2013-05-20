\name{spprobitml }
\alias{spprobitml }
\title{Maximum Likelihood Estimation of a Spatial Probit Model}
\description{Probit estimation for a model with an underlying latent variable of the form \eqn{Y^* = \rho WY^* + X \beta +u} }
\usage{
  spprobitml(form,wmat=NULL,shpfile=NULL,blockid=NULL,
  minblock=NULL,maxblock=NULL,stdprobit=TRUE,data=NULL)
}

\arguments{
  \item{form }{Model formula}
  \item{wmat}{Directly enter \emph{wmat} rather than creating it from a shape file.  Default:   not specified.  
One of the \emph{wmat} or \emph{shpfile} options must be specified.}
  \item{shpfile}{Shape file to be used for creating the \emph{W} matrix.  
Default:   not specified.  One of the \emph{wmat} or \emph{shpfile} options must be specified.
The order of the observations in \emph{wmat} must be the same as the order of observations in \emph{data}.}
  \item{blockid}{A variable identifying groups used to specify a block diagonal structure for the \emph{W} matrix, 
e.g., \emph{blockid=state} or \emph{blockid=region}.  Calculates a separate \emph{W} matrix for each block.  
The \emph{shpfile} option must be specified; \emph{wmat} is ignored.}
  \item{minblock}{Groups with fewer than \emph{minblock} observations are omitted.  Default is the number of explanatory variables, including the spatial lag term.}
 \item{maxblock}{Groups with more than \emph{maxblock} observations are omitted.  Unlimited by default.  This option may be useful for very large
data sets as full nblock x nblock matrices must be constructed for each block, where nblock is the number of observations in the block.}
  \item{stdprobit}{If \emph{TRUE}, also prints standard probit model results.  Default:  \emph{stdprobit=TRUE}.}
  \item{data}{A data frame containing the data. Default:  use data in the current working directory}

}

\value{
  \item{coef}{Coefficient estimates.}
  \item{logl}{The log-likelihood value.}
  \item{vmat1}{The covariance matrix for \eqn{\hat{\beta}}, conditional on \eqn{\hat{\rho}}.}
  \item{vmat2}{The unconditional covariance matrix for \eqn{\hat{\theta} = (\hat{\beta}, \hat{\rho})}.}
}
\details{
Estimation is based on the reduced form of the spatial AR model, \eqn{Y^* = (I - \rho W)^{-1}(X \beta + u)}.  
The model structure typically implies heteroskedasticity:  the variance of the reduced form error term, \eqn{(I - \rho W)^{-1}u}, is
\eqn{\sigma^2 diag \{ (I - \rho W)^{-1}(I - \rho W')^{-1} \} }.  For probit estimation, \eqn{\sigma^2} is normalized to one.
Let \eqn{s_i^2} denote the variance for observation \emph{i}, and define \eqn{X^* =  (I - \rho W)^{-1}X}. 
Then the probability that \eqn{Y_i^* > 0} is \eqn{\Phi (X_i^* \beta / s_i) }, and 
the log-likelihood function is \eqn{\sum_i \{ y_i ln (\Phi_i ) + (1-y_i ) ln(1-\Phi_i) \} }.  
The \emph{spprobitml} commands estimates the model by maximizing this log-likelihood function with respect to \eqn{\beta} and \eqn{\rho}.

Variants of this approach -- maximizing the log-likelihood function implied by the reduced form of the model -- were proposed by
Case (1992) and McMillen (1992).  Case's estimation procedure relies on a simple form of the spatial weight matrix
in which each observation within a district is affected equally by the other observations in the district.
McMillen's (1992) approach is equivalent to the one used here, but he suggested using an EM algorithm to estimate the model.
Neither author suggested a covariance matrix:  Case (1992) appears to have relied on the standard probit estimate which applies when 
the model is estimated conditional on \eqn{\rho}, while McMillen (1992) proposed a bootstrap approach. 

A consistent covariance matrix can be calculated using the gradient terms:

\deqn{V(\hat{\theta})^{-1} = \left( \sum_i \partial lnL_i / \partial \hat{\theta} \right)\left( \sum_i \partial lnL_i / \partial \hat{\theta}' \right) }

The gradient term for \eqn{\hat{\rho}} is calculated using numeric derivatives.
The covariance matrix, \eqn{V(\hat{\theta})}, is not fully efficient because the estimation procedure only indirectly takes into account the autocorrelation structure.  
An analogous approach is used to calculate standard errors conditional on \eqn{\hat{\rho}}.  
In the conditional case, only the gradient terms for \eqn{\hat{\beta}} are used; they are evaluated using the estimated values of \eqn{\rho}.

Estimation can be very slow because each iteration requires the inversion of an \emph{nxn} matrix.  
To speed up the estimation process and to reduce memory requirements, it may be desirable to impose a block diagonal structure on \emph{W}.  
For example, it may be reasonable to impose that each state or region has its own error structure, with no correlation of errors across regions.  
The \emph{blockid} option specifies a block diagonal structure such as \emph{blockid=region}.  
If there are \emph{G} groups, estimation requires \emph{G} sub-matrices to be inverted rather than one \emph{nxn} matrix, 
which greatly reduces memory requirements and significantly reduces the time required in estimation. 
The weight matrix must be calculated from \emph{shpfile} if the \emph{blockid} option is specified; the \emph{wmat} option should be set to \emph{NULL}.
}

\references{
Case, Anne C., "Neighborhood Influence and Technological Change," \emph{Regional Science and Urban Economics} 22 (1992), 491-508.

McMillen, Daniel P., "Probit With Spatial Autocorrelation," \emph{Journal of Regional Science} 32 (1992), 335-348.
 }

\examples{
set.seed(9947)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
lmat <- coordinates(cmap)
dnorth <- geodistance(lmat[,1],lmat[,2], lotarget=-87.627800, 
	latarget=41.881998,dcoor=TRUE)$dnorth
cmap <- cmap[dnorth>1,]
wmat <- makew(cmap)$wmat
n = nrow(wmat)
rho = .4
x <- runif(n,0,10)
ystar <- as.numeric(solve(diag(n) - rho*wmat)\%*\%(x + rnorm(n,0,2)))
y <- ystar>quantile(ystar,.4)
fit <- spprobitml(y~x,  wmat=wmat)

}


\seealso{
  \link{cparlogit}

  \link{cparprobit}

  \link{cparmlogit}

  \link{gmmlogit}

  \link{gmmprobit}

  \link{makew}

  \link{splogit}

  \link{spprobit}

}

\keyword{Spatial AR Model}
\keyword{Parametric Models}
\keyword{Discrete Choice Models}
