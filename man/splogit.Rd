\name{splogit }
\alias{splogit }
\title{Linearized GMM Spatial Logit }
\description{Implements the Klier-McMillen (2008) linearized GMM logit model for a 0-1 dependent variable 
and an underlying latent variable of the form \eqn{Y^* = \rho WY^* + X \beta +u} }
\usage{ 
splogit(form,inst=NULL,winst=NULL,wmat=NULL,shpfile=NULL,data=NULL,silent=FALSE) }

\arguments{
  \item{form }{Model formula}
  \item{inst}{List of instruments \emph{not} to be pre-multiplied by \emph{W}.  Entered as \emph{inst=~w1+w2 ...}
   Default:  \emph{inst=NULL}.  See \emph{details} for more information.}
  \item{winst}{List of instruments to be pre-multiplied by \emph{W} before use.  Entered as \emph{winst=~w1+w2 ...}  
Default:  \emph{inst=NULL}.  See \emph{details} for more information.}
  \item{wmat}{Directly enter \emph{wmat} rather than creating it from a shape file.  Default:   not specified.  
One of the \emph{wmat} or \emph{shpfile} options must be specified.}
  \item{shpfile}{Shape file to be used for creating the \emph{W} matrix.  Default:   not specified.  One of the \emph{wmat} or \emph{shpfile} options must be specified.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
  \item{silent }{If \emph{silent=T}, no output is printed}
}


\value{
  \item{coef }{Coefficient estimates }
  \item{se }{Standard error estimates }
}
\details{
The linearized model is a three-step estimation procedure.  Let \emph{y} be the indicator value:  \emph{y} = 1 when \eqn{y^* > 0}{y* > 0} and \emph{y} = 0 
when \eqn{y^* < 0}{y* < 0}.  
The first stage is standard logit of \emph{y} on \emph{X}.  The probability estimates from this regression are \eqn{p = exp(X \beta)/(1+exp(X \beta))}.    
The second/third stage of the procedure is standard 2SLS estimation of \eqn{u = y-p+gX \hat{\beta}}{u = y-p+gX \beta } 
on \eqn{gX} and \eqn{gWX \hat{\beta}}{gWX \beta} using \emph{Z} as  instruments. 
\emph{g} is the gradient vector, \eqn{dp/d\beta}. 
The covariance matrix (equation 3 in Klier-McMillen, 2008)  is estimated using the \emph{car} package.  
The final estimates minimize \eqn{(y-p)'Z(Z'Z)^{-1}Z'(y-p)} with \emph{p} linearized around \eqn{\hat{\beta}}{\beta-logit} and \emph{p} = 0.  \cr

\emph{splogit} provides flexibility in specifying the list of instruments.  
By default, the instrument list includes \emph{X} and \emph{WX}, where \emph{X} is the original explanatory variable list and \emph{W} is the spatial weight matrix.  
It is also possible to directly specify the full instrument list or to include only a subset of the \emph{X} variables in the list that is to be pre-multiplied by \emph{W}.  \cr

Let \emph{list1} and \emph{list2} be user-provided lists of the form \emph{list=~z1+z2}.  
The combinations of defaults (\emph{NULL}) and lists for \emph{inst} alter the final list of instruments as follows: \cr

\emph{inst = NULL},  \emph{winst = NULL}:   \emph{Z = (X, WX)} \cr
\emph{inst = list1}, \emph{winst = NULL}:   \emph{Z = list1} \cr
\emph{inst = NULL},  \emph{winst = list2}:  \emph{Z = (X, W*list2)} \cr
\emph{inst = list1}, \emph{winst = list2}:  \emph{Z = (list1, W*list2)} \cr

Note that when \emph{inst=list1} and \emph{winst=NULL} it is up to the user to specify at least one variable in \emph{list1} that is not also included in \emph{X}.
 }
\references{
Klier, Thomas and Daniel P. McMillen, "Clustering of Auto Supplier Plants in the United States:  Generalized Method of Moments Spatial Logit for Large Samples," \emph{Journal of 
Business and Economic Statistics} 26 (2008), 460-471.
 }

\seealso{
  \link{cparlogit}

  \link{cparprobit}

  \link{cparmlogit}

  \link{gmmlogit}

  \link{gmmprobit}

  \link{spprobit}

  \link{spprobitml}


}

\examples{
set.seed(9947)
library(maptools)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
wmat <- makew(cmap)$wmat
n = nrow(wmat)
rho = .4
x <- runif(n,0,10)
ystar <- as.numeric(solve(diag(n) - rho*wmat)\%*\%(x + rnorm(n,0,2)))
y <- ystar>quantile(ystar,.4)
fit <- splogit(y~x,  wmat=wmat)
 }


\keyword{Spatial AR Model}
\keyword{Parametric Models}
\keyword{Discrete Choice Models}
