\name{gmmprobit }
\alias{gmmprobit }
\title{GMM Spatial Probit }
\description{Estimates a GMM probit model for a 0-1 dependent variable 
and an underlying latent variable of the form \eqn{Y^* = \rho WY^* + X \beta +u} }
\usage{
gmmprobit(form,inst=NULL,winst=NULL,wmat=NULL,shpfile,
  startb=NULL,startrho=0,blockid=0,cvcrit=.0001,data=NULL,silent=FALSE) 
  }

\arguments{
  \item{form}{Model formula}
  \item{inst}{List of instruments \emph{not} to be pre-multiplied by \emph{W}.  Entered as \emph{inst=~w1+w2 ...}
   Default:  \emph{inst=NULL}.  See \emph{details} for more information.}
  \item{winst}{List of instruments to be pre-multiplied by \emph{W} before use.  Entered as \emph{winst=~w1+w2 ...}  
Default:  \emph{inst=NULL}.  See \emph{details} for more information.}
  \item{wmat}{Directly enter \emph{wmat} rather than creating it from a shape file.  Default:   not specified.  
One of the \emph{wmat} or \emph{shpfile} options must be specified.}
  \item{shpfile}{Shape file to be used for creating the \emph{W} matrix.  Default:   not specified.  One of the \emph{wmat} or \emph{shpfile} options must be specified.}
  \item{startb}{Vector of starting values for \emph{B}.  Default:  use estimates from \emph{spprobit}, the linearized version of the model.  
Specified as \emph{startb=0.}}
  \item{startrho}{Vector of starting values for \eqn{\rho}.  Default:  use estimates from \emph{spprobit}, the linearized version of the model.  
Specified as \emph{startrho=0.}}
  \item{blockid}{A variable identifying groups used to specify a block diagonal structure for the \emph{W} matrix, e.g., \emph{blockid=state} or \emph{blockid=region}.  
Imposes that all elements outside of the blocks equal zero and then re-standardizes \emph{W} such that the rows sum to one. 
By default, \emph{blockid = 0}, and a block diagonal structure is \emph{not} imposed.}
  \item{cvcrit}{Convergence criterion.  Default:  \emph{cvcrit} = 0.0001.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
  \item{silent }{If \emph{silent=T}, no output is printed}
}


\value{
  \item{coef }{Coefficient estimates }
  \item{se }{Standard error estimates }
}
\details{
The underlying latent variable for the model is \eqn{Y^* = \rho WY^* + X \beta + u}{Y* = \rho WY* + X \beta + u} 
or \eqn{Y^* = (I - \rho W)^{-1}(X \beta + u)}{Y* = (I - \rho W)^{-1}(X \beta + u)}.  
The covariance matrix is \eqn{Euu' = \sigma^2 ((I - \rho W)(I - \rho W)')^{-1}}, with \eqn{\sigma^2} normalized to unity.  
Typical specifications imply heteroskedasticity, i.e., the diagonal elements of \eqn{Euu'}, denoted by \eqn{\sigma_i^2}, vary across observations.  
Heteroskedasticity makes standard probit estimates inconsistent.  Letting \eqn{X_i^* = X_i/ \sigma_i} and \eqn{H = (I - \rho W)^{-1} X^*}, 
the probit probabilities implied by the latent variable are \eqn{p = \Phi(H \beta)} and the generalized error term is \eqn{e = (y - p)\phi(H \beta)/(p(1-p)) }, 
where \eqn{y = 1} if \eqn{Y^* >0} and \eqn{y = 0} otherwise. \cr

The GMM estimator chooses \eqn{\beta} and \eqn{\rho} to minimize \eqn{e'Z(Z'Z)^{-1}Z'e}, 
where \emph{Z} is a matrix of instruments specified using the \emph{inst} and \emph{winst} options. 
Unless specified otherwise using the \emph{startb} and \emph{startrho} options, initial estimates are obtained using \emph{spprobit}, 
which implements the simple (and fast) linearized version of the GMM probit model proposed by Klier and McMillen (2008).  
Convergence is defined by \eqn{abs(change) < cvcrit}, where \emph{change} is the gradient vector implied by applying a standard Gauss-Newton algorithm to the objective function.  
The covariance matrix (equation 3 in Klier-McMillen, 2008)  is estimated using the \emph{car} package.  \cr

Estimation can be very slow because each iteration requires the inversion of an \emph{nxn} matrix.  
To speed up the estimation process and to reduce memory requirements, it may be desirable to impose a block diagonal structure on \emph{W}.  
For example, it may be reasonable to impose that each state or region has its own error structure, with no correlation of errors across regions.  
The \emph{blockid} option specifies a block diagonal structure such as \emph{blockid=region}.  
The option leads the program to re-calculate the \emph{W} matrix, imposing the block diagonal structure and re-normalizing the matrix to again have each row sum to one.  
If there are \emph{G} groups, estimation requires \emph{G} sub-matrices to be inverted rather than one \emph{nxn} matrix, 
which greatly reduces memory requirements and significantly reduces the time required in estimation. \cr


\emph{gmmprobit} provides flexibility in specifying the list of instruments.  
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
Klier, Thomas and Daniel P. McMillen, "Clustering of Auto Supplier Plants in the United 
States:  Generalized Method of Moments Spatial probit for Large Samples," \emph{Journal of 
Business and Economic Statistics} 26 (2008), 460-471. \cr

Pinkse, J. and M. E. Slade, "Contracting in Space:  An Application of Spatial Statistics to 
Discrete-Choice Models," \emph{Journal of Econometrics} 85 (1998), 125-154. \cr
 }

\seealso{
  \link{cparlogit}

  \link{cparprobit}

  \link{cparmlogit}

  \link{gmmlogit}

  \link{splogit}

  \link{spprobit}

}

\examples{
set.seed(9947)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
lmat <- coordinates(cmap)
dnorth <- geodistance(lmat[,1],lmat[,2], -87.627800, 
	41.881998, dcoor=TRUE)$dnorth
cmap <- cmap[dnorth>0,]
wmat <- makew(cmap)$wmat
n = nrow(wmat)
rho = .4
x <- runif(n,0,10)
ystar <- as.numeric(solve(diag(n) - rho*wmat)\%*\%(x + rnorm(n,0,2)))
y <- ystar>quantile(ystar,.4)
fit <- gmmprobit(y~x,  wmat=wmat)
}


\keyword{Parmametric Models}
\keyword{Spatial AR Model}
\keyword{Discrete Choice Models}

