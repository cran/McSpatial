\name{spprobit }
\alias{spprobit }
\title{Linearized GMM spatial probit }
\description{Implements the Klier-McMillen (2008) linearized GMM probit model for a 0-1 dependent variable 
and an underlying latent variable of the form \eqn{Y^* = \rho WY^* + X \beta +u} }
\usage{ 
spprobit(form,inst=NULL,winst=NULL,wmat=NULL,shpfile=NULL,blockid=NULL,
         minblock=NULL,maxblock=NULL,data=NULL,silent=FALSE,minp=NULL)
 }

\arguments{
  \item{form }{Model formula}
  \item{inst}{List of instruments \emph{not} to be pre-multiplied by \emph{W}.  Entered as \emph{inst=~w1+w2 ...}
   Default:  \emph{inst=NULL}.  See \emph{details} for more information.}
  \item{winst}{List of instruments to be pre-multiplied by \emph{W} before use.  Entered as \emph{winst=~w1+w2 ...}  
Default:  \emph{inst=NULL}.  See \emph{details} for more information.}
  \item{wmat}{Directly enter \emph{wmat} rather than creating it from a shape file.  Default:   not specified. }
  \item{shpfile}{Shape file to be used for creating the \emph{W} matrix.  
The order of the observations in \emph{wmat} must be the same as the order of observations in \emph{data}.}
  \item{blockid}{A variable identifying groups used to specify a block diagonal structure for the \emph{W} matrix, 
e.g., \emph{blockid=state} or \emph{blockid=region}.  Calculates a separate \emph{W} matrix for each block.  
The \emph{shpfile} option must be specified; \emph{wmat} is ignored.}
  \item{minblock}{Groups with fewer than \emph{minblock} observations are omitted.  Default is the number of explanatory variables, including \emph{WXB}.
This option helps to avoid singularity since the instrumental variables are constructed by a separate regression for each block.}
 \item{maxblock}{Groups with more than \emph{maxblock} observations are omitted.  Unlimited by default.  This option may be useful for very large
data sets as full nblock x nblock matrices must be constructed for each block, where nblock is the number of observations in the block.}

  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
  \item{silent }{If \emph{silent=T}, no output is printed}
  \item{minp}{Specifies a limit for the estimated probability.  Any estimated probability lower than \emph{minp} will be set to \emph{minp} and
any probability higher than 1-\emph{minp} will be set to 1-\emph{minp}.  By default, the estimated probabilities are bounded by 0 and 1.}
}


\value{
  \item{coef}{Coefficient estimates.}
  \item{se}{Standard error estimates.}
  \item{u}{The generalized error term.}
  \item{gmat}{The matrix of gradient terms, G.}
}
\details{
The linearized model is a three-step estimation procedure.  Let \emph{y} be the indicator value:  \emph{y} = 1 when \eqn{y^* > 0}{y* > 0} and \emph{y} = 0 
when \eqn{y^* < 0}{y* < 0}.  
The first stage is standard probit of \emph{y} on \emph{X}.  The probability estimates from this regression are \eqn{p = \Phi(X \hat{\beta)}}{p = \Phi(X \beta)} 
and the generalized error is \eqn{e = (y-p)\phi(X \hat{\beta})/(p(1-p))}{e = (y-p)*\phi(X \beta)/(p(1-p))}.      
The second/third stage of the procedure is standard 2SLS estimation of \eqn{u = e + gX \hat{\beta}}{u = e + gX \beta} on 
  \eqn{gX} and \eqn{gWX \hat{\beta}}{gWX \beta} using \emph{Z }as instruments, where \emph{g} is the gradient vector, \eqn{-de/d \hat{\beta}}{-de/d \beta}.  
The covariance matrix (equation 3 in Klier-McMillen, 2008) is estimated using the \emph{car} package.  
The final estimates minimize \eqn{e'Z(Z'Z)^{-1}Z'e} with \emph{e} linearized around \eqn{\hat{\beta}}{\beta-probit} and \emph{p} = 0.     


\emph{spprobit} provides flexibility in specifying the list of instruments.  
By default, the instrument list includes \emph{X} and \emph{WX}, where \emph{X} is the original explanatory variable list and \emph{W} is the spatial weight matrix. 
Either \emph{wmat} or \emph{shpfile} must be specified if \emph{inst} and \emph{winst} are set to their default values.


It is also possible to directly specify the full instrument list or to include only a subset of the \emph{X} variables 
in the list that is to be pre-multiplied by \emph{W}.  Let \emph{list1} and \emph{list2} be user-provided lists of the form \emph{list=~z1+z2}.  
The combinations of defaults (\emph{NULL}) and lists for \emph{inst} produce the following results for \emph{Z}: 



1. \emph{inst = NULL},  \emph{winst = NULL}, and either \emph{shpfile} or \emph{wmat} specified:   \emph{Z = (X, WX)}

2. \emph{inst = list1}, \emph{winst = NULL}, and either \emph{shpfile} or \emph{wmat} specified:   \emph{Z = list1} 

3. \emph{inst = NULL},  \emph{winst = list2}, and either \emph{shpfile} or \emph{wmat} specified:  \emph{Z = (X, W*list2)} 

4. \emph{inst = list1}, \emph{winst = list2}, and either \emph{shpfile} or \emph{wmat} specified:  \emph{Z = (list1, W*list2)} 

5. \emph{inst = list1}, \emph{winst = list2}, and both \emph{shpfile} and \emph{wmat} NOT specified:  \emph{Z = (list1, list2)} 

Note that when \emph{inst=list1} and \emph{winst=NULL} it is up to the user to specify at least one variable in \emph{list1} that is not also included in \emph{X}.

The difference between cases (4) and (5) is that the \emph{list2} variables are left unaltered in case (5) rather than being pre-multiplied by \emph{W}.
The case (5) option makes it possible to avoid manipulations of large matrices from within \emph{spprobit}.  The idea is that 
W*list2 should be calculated prior to running \emph{spprobit}, with the variables implied by W*list2 being provided directly to \emph{spprobit} using 
the \emph{winst} option.

 }
\references{
Klier, Thomas and Daniel P. McMillen, "Clustering of Auto Supplier Plants in the United States:  Generalized Method of Moments Spatial Logit for Large Samples," \emph{Journal of 
Business and Economic Statistics} 26 (2008), 460-471.
 }

\examples{
set.seed(9947)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
wmat <- makew(cmap)$wmat
n = nrow(wmat)
rho = .4
x <- runif(n,0,10)
ystar <- as.numeric(solve(diag(n) - rho*wmat)\%*\%(x + rnorm(n,0,2)))
y <- ystar>quantile(ystar,.4)
fit <- spprobit(y~x,  wmat=wmat)
}


\seealso{
  \link{cparlogit}

  \link{cparprobit}

  \link{cparmlogit}

  \link{gmmlogit}

  \link{gmmprobit}

  \link{splogit}

  \link{spprobitml}

}

\keyword{Spatial AR Model}
\keyword{Parametric Models}
\keyword{Discrete Choice Models}
