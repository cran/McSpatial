\name{kdensity }
\alias{kdensity }
\title{K-density functions for distances between geographic coordinates }
\description{
Calculates \emph{K}-density functions for lat-long coordinates.  
Calculates the distance, \emph{d}, between every pair of observations and plots the density, \eqn{f(d_0)}, at a set of target distances, \eqn{d_0}.  
The kernel density functions are calculated using the \emph{density} function. 
}

\usage{
kdensity(longitude,latitude,kilometer=FALSE,noplot=FALSE,
dmin=0,dmax=0,dlength=512,h=0,kern="gaussian",nsamp=0,
confint=TRUE,pval=.05) 
}

\arguments{
  \item{longitude }{Longitude variable, in degrees. }
  \item{latitude }{Latitude variable, in degrees. }
  \item{kilometer }{If \emph{kilometer = T}, measurements are in kilometers rather than miles.  Default:  \emph{kilometer = F}. }
  \item{noplot }{If \emph{noplot = T}, does not show the graph of the \emph{K-density} function. }
  \item{dmin }{Minimum value for target distances.  Default:  \emph{dmin}=0. }
  \item{dmax }{Maximum value for target distances.  Default:  \emph{dmin} = max(distance), specified by setting \emph{dmin}=0. }
  \item{dlength }{Number of target values for density calculations.  Default:  \emph{dlength} = 512. }
  \item{h }{Bandwidth.  Default:  
(.9*(quantile(distance,.75)-quantile(distance,.25))/1.34)*(n^(-.20)), where n = 2*length(dvect).
 }
  \item{kern }{Kernel.  Default:  "gaussian ".  
Other options from the \emph{density} function are also available, 
including "epanechnikov", "rectangular", "triangular", "biweight", and "optcosine".  
The "cosine" kernel is translated to "optcosine".
 }
  \item{nsamp }{If \emph{nsamp}>0, draws a random sample of lat-long pairs for calculations rather than the full data set.  
Can be much faster for large samples.  Default:  use full sample. }
  \item{confint}{If \emph{TRUE}, adds local confidence intervals to the graph.  Default:  \emph{confint=TRUE}.}
  \item{pval}{p-value for confidence intervals.  Default:  \emph{pval=.05}.}
}


\value{
  \item{distance}{The vector of target distances.}
  \item{dhat}{The vector of densities for the target distances.}
  \item{dvect}{The full vector of distances between observation pairs.  Length is n(n-1)/2.}
  \item{h}{The bandwidth.}
  \item{se}{The vector of standard errors.}
}
\details{
The \emph{kdensity} function uses Silverman's (1986) reflection method to impose zero densities at negative densities.  
This method involves supplementing each distance observation with its negative value to form a pseudo data set with twice the original number of observations.  
The following commands are the core of the function:


dfit1 <- density(dvect,from=dmin,to=dmax,n=dlength,kernel=kern,bw=h) \cr
dfit2 <- density(-dvect,from=dmin,to=dmax,n=dlength,kernel=kern,bw=h) \cr
distance <- dfit1$x \cr
dhat <- dfit1$y + dfit2$y \cr


Local standard errors are calculated using the following asymptotic formula:

\eqn{(nh)^{-.5} (f(x) \int K^2(\psi)d \psi )^{.5} }
 }

\references{
Duranton, Gilles and Henry G. Overman, "Testing for Localisation using Microgeographic Data", \emph{Review of Economic Studies} 72 (2005), 1077-1106. 

Klier Thomas and Daniel P. McMillen, "Evolving Agglomeration in the U.S. Auto Industry," \emph{Journal of Regional Science} 48 (2008), 245-267. 

Silverman, A. W., \emph{Density Estimation for Statistics and Data Analysis}, Chapman and Hall, New York (1986). 

 }

\seealso{ \link{ksim}}

\examples{
data(matchdata)
lmat <- cbind(matchdata$longitude,matchdata$latitude)
# Smaller sample to reduce computation time for example
set.seed(18493)
obs <- sample(seq(1,nrow(lmat)),400)
lmat <- lmat[obs,]
fit95 <- kdensity(lmat[,1],lmat[,2],noplot=FALSE)
}


\keyword{Density Functions}
\keyword{Nonparametric}

