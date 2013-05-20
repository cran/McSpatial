\name{ksim }
\alias{ksim }
\title{
Estimates K-density functions with local or global confidence intervals for counter-factual locations
 }
\description{
Calculates \emph{K}-density functions for lat-long coordinates.  
Calculates the distance, \emph{d}, between every pair of observations and plots the density, \eqn{f(d_0)} at a set of target distances, \eqn{d_0}.  
Also uses the Duranton-Overman bootstrap method to construct local or global confidence intervals for the density of distances 
between pairs of observations if the same number of points were allocated across another set of possible locations. 
}

\usage{
ksim(long1,lat1,long2,lat2,kilometer=FALSE,noplot=FALSE,
  dmin=0,dmax=0,dlength=512,h=0,kern="gaussian",
  nsim=2000,nsamp=0,pval=.05,cglobal=FALSE) 
}

\arguments{
  \item{long1 }{Longitude variable, in degrees. }
  \item{lat1 }{Latitude variable, in degrees. }
  \item{long2 }{Longitude variable for counter-factual locations, in degrees. }
  \item{lat2 }{Latitude variable for counter-factual locations, in degrees.}
  \item{kilometer }{If \emph{kilometer = T}, measurements are in kilometers rather than miles.  Default:  \emph{kilometer = F}. }
  \item{noplot }{If \emph{noplot = T}, does not show the \emph{K}-density graph. }
  \item{dmin }{Minimum value for target distances.  Default:  \emph{dmin}=0. }
  \item{dmax }{Maximum value for target distances.  Default:  \emph{dmin} = max(distance). }
  \item{dlength }{Number of target values for density calculations.  Default:  \emph{dlength} = 512. }
  \item{h }{Bandwidth.  Default:  
(.9*(quantile(distance,.75)-quantile(distance,.25))/1.34)*(n^(-.20)), where n = 2*length(dvect).
 }
  \item{kern }{Kernel.  Default:  "gaussian".  
Other options from the \emph{density} function are also available, 
including "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine".
 }
  \item{nsim}{Number of simulations for constructing the confidence intervals.  Default:  \emph{nsim}=2000. }
  \item{nsamp }{If \emph{nsamp}>0, uses a random sample of lat-long pairs for calculations rather than full data set.  
Takes random draws from \emph{long1, lat1} pairs; the \emph{long2, lat2} remain as specified by the user.
Can be much faster for large samples.  Default:  use full sample. }
  \item{pval}{Significance level for confidence intervals.  Default: \emph{pval} = .05, i.e., a 95 percent confidence interval. }
  \item{cglobal}{If \emph{cglobal=T}, calculates global confidence intervals.  Default:  \emph{cglobal=F}, calculates local confidence intervals. }
}


\value{
  \item{distance }{The vector of target distances }
  \item{dhat }{The vector of densities for the target distances }
  \item{h }{The bandwidth }
  \item{local.lo }{The local confidence interval at each target distance, if calculated. }
  \item{local.hi }{The local confidence interval at each target distance, if calculated. }
  \item{global.lo }{The global confidence interval at each target distance, if calculated. }
  \item{global.hi }{The global confidence interval at each target distance, if calculated. }

}
\details{
Let \emph{n} be the number of observations in the \emph{long1, lat1} data set.  
\emph{ksim} draws \emph{n} observations from the \emph{long2, lat2} pairs and re-calculates the \emph{K}-density function using the new, simulated data set.  
The process is repeated \emph{nsim} times, producing \emph{nsim} bootstrap \emph{K}-density functions.  
The local confidence interval treats each target distance as a separate observation, 
and calculates the confidence interval at each distance using the standard bootstrap percentile method.  
In contrast, the global confidence interval treats the full \emph{K}-density function as an observations and shifts the interval outward at each data point until 95 percent of 
the density functions lie within the interval.  Large values of \emph{nsim} - perhaps greater than the default of 2000 - 
are necessary to get accurate global confidence intervals. 


The \emph{ksim} function is intended for cases where the counterfactual data set has more observations than the base, i.e., n2>n1. 
In this case, observations are drawn without replacement from the counterfactual data set.  
When the counterfactual data set has fewer observations than the base (i.e., n2<=n1), n1 observations are drawn with replacement from the counterfactual data set.
 

Duranton and Overman (2005) proposed this method for constructing global confidence intervals for \emph{K}-density functions.  
See Klier and McMillen (2008) for a description of the procedures used here.  
See the description of the \emph{kdensity} function for more details on the estimation procedure of the \emph{K}-density functions.

 }

\references{
Duranton, Gilles and Henry G. Overman, "Testing for Localisation using Microgeographic Data", \emph{Review of Economic Studies} 72 (2005), 1077-1106. 

Klier Thomas and Daniel P. McMillen, "Evolving Agglomeration in the U.S. Auto Industry," \emph{Journal of Regional Science} 48 (2008), 245-267. 

Silverman, A. W., \emph{Density Estimation for Statistics and Data Analysis}, Chapman and Hall, New York (1986). 

 }

\seealso{ \link{kdensity} }

\examples{
data(matchdata)
lmat <- cbind(matchdata$longitude,matchdata$latitude)
lmat1 <- lmat[matchdata$carea=="Rogers Park"|matchdata$carea=="Albany Park",]
lmat2 <- lmat[matchdata$carea!="Rogers Park"&matchdata$carea!="Albany Park",]
# smaller samples to reduce time for examples
set.seed(4941)
obs <- sample(seq(1,nrow(lmat1)),200)
lmat1 <- lmat1[obs,]
obs <- sample(seq(1,nrow(lmat2)),400)
lmat2 <- lmat2[obs,]

fit <- ksim(lmat1[,1],lmat1[,2],lmat2[,1],lmat2[,2],dmax=9,nsim=100,
  nsamp=100,noplot=TRUE,cglobal=FALSE)
ymin = min(fit$dhat,fit$local.lo)
ymax = max(fit$dhat,fit$local.hi)
plot(fit$distance, fit$dhat, xlab="Distance", ylab="Density", ylim = c(ymin,ymax), 
  type="l", main="Albany Park & Rogers Park v. Other Areas")
lines(fit$distance, fit$local.lo, col="red")
lines(fit$distance, fit$local.hi, col="red")
}

\keyword{Density Functions}
\keyword{Nonparametric}

