\name{subnp}
\alias{subnp}
\title{
Nonparametric Employment Subcenter Identification
}
\description{
Subcenter identification based on significant residuals from \emph{LWR} employment density estimates
}
\usage{
subnp(ydens, long, lat, window = 0.5, pval = 0.10)
}
\arguments{
  \item{ydens}{The dependent variable. Generally employment density or its natural log.}
  \item{long}{Longitude variable, in degrees}
  \item{lat}{Latitude variable, in degrees}
  \item{window}{The window size for the nonparametric regression.  Default:  .50}
  \item{pval}{Significance level for identifying significant residuals.  Default:  .10, corresponding to a \emph{t}-value of 1.64.}
}

\value{
\item{subobs }{\emph{subobs} = 1 if observation is identified as being part of a subcenter.}
}
\references{
McMillen, Daniel P., "Nonparametric Employment Subcenter Identification," \emph{Journal of Urban Economics} 50 (2001), 448-473.

McMillen, Daniel P., "Employment Densities, Spatial Autocorrelation, and Subcenters in Large Metropolitan Areas," \emph{Journal of Regional Science} 44 (2004), 225-243.

McMillen, Daniel P. and Stefani C. Smith, "The Number of Subcenters in Large Urban Areas," \emph{Journal of Urban Economics} 53 (May 2003), 321-338. 

}


\examples{

library(maptools)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap$acres <- cmap$AREA*640
lmat <- coordinates(cmap)
fit <- geodistance(longvar=lmat[,1], latvar=lmat[,2],  
  lotarget=-87.627800, latarget=41.881998)
cmap$dcbd <- fit$dist
dnorth <- fit$dnorth
deast <- fit$deast

cmap$empdens <- ifelse(cmap$dcbd<=10, exp(log(50) - .16*cmap$dcbd), 10)
dsub1 <- sqrt((dnorth-10)^2 + (deast+10)^2)
dsub2 <- sqrt((dnorth+15)^2 + (deast+5)^2)
cmap$empdens <- cmap$empdens + ifelse(dsub1<=3, exp(log(20) - .19*dsub1), 0)
cmap$empdens <- cmap$empdens + ifelse(dsub2<=3, exp(log(20) - .19*dsub2), 0)
cmap$emp <- round(cmap$empdens*cmap$acres)
fit <- subnp(cmap$empdens,lmat[,1],lmat[,2])
cmap$subobs <- fit$subobs

mapplot(cmap,"subobs")
}
\seealso{
  \code{\link{subgs}}
}

\keyword{Subcenter Identification}


