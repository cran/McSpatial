\name{geodensity}
\alias{geodensity}
\title{Kernel density functions for geo-coded data }
\description{Calculates kernel density functions for geo-coded data based on straight-line distances between observations }
\usage{
geodensity(longvar,latvar,window=.25,kern="tcub",alldata=FALSE) 
}
\arguments{
  \item{longvar}{Longitude variable, in degrees.}
  \item{latvar}{Latitude variable, in degrees.}
  \item{window }{Window size.   Default:  0.25.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", and "trwt".}
  \item{alldata }{If \emph{alldata=T}, each observation is used as a target value for \emph{x}.  
When \emph{alldata=F}, the function is estimated at a set of points chosen by the \emph{locfit} program using an adaptive decision tree approach, 
and the \emph{smooth12} command is used to interpolate to the full set of observations.  
Specifying \emph{alldata=T} can lead to long estimation times.}
}

\value{

  \item{target }{The matrix of target values.  Dimensions = ntx2, where nt is the number of target points.
  First column = longitude, second column = latitude.}
  \item{dens.target }{The estimated densities at the target coordinates.}
  \item{denshat }{The estimated densities at the original data points.}
}

\seealso{
  \link{geodistance}

  \link{geogravity}
}


\examples{
\dontrun{
library(spdep)
library(RColorBrewer)
cook <- readShapePoly(system.file("maps/CookCensusTracts.shp",package="McSpatial"))
# measure distance to Chicago city center
lmat <- coordinates(cook)
cook$longitude <- lmat[,1]
cook$latitude  <- lmat[,2]
fit <- geodensity(cook$longitude, cook$latitude)
cook$denshat <- fit$denshat
brks <- seq(min(cook$denshat,na.rm=TRUE),max(cook$denshat,na.rm=TRUE),length=9)
spplot(cook,"denshat",at=brks,col.regions=rev(brewer.pal(9,"RdBu")),
  main="Census Tract Densities")
}
}




\details{
  The \emph{geodistance} function is used to calculate straight-line distances between all observations and each target point.
The vector of distances for a given target value is \emph{d}.  The window is determined by finding \eqn{dmax = quantile(d, window)}.
The estimated density at the target point is simply:

\deqn{f = \frac{1}{dmax*n} \sum_i K ( \frac{d_i}{dmax} ) }{f = (1/(dmax*n)) \sum_i K ( d_i/dmax ) }

Available kernel weighting functions include the following:

\tabular{lll}{
Kernel \tab Call abbreviation \tab Kernel function K(z) \cr
Rectangular \tab ``rect'' \tab \eqn{\frac{1}{2} I(|z| <1)}{1/2 * I(|z|<1)} \cr
Triangular \tab ``tria'' \tab \eqn{(1-|z|)I(|z|<1)}{(1-|z|) * I(|z|<1)}\cr
Epanechnikov \tab ``epan'' \tab \eqn{\frac{3}{4} (1-z^2) * I(|z| <1)}{3/4 * (1-z^2)*I(|z| < 1)} \cr
Bi-Square \tab ``bisq'' \tab \eqn{\frac{15}{16} (1-z^2)^2 * I(|z| <1)}{15/16 * (1-z^2)^2 * I(|z| < 1)} \cr
Tri-Cube \tab ``tcub'' \tab \eqn{\frac{70}{81} (1-|z|^3)^3 * I(|z| <1)}{70/81 * (1-|z|^3)^3 * I(|z| < 1)} \cr
Tri-Weight \tab ``trwt'' \tab \eqn{\frac{35}{32} (1-z^2)^3 * I(|z| <1)}{35/32 * (1-z^2)^3 * I(|z| < 1)} \cr
}

The gaussian kernel is not available.

If \emph{alldata=T}, each data point in turn is used as a target point.
If \emph{alldata=F}, \emph{locfit} is used to find a set of target points, and the \emph{smooth12} command is used to interpolate to the full set of observations.
The matrix of target coordinates is stored in \emph{target}, and the estimated densities at the target points are stored in \emph{dens.target}.  
If \emph{alldata=T}, \emph{target} contains the full set of values for longitude and latitude, and \emph{dens.target} = \emph{denshat}.

}


\keyword{Density Functions}
\keyword{Nonparametric}

