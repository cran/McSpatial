\name{geoshape}
\alias{geoshape}
\title{Calculates distances between a shape file and a set of points.}
\description{Calculates distances in miles between a shape file and a set of geographic coordinates.}
\usage{
geoshape(longvar,latvar,linefile=NULL,pointfile=NULL,coormatrix=NULL)
}
\arguments{
  \item{longvar}{Longitude variable, in degrees.}
  \item{latvar}{Latitude variable, in degrees.}
  \item{linefile}{A shape file with lines. In this case, \emph{geoshape} calculates the distance between each observation and the nearest line in the shape file.
Default:  \emph{linefile=NULL}.}
  \item{pointfile}{A shape file with points. In this case, \emph{geoshape} calculates the distance between each observation and the nearest point in the shape file.
Also used for polygon files, in which case the centroids are considered points.
Default:  \emph{shapefile=NULL}.}
  \item{coormatrix}{A matrix of geographic coordinates.  The first column must be the longitude and the second column the latitude, both in degrees. 
In this case, \emph{geoshape} calculates the distance between each observation and the nearest point in the matrix.  
The class of the object sent to \emph{coormatrix} must actually be \emph{matrix}.  Default:  \emph{coormatrix=NULL}.}
}

\value{Returns the calculated distances, in miles.}

\details{
  Uses the \emph{nncross} command from the \emph{spatstat} package to calculate distances between the set of points given by \emph{(longvar, latvar)} and the 
shape file provided by one of the \emph{linefile}, \emph{pointfile}, or \emph{coormatrix} options.  
Only one of the three shape file options should be specified, and it is critical that the appropriate option is matched to the shape file.  
A polygon file is considered a point file for purposes of the \emph{geoshape} command; in this case, \emph{nncross} calculates distances to the polygon centroids.  
For the \emph{coormatrix} option, it is critical that the object is of class \emph{matrix},
the first column is longitude, and the second is latitude.
}

\seealso{
  \link{geodistance}
}



\examples{
data(matchdata)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1,]
lmat <- coordinates(cmap)
# Calculate distance between homes in matchdata and the census tract centroids
matchdata$dist1 <- geoshape(matchdata$longitude,matchdata$latitude,pointfile=cmap)
# Alternative method using coormatrix option
matchdata$dist2 <- geoshape(matchdata$longitude,matchdata$latitude,
  coormatrix=coordinates(cmap))

# measure distance from census tract centroids to Chicago city center
longitude <- lmat[,1]
latitude  <- lmat[,2]
cmat <- t(as.matrix(c(-87.627800, 41.881998)))
dcbd <- geoshape(longitude, latitude, coormatrix=cmat)
summary(dcbd)


}

\keyword{Distance}
\keyword{Maps}
