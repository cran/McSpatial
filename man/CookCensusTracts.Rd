\name{CookCensusTracts}
\alias{CookCensusTracts}
\docType{data}
\title{
Shapefile of Census Tracts in Cook County for 2000
}
\description{
A map of census tracts in Cook County, Illinois with data on population, median household income, average floor area ratios, and average home ages.
}

\format{
  A shape file with 1343 census tracts.  Can be read directly into a GIS program.
  \describe{
    \item{\code{SP_ID}}{An identification variable}
    \item{\code{AREA}}{Area of the tract in square miles}
    \item{\code{TRACT}}{Census tract number, including state and county codes}
    \item{\code{POPULATION}}{Population in 2000}
    \item{\code{HHMEDINC}}{Median household income in the census tract in 2000}
    \item{\code{FAR}}{Average floor area ratio for census tract in 2000.  Calculated from Cook County assessment file}
    \item{\code{AGE}}{Average age of homes in 2000.  Calculated from Cook County assessment file.}
    \item{\code{CHICAGO}}{Indicates whether the tract is part of the City of Chicago.  O'Hare airport is NOT included in Chicago. 
It is identified in the variable \emph{CAREA}.}
    \item{\code{CAREA}}{For Chicago observations, indicates the community area.  Missing for suburban observations. }
  }
}

\source{
Daniel McMillen.  Data are drawn from the U.S. Census and the Cook County Assessment File.
}
\examples{
library(maptools)
cook <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
# to add O'Hare back in:
sampvar <- cook$CHICAGO==1|(!is.na(cook$CAREA)&cook$CAREA=="O'Hare")
cook <- cook[sampvar==TRUE,]

}
\keyword{Datasets}
\keyword{Maps}

