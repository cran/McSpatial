\name{cookdata}
\alias{cookdata}
\docType{data}
\title{Data set associated with CookCensusTracts shape file}
\description{Census tract data for Cook County, IL for 2000}
\usage{data(cookdata)}
\format{
  A data frame with 1343 observations on the following variables:
  \describe{
    \item{\code{AREA}}{Area of the census tract in square miles}
    \item{\code{TRACT}}{Census tract number}
    \item{\code{POPULATION}}{Total population}
    \item{\code{HHMEDINC}}{Median household income}
    \item{\code{FAR}}{Average floor-area ratio for single-family residential homes}
    \item{\code{AGE}}{Average age of single-family residential homes}
    \item{\code{CHICAGO}}{Indicates whether the tract is in Chicago. 
Note that the "O'Hare" community area is not counted as part of Chicago when defining this variable.}
    \item{\code{CAREA}}{Community area name for tracts in Chicago}
    \item{\code{LONGITUDE}}{Longitude of the tract centroid}
    \item{\code{LATITUDE}}{Latitude of the tract centroid}
    \item{\code{DCBD}}{Distance from the Chicago CBD, measured in miles}
    \item{\code{LNFAR}}{Natural log of FAR}
    \item{\code{LNDENS}}{Natural Log of Dens}
  }
}

\keyword{Datasets}
\seealso{
  \link{CookCensusTracts}
}
