\name{geodistance}
\alias{geodistance}
\title{Calculates distances using the great circle formula}
\description{
\emph{geodistance} calculates distances in miles between a set of observations and a location.  Distances are calculated in miles using the great circle formula.  
Geographic coordinates must be expressed in latitudes and longitudes.
}
\usage{geodistance(longvar, latvar, lotarget, latarget, dcoor = FALSE)}
\arguments{
  \item{longvar}{Longitude variable, in degrees.}
  \item{latvar}{Latitude variable, in degrees.}
  \item{lotarget}{Target longitude.}
  \item{latarget}{Target latitude.}
  \item{dcoor}{If dcoor = T, also calculates the distance east \emph{(deast)} and north \emph{(dnorth)} of the target point.}
}

\value{

  \item{dist}{A vector with the distance in miles between each data point and the target point.}
  \item{dnorth}{A vector with the number of miles north of the target point for each data point. \cr
\emph{dnorth} < 0 for observations that are south of the target point.
}
  \item{deast}{A vector with the number of miles east of the target point for each data point. \cr
\emph{deast} < 0 for observations that are west of the target point.
}
}


\seealso{
  \link{geodensity}

  \link{geoshape}
}


\examples{
data(cookdata)
dcbd <- geodistance(cookdata$LONGITUDE, cookdata$LATITUDE, -87.627800, 41.881998)$dist
}

\keyword{Distance}
\keyword{Great-Circle Distance Formula}

