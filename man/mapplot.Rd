\name{mapplot}
\alias{mapplot}
\title{Maps for Shape File Data}
\description{A simple interface to the SP package's \emph{spplot} command}
\usage{
mapplot(shpfile,varname,brks=8,col="RdBu",sampvar=NULL,print=TRUE,title=NULL) 
}
\arguments{
  \item{shpfile}{The shape file}
  \item{varname}{The name of the variable to be plotted. Must be in quotation marks.}
  \item{brks}{Either the number of break points or a vector of break points for the intervals shown on the map.}
  \item{col}{Color scheme for the map.  All of \emph{RColorBrewer's} schemes are available, e.g., "Blues", "BuGn", "BuPu", and "GnBu".
Default:  "RdBu"}
  \item{sampvar}{A variable indiciating which observations to include in the map.  Tracts for which \emph{sampvar=F} will not be drawn.}
  \item{print}{If \emph{print=T}, prints the endpoints for the intervals shown on the map.  Default:  \emph{print=T}.}
  \item{title}{A title for the map.}
}


\examples{
library(maptools)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",package="McSpatial"))
cmap$lnpopdens <- log(cmap$POPULATION/cmap$AREA)
mapplot(cmap,"lnpopdens",sampvar=cmap$AREA>0&cmap$POPULATION>0,col="Blues",
  title="Log of Population Density")
}
\seealso{
  \link{mapoverlay} 

  \link{mappoints}

}
\keyword{Maps}
