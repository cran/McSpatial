\name{mappoints}
\alias{mappoints}
\title{Maps showing points on a shape file}
\description{Uses the \emph{classInt} package to show a color ramp of points placed on a shape file map}
\usage{
mappoints(shpfile,longvar,latvar,pointvar,sampvar=NULL,nclass=8,
  col="RdBu",pointsize=.6,colreverse=TRUE,xlimfact=c(1,1),ylimfact=c(1,1),
  legendloc=NULL,legendsize=.75,legenddigits=4,title=NULL)
}
\arguments{
  \item{shpfile}{The shape file, a polygon map showing areas.}
  \item{longvar}{Longitude variable for the point data, in degrees.}
  \item{latvar}{Latitude variable for the point data, in degrees.}
  \item{pointvar}{The variable to be mapped.}
  \item{sampvar}{A variable indiciating which observations to include in the map.  Tracts for which \emph{sampvar=F} will not be drawn.}
  \item{nclass}{The number of classes used when mapping \emph{pointvar}.  Default:  \emph{nclass=8}.}
  \item{col}{Color scheme for the map.  All of \emph{RColorBrewer's} schemes are available, e.g., "Blues", "BuGn", "BuPu", and "GnBu".
Default:  "RdBu"}
  \item{pointsize}{The size of the points drawn on the map.  Fed directly to the \emph{cex} option for the \emph{points} command.  Default:  \emph{pointsize=.6}.}
  \item{colreverse}{If \emph{TRUE}, reverses the order of the color ramp for the \emph{col} options.  Default:  \emph{colreverse=T}.}
  \item{xlimfact}{A vector of factors used to adjust the limits for the box around the map.  The limits for the x-axis (i.e., the longitudes) are 
xlim=c(min(longvar),max(longvar))*xlimfact.  Since the geographic coordinates are expressed in degrees, small adjustments will produce large
changes in the box.  Values like \emph{xlimfact} = c(.999,1.001) can produce large changes.  Default:  \emph{xlimfact} = c(1,1).}
  \item{ylimfact}{The counterpart to \emph{xlimfact} for the \emph{y} (i.e., latitude) axis.  Default:  \emph{ylimfact} = c(1,1).}
  \item{legendloc}{The location for the legend, e.g., "topleft" or "bottom".  By default, the legend is omitted.}
  \item{legendsize}{The font size for the legend.  Fed directly to the \emph{cex} option for the \emph{legend} command.  Default:  \emph{legendsize=.75}.}
  \item{legenddigits}{Number of significant digits for the numbers shown in the legend.  Default:  \emph{legenddigits} = 4.}
  \item{title}{A title for the map.}
}

\examples{
library(maptools)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
carea <- cmap$CAREA=="Albany Park"|cmap$CAREA=="Edgewater"|
  cmap$CAREA=="Edison Park"|cmap$CAREA=="Forest Glen"|
  cmap$CAREA=="Jefferson Park"|cmap$CAREA=="Lincoln Square"|
  cmap$CAREA=="North Park"|cmap$CAREA=="Norwood Park"|
  cmap$CAREA=="Rogers Park"|cmap$CAREA=="Uptown"|
  cmap$CAREA=="West Ridge"
carea[is.na(carea)] <- FALSE

data(matchdata)
mappoints(cmap,matchdata$longitude,matchdata$latitude,matchdata$lnprice,
  sampvar=carea,legendloc="bottomleft",title="Location of Home Sales")

}
\seealso{
  \link{mapoverlay} 

  \link{mapplot}
}
\keyword{Maps}
