\name{makew}
\alias{makew}
\title{Calculation of spatial weight matrices}
\description{Constructs a spatial weight matrix from a shape file or a matrix of geographic coordinates}
\usage{
makew(shpfile=NULL,coormat=NULL,method="queen",knum=10,
  ringdist=.25,kern="tcub",window=.10,eigenvalues=FALSE)
}
\arguments{
  \item{shpfile}{A shape file. }
  \item{coormat}{A matrix of geographic coordinates.  The first column should be the longitude and the second latitude, in degrees.}
  \item{method}{Options using shape files to identify first-order contiguity are "queen" and "rook".
Options requiring a matrix of geographic coordinates are "knear", "ring", and "kernel".  
The coordinate matrix can be inputted directly with the \emph{coormat} option, or it can be calculated directly from the shape file.
The shape file takes precedence if both \emph{shapefile} and \emph{coormat} are specified.
 }
  \item{knum}{The number of nearest neighbors for the \emph{knear} option. Default:  \emph{knum} = 10 }
  \item{ringdist}{The maximum distance for the \emph{ring} option. Default: \emph{ringdist} = .25 }
  \item{kern}{The kernel function for the \emph{kernel} option.  Options include \emph{"rect", "tria", "epan", "bisq", "tcub"} and \emph{"trwt"}.}
  \item{window}{Window size for the \emph{kernel} option.  Default:  \emph{window} = .25 }
  \item{eigenvalues}{If \emph{TRUE}, calculates eigenvalues.  Default:  \emph{eigenvalues=FALSE}.}
}

\value{
  \item{wmat}{The nxn spatial weight matrix.  The matrix is row-normalized.}
  \item{eigvar}{The eigvenvalues of \emph{wmat}.  Calculated if \emph{eigenvalues=TRUE}. }
}


\details{
If \emph{method=rook} or \emph{method=queen}, an nxn contiguity matrix is defined using the \emph{spdep} package.
A \emph{queen} definition of contiguity means that two tracts are defined as contiguous if they share at least one point. 
A better definition of \emph{rook} might be "non-queen":  two tracts are defined as contiguous if they share two or more points. 
The \emph{rook} and \emph{queen} options require a shape file.  The contiguity matrix is row-normalized to form the weight matrix.

If \emph{method=knear}, the k nearest neighbors are each given a weight of 1/k to form W.  The calculations are made using the \emph{spdep} package.
Either a shape file or a matrix of geographic coordinates can be provided to the \emph{knear} option.

If \emph{method=ring}, each observation within a distance of \emph{ringdist} from observation i is given equal weight in row i.
More distant observations receive a weight of zero.

If \emph{method=kernel}, a kernel weight function is used to define W, with the window size determined by the \emph{window} option.
The kernel weight function is defined by the \emph{kern} option. 
The weights are \eqn{W_{ij} = K(d_{ij}/h)/h} for \eqn{d_{ij}<h} and \eqn{W_{ij} = 0} for i = j and \eqn{d_{ij}>h}.  The matrix is then row-normalized. 

Eigenvalues are returned if the \emph{eigenvalues=T} is specified.

}


\examples{
library(maptools)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$POPULATION>0&cmap$AREA>0,]
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
lmat <- coordinates(cmap)
fit <- makew(shpfile=cmap,method="queen")
# fit <- makew(coormat=lmat,method="ring",ringdist=1)


}
\seealso{
  \code{\link{sarml}}

  \code{\link{qregspiv}}

}

\keyword{Contiguity Matrix}
\keyword{Spatial AR Model}

