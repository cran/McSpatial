\name{geogravity}
\alias{geogravity}
\title{Gravity matrix and gravity variable calculations}
\description{
Calculates a variable showing the average value for each observation of a gravity measure of the spatial interaction between the observation 
and the other points in a data set. 
}
\usage{
geogravity(x,longvar,latvar,alpha=1,maxd=NULL,alldata=FALSE,
  window=.10,outmatrix=FALSE)
}
\arguments{
  \item{x}{The variable of interest, e.g., population or employment.}
  \item{longvar}{Longitude variable, in degrees.}
  \item{latvar}{Latitude variable, in degrees.}
  \item{alpha}{The \eqn{\alpha} parameter for the distance variables.  Default:  \emph{alpha} = 1.}
  \item{maxd}{Maximum distance, beyond which observations get zero weight. Default:  all observations are included in the calculations.}
  \item{alldata}{If \emph{FALSE}, interpolates between target points rather than making the calculations at every observation.  Default:  \emph{alldata=FALSE}.}
  \item{window}{Window size used to determine a set of target points when \emph{alldata=FALSE}.  Default:  \emph{window} = .10.}
  \item{outmatrix}{If \emph{TRUE} and \emph{alldata==TRUE}, stores the full matrix of gravity values.  Default:  \emph{outmatrix=FALSE}.}
}

\value{
  \item{targetobs}{The vector of target observation numbers.}
  \item{gtarget}{The gravity variable at the target points.}
  \item{gravity}{The gravity variable for the full data set.}
  \item{dmat}{The full nxn matrix of results, if \emph{outmatrix=TRUE} and \emph{alldata=TRUE}.  The diagonal elements are set to zero.}
}


\examples{
library(spdep)
cook <- readShapePoly(system.file("maps/CookCensusTracts.shp",package="McSpatial"))
cook <- cook[cook$POPULATION>0&cook$AREA>0,]
cook$lndens <- log(cook$POPULATION/cook$AREA)
lmat <- coordinates(cook)
longitude <- lmat[,1]
latitude  <- lmat[,2]
fit <- geogravity(cook$lndens,longitude,latitude)
cook$gravity <- fit$gravity

}

\seealso{
  \link{geodistance}

  \link{geodensity}
}

\details{
  The gravity measure of the spatial interaction between two points \emph{i} and \emph{j} is
\eqn{ g_{ij} =  P_i P_j / d_{ij}^\alpha} if \eqn{d_{ij} \leq maxd}{d_{ij} < maxd} and \eqn{g_{ij} = 0} if \eqn{d_{ij} > maxd},
where \emph{P} is a variable such as population or employment that measures the importance of the observation,
\eqn{d_{ij}} is the straight-line distance between observations \emph{i} and \emph{j}, and \eqn{\alpha} is a parameter.
The variable \emph{P} is provided by the argument \emph{x}.
The full \eqn{n \times n}{n x n} matrix of values is stored in \emph{dmat} if \emph{alldata=TRUE} and \emph{outmatrix=TRUE}, with the diagonal elements set to zero.

The "gravity variable" that is stored in \emph{gravity} is the average value for each observation:

\deqn{g_i = \frac{1}{n-1} \sum_{j \neq i} g_{ij} }{g_i = (1/(n-1)) \sum g_{ij}}

By default, the \emph{locfit} program is used to find a set of target values for calculating \eqn{ g_{ij}}.
The full set of observations is used for the \emph{j} index, but a smaller set of observations is used for the target values, \emph{i}.
The abbreviated set of \eqn{g_i} values is then interpolated to the full set of data points using the \emph{smooth12} command.
The vector of target observation numbers is stored in \emph{targetobs}. 

}
