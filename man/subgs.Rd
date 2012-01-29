\name{subgs}
\alias{subgs}
\title{
Giuliano-Small Employment Subcenter Identification Procedure
}
\description{
Uses spatial weight matrices to find the set of contiguous tracts in a shape file that each 
have at least \emph{mind} employees and that together have total employment of at least \emph{totemp}.

}
\usage{
subgs(shpfile, dens, emp, mind = 10, totemp = 10000, wmat = 0)
}
\arguments{
  \item{shpfile}{A shape file to be used to calculate the spatial contiguity matrix.}
  \item{dens}{The employment density variable}
  \item{emp}{The total employment variable}
  \item{mind}{The minimum density threshold for employment density.  Default:  \emph{mind}=10}
  \item{totemp}{The threshold for total employment across the set of contiguous tracts.  Default:  \emph{totemp}=10000.
}
  \item{wmat}{The spatial contiguity matrix.  Can be inputted directly to program.  
By default, \emph{wmat}=0, leading the program to calculate \emph{wmat} from the shape file specified using the \emph{shpfile} option.
}
}

\value{
\item{subemp}{Total employment in subcenter, 0 if not part of a subcenter.  
Can be used to identify individual subcenters unless two subcenters happen to have identical total employment.
}
\item{subobs}{\emph{subobs} = 1 if observation is identified as being part of a subcenter.}
}
\references{
Giuliano, Genevieve and Kenneth A. Small, Subcenters in the Los Angeles Region, \emph{Regional Science and Urban Economics} 21 (1991), 163-182. 

McMillen, Daniel P., "Identifying Subcenters using Contiguity Matrices," \emph{Urban Studies} 40 (2003), 57-69.

McMillen, Daniel P. and T. William Lester, "Evolving Subcenters:  Employment and Population Densities in Chicago, 1970-2020," \emph{Journal of Housing Economics} 12 (2003), 60-81. 

McMillen, Daniel P. and John F. McDonald, "A Nonparametric Analysis of Employment Density in a Polycentric City," \emph{Journal of Regional Science} 37 (Nov. 1997), 591-612.

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
fit <- subgs(cmap,cmap$empdens,cmap$emp,mind=20,totemp=10000)
cmap$subobs <- fit$subobs
cmap$subemp <- fit$subemp

mapplot(cmap,"subemp")


}

\seealso{
  \code{\link{subnp}}
}

\keyword{Subcenter Identification}
\keyword{Giuliano-Small}

