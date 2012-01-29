\name{mapoverlay}
\alias{mapoverlay}
\title{Merges variables across shape and point files}
\description{Associates values for a shape file variable with the observations in a point file. Also can match means, sums, or other functions of 
point file observations to a shape file.}
\usage{
mapoverlay(shpfile,pointfile,shpvar,pointvar=NULL,func="sum")
}
\arguments{
  \item{shpfile}{A shape file.}
  \item{pointfile}{A data frame with point data.  Must include the geographic coordinates for the point data.  
The geographic coordinates must be stored in  variables named \emph{longitude} and \emph{latitude}.}
  \item{shpvar}{The name of the variable in the shape file whose values are to be added to the point data.  Must be in quotes.
Example:  \emph{shpvar}="TRACT" or \emph{shpvar}="city".}
  \item{pointvar}{Either NULL or the name of the variable in \emph{pointfile}.  In the latter case, the name must be in quotes, e.g., \emph{pointvar}="price".}
  \item{func}{The function to be applied to \emph{pointvar}, such as "mean" or "sd".  If \emph{pointvar}=NULL, \emph{func} is ignored.  
Default:  \emph{func} = "sum".  }
}
\value{
  \item{pointout}{A variable with the values of \emph{shpvar}, matched to the observations in \emph{pointfile}.  
Stored in the same order as the \emph{pointfile} observations, so it can be added directly to \emph{pointfile}. }
  \item{shpout}{A variable with the values of \emph{pointvar}, matched to the observations in \emph{shpfile}.  
Stored in the same order as the observations in \emph{shpfile}.}
}

\details{
  The \emph{mapoverlay} command serves as a simple interface to the \emph{sp} package's \emph{overlay} command.  
As an example, \emph{shpfile} might represent census tracts while \emph{pointfile} includes data on home sales, including the geographic coordinates.
The geographic coordinates for \emph{pointfile} must be named \emph{longitude} and \emph{latitude},  
\emph{shpfile} must be a polygon file, and \emph{pointfile} must be a data frame.

In this example, suppose that the name of the census tract variable in \emph{shpfile} is "TRACT".  Also, assume there are n1 observation in
\emph{shpfile} and n2 observations in \emph{pointfile}.
Then the output from \emph{fit <- mapoverlay(shpfile,pointfile,shpvar="TRACT")} is a vector with n2 entries indicating the tract in which each
\emph{pointfile} observation is located.  The variable can be added to \emph{pointfile} data set using the command
\emph{pointfile$TRACT <- fit$pointout}. 

While the \emph{shpvar} option is required, \emph{pointvar} is only needed to aggregate a \emph{pointfile} variable to the shape file.
For example, suppose we want to calculate the average home sale price in each census tract.  The command
\emph{fit <- mapoverlay(shpfile, pointfile, shpvar="TRACT", pointvar="price", func="mean")} calculates the mean of \emph{pointfile$price} for each
value of \emph{shpfile$TRACT}.  In this case, \emph{fit$pointout} is the same as before, but \emph{fit$shpout} is a variable with n1 entries 
with the average values of \emph{pointfile$price}.  It can be added directly to the shape file with the command
\emph{shpfile$price <- fit$shpout}.

}


\examples{
library(maptools)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
data(matchdata)
matchdata$price <- exp(matchdata$lnprice)
fit <- mapoverlay(cmap,matchdata,shpvar="TRACT",
  pointvar="price",func="mean")
matchdata$TRACT <- fit$pointout
table(matchdata$TRACT)
cmap$avgprice <- fit$shpout
psamp <- !is.na(cmap$avgprice)
mapplot(cmap,"avgprice",sampvar=psamp)

}

\seealso{
  \link{mapplot}

  \link{mappoints}
}

\keyword{Maps}
