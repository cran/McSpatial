\name{repsaledata}
\alias{repsaledata}
\title{
Preparation of a Repeat Sales Data Set
}
\description{
Identifies repeat sales from a data set with observations on sale price, time of sale, and a property id.  
Returns a data frame in which each observation is a repeat sales pair.
}

\usage{
  repsaledata(price,timevar,id)
}

\arguments{
  \item{price}{Variable representing sale price }
  \item{timevar}{Variable representing date of sale}
  \item{id}{Property id}
}

\value{
  \item{id}{Property id}
  \item{price0}{Sale price at earlier date}
  \item{time0}{Earlier sales date}
  \item{price1}{Sale price at later date}
  \item{time1}{Later sales date}
}

\examples{
id <-      c(1,1,1, 2,2,2, 3,3,3, 4,4,4, 5,5,5)
timevar <- c(1,2,3, 1,2,2, 3,1,1, 1,1,2, 2,2,3)
price <- seq(1:15)
basedata <- data.frame(id,timevar,price)
basedata

rdata <- repsaledata(price,timevar,id)
rdata

}


\details{
The input to \emph{repsaledata} is meant to be a set of variables from a standard hedonic data set -- 
the sale price, date, and an \emph{id} number for the individual property.  The function identifies the subset of properties that sold 
at least twice and forms a new data set in which each observation is a repeat sales pair, with "0" denoting the earlier time and
"1" denoting the later date in the variable names.  The observations are ordered first by id, then by timevar, and then by price.
A repeat sales pair is formed by matching an observation for which id(t)=id(t-1) and timevar(t)=timevar(t-1).  Thus, a property that
sold in times 1, 2, and 3 will produce 2 repeat sales pairs:  (1) t = 1 and 2, and (2) t = 2 and 3.

The output of \emph{repsaledata} is a data frame with 5 variables.  If some of the original hedonic data set variables need to be 
included in the repeat sales data set, the original hedonic data set and the \emph{repsaledata} data frame can be merged by the \emph{id} variable.


}


\seealso{
  \link{repsale}

  \link{repsalefourier}

  \link{repsaleqreg}
}

\keyword{Repeat Sales}
