\name{maketarget}
\alias{maketarget}
\title{Target Points for Nonparametric Models}
\description{Identifies target points at which to evaluate nonparametric models with one or two explanatory variables.}
\usage{
  maketarget(form,window=.25,bandwidth=0,kern="tcub",actualobs=FALSE,data=NULL)
}

\arguments{
  \item{form}{The model formula to be used to determine target points.  The dependent variable is irrelevant, so form=y~x and form=~x are equivalent.  
No more than two explanatory variables should be listed.}
  \item{window}{The window size.}
  \item{bandwidth}{The bandwidth.  The \emph{window} option is ignored if \emph{bandwidth} >0.}
  \item{kern}{The kernel weight function.}
  \item{actualobs}{If FALSE, the points identifed as targets are not constrained to be actual data points.  
Specifying \emph{actualobs = TRUE} produces target points that are drawn from the original data matrix implied by \emph{form}.}
  \item{data}{A data frame containing the data. Default:  use data in the current working directory}
}

\value{
  \item{target}{The set of target points.}
  \item{obs}{If \emph{actualobs = TRUE}, \emph{obs} is the list of observation numbers in the original data set from which the target points are drawn.
\emph{obs = NULL} if \emph{actualobs = FALSE}.}
}



\details{
The \emph{maketarget} function uses the \emph{locfit} package's adaptive decision tree approach to identify target locations.
If \emph{actualobs = TRUE}, the output of \emph{maketarget} is set of actual data points closest to these target locations, along with the 
convex hull identified by the \emph{chull} command.  The variable \emph{obs} can be used to indicate the target observations for \emph{actualobs} models.
This \emph{actualobs = TRUE} option is required by the following commands:    \link{lwr},  \link{lwrgrid},   \link{cparlwr},   \link{cparlwrgrid},
and   \link{semip}.
}

\references{
Loader, Clive.  \emph{Local Regression and Likelihood.}  New York:  Springer, 1999.  Section 12.2.
}

\examples{
data(cookdata)
target <- maketarget(~LONGITUDE+LATITUDE,window=.25,data=cookdata)$target
}

\seealso{
  \link{cparlogit}
  \link{cparlwr}
  \link{cparlwrgrid}
  \link{cparmlogit}
  \link{cparprobit}
  \link{lwr}
  \link{lwrgrid}
  \link{qregcpar}
  \link{qreglwr}
  \link{semip}
}

