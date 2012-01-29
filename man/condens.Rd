\name{condens}
\alias{condens}
\title{
Conditional density estimation
}
\description{Estimates conditional density functions of the form f(y| x) = f(x,y)/f(x).  
Kernel density estimators are used to estimate f(x,y) and f(x).  
The conditional density function can be plotted as a three-dimensional surface or as a contour map.  
Alternatively, the conditional density of y can be graphed for as many as five target values of x.  

}

\usage{

condens(form,window=.7,bandwidth=0,kern="tcub",
  mingrid.x=NULL,maxgrid.x=NULL,mingrid.y=NULL,maxgrid.y=NULL,ngrid=50,
  xlab="x",ylab="y",zlab="fxy/fx",contour=TRUE,level=TRUE,wire=TRUE,dens=TRUE,
  targetx.dens=NULL,quantile.dens=c(.10,.25,.50,.75,.90),data=NULL) 

}

\arguments{
  \item{form}{Model formula}
  \item{window }{Window size.   Default:  0.25. }
  \item{bandwidth }{Bandwidth.   Default: not used. }
  \item{kern}{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{mingrid.x, maxgrid.x, mingrid.y, maxgrid.y,
ngrid}{The mingrid and maxgrid values are the boundaries for the \emph{ngrid x ngrid} lattice used in the graphs produced by \emph{condens}.  
By default, mingrid.x = min(x), maxgrid.x = max(x), mingrid.y = min(y), maxgrid.y = max(y), and ngrid=50.  
}
  \item{xlab}{Label for the \emph{x}-axis in graphs.  Default:  "x"}
  \item{ylab}{Label for the \emph{y}-axis in graphs.  Default:  "y"}
  \item{zlab}{Label for the \emph{z}-axis in graphs.  Default:  "fxy/fx"}
  \item{contour}{If \emph{contour=T}, produces a two-dimensional contour plot of the conditional density estimates.  
Evaluated for an \emph{ngrid} x \emph{ngrid} lattice.  Default is \emph{contour=T}.}
  \item{level}{If \emph{level=T}, produces a two-dimensional level plot of the conditional density estimates.  
Evaluated for an \emph{ngrid} x \emph{ngrid} lattice.  Default is \emph{level=F}.}
  \item{wire}{If \emph{wire=T}, produces a three-dimensional plot of the conditional density estimates.  
Evaluated for an \emph{ngrid} x \emph{ngrid} lattice.  Default is \emph{wire=T}.}
  \item{dens}{If \emph{dens=T}, produces a plot showing how \emph{f(y|x)} varies over \emph{y} for given target values of \emph{x}.  
Target values of \emph{x} are provided using the \emph{targetx.dens} or \emph{quantile.dens} options.  Default is dens=F.}
  \item{targetx.dens}{Target values for \emph{x} in the density plots, e.g, \emph{targetx.dens} = c(200,400,600).  Maximum number of entries is 5.  
If \emph{targetx.dens} has more than 5 entries, only the first 5 will be used.  Default is \emph{targetx.dens = NULL}, 
meaning that the target values for \emph{x} are determined by the \emph{quantile.dens} option.}
  \item{quantile.dens}{Quantiles for the target values for \emph{x} in the density plots, e.g, \emph{quantile.dens} = c(.25,.50,.75).  Maximum number of entries is 5.  
If quantile.dens has  more than 5 entries, only the first 5 will be used.  Default is \emph{quantile.dens} = c(.10,.25,.50,.75,.90).}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}


\value{

\item{fx}{The values of \emph{f(x)}, one for each data point.}
\item{fy}{The values of  \emph{f(y)}, one for each data point.}
\item{fxy}{The values of  \emph{f(x,y)}, one for each data point.  The conditional densities are fxy/fx for \emph{x} and fxy/fy for \emph{y}.}

\item{gridmat}{An (\emph{ngrid*ngrid})x3 matrix used to produce the contour, level, and wire maps.  
The first column contains the lattice values for  \emph{x}, the second column contains the lattice values for  \emph{y}, 
and the third column has the estimated values of  \emph{f(y|x)} at the target values for  \emph{x} and  \emph{y}.}
\item{densmat}{The estimated values of  \emph{f(y|x)} for the two-dimensional density graphs produced when \emph{dens = TRUE}.  
If the number of observations in the call to  \emph{condens} is  \emph{n} and the number of entries in  \emph{quantile.dens} is  \emph{nq}, 
then densmat is an  \emph{n} x  \emph{nq} matrix.}
}

\references{
Li, Oi and Jeffrey Scott Racine.  \emph{Nonparametric Econometrics:  Theory and Practice.}  Princeton, NJ:  Princeton University Press, 2007.

Loader, Clive.  \emph{Local Regression and Likelihood.}  New York:  Springer, 1999.

Pagan, Adrian and Aman Ullah.  \emph{Nonparametric Econometrics.} New York:  Cambridge University Press, 1999.
}

\seealso{
  \link{qregcdf}
}
\examples{
data(dupage99)
dupage99$ratio <- dupage99$av/dupage99$price
dupage99$price <- dupage99$price/1000
par(ask=TRUE)
fit <- condens(ratio~price,contour=TRUE,level=TRUE,wire=TRUE,dens=TRUE, 
  targetx.dens=seq(100,500,100), data=dupage99)

}


\details{
The \emph{locfit} package is used to find the target values of \emph{x} for \emph{f(x)} and \emph{y} for \emph{f(y)}. 
The \emph{expand.grid} command is then used to determine the target values of \emph{x} and \emph{y} for \emph{f(x,y)}.  
The \emph{akima} package is used to interpolate \emph{f(x)}, \emph{f(y)}, and \emph{f(x,y)} to the 
full data set and to the grid of target values for the contour, level, and wire plots.

The density functions \emph{f(x)} and \emph{f(y)} are as follows:

\deqn{f(x) = \frac{1}{sd(x)*b*n} \sum_i K ( \frac{ x_i - x}{sd(x)*b} )}{(sd(x)*b)^{-1} mean((x_i - x)/(sd(x)*b)) }
\deqn{f(y) = \frac{1}{sd(y)*b*n} \sum_i K ( \frac{ y_i - y}{sd(y)*b} )}{(sd(y)*b)^{-1} mean((y_i - y)/(sd(y)*b)) }

A product kernel is used for \emph{f(x,y)}:

\deqn{f(x,y) = \frac{1}{sd(x)*b*sd(y)*b*n}\sum_i K ( \frac{ x_i - x}{sd(x)*b} ) K ( \frac{ y_i - y}{sd(y)*b}  ) }{(sd(x)*b*sd(y)*b)^{-1} mean (K ( (x_i - x)/(sd(x)*b) ) K ( (y_i - y)/(sd(y)*b)  ) }

where \emph{b} is the bandwidth and the target points are \emph{x} and \emph{y}.  
The bandwidth, \emph{b}, can be set using the \emph{bandwidth} option.
If \emph{b} = 0 (the default), \emph{sd(x)*b} and \emph{sd(y)*b} are replaced by window values, \eqn{h = quantile(dist, window)}, 
where \eqn{dist = |x_i - x|} or \eqn{dist = |y_i - y|}.
The window size is set using the \emph{window} option.  By default, \emph{window} = .7 and \emph{bandwidth} = 0.
Available kernel weighting functions include the following:

\tabular{lll}{
Kernel \tab Call abbreviation \tab Kernel function K(z) \cr
Rectangular \tab ``rect'' \tab \eqn{\frac{1}{2} I(|z| <1)}{1/2 * I(|z|<1)} \cr
Triangular \tab ``tria'' \tab \eqn{(1-|z|)I(|z|<1)}{(1-|z|) * I(|z|<1)}\cr
Epanechnikov \tab ``epan'' \tab \eqn{\frac{3}{4} (1-z^2) * I(|z| <1)}{3/4 * (1-z^2)*I(|z| < 1)} \cr
Bi-Square \tab ``bisq'' \tab \eqn{\frac{15}{16} (1-z^2)^2 * I(|z| <1)}{15/16 * (1-z^2)^2 * I(|z| < 1)} \cr
Tri-Cube \tab ``tcub'' \tab \eqn{\frac{70}{81} (1-|z|^3)^3 * I(|z| <1)}{70/81 * (1-|z|^3)^3 * I(|z| < 1)} \cr
Tri-Weight \tab ``trwt'' \tab \eqn{\frac{35}{32} (1-z^2)^3 * I(|z| <1)}{35/32 * (1-z^2)^3 * I(|z| < 1)} \cr
Gaussian \tab ``gauss'' \tab \eqn{(2\pi)^{-.5} e^{-z^2/2}}{2pi^{-.5} exp(-z^2/2)} \cr }

The contour, level, and wire plots are produced from the values in \emph{gridmat} using the \emph{lattice} package.
The two-dimensional density graphs produced when \emph{dens=TRUE} are plots of \emph{f(y,x)/f(x)} at given values of \emph{x}.  
By default, the values for \emph{x} are the quantiles given in \emph{quantile.dens}.  
Alternatively, the values of \emph{x} can be specified directly using the \emph{targetx.dens} option.  
The values used to construct the density graphs are stored in \emph{densmat}.
Both \emph{gridmat} and \emph{densmat} are stored by \emph{condens} even if the printing of the graphs is suppressed.
}

\keyword{Nonparametric}
\keyword{Conditional Density}
\keyword{Density Functions}

