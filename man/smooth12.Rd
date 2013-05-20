\name{smooth12}
\alias{smooth12}
\title{Univariate or Bivariate Interpolation }
\description{Uses the Akima (1970) method for univariate interpolation and the Modified Shephard Algorithm for bivariate interpolation.}
\usage{ smooth12(x,y,xout,knum=16,std=TRUE) }

\arguments{
  \item{x}{The actual values of the x-variable(s).  A simple numeric variable for univariate interpolation and a matrix of locations for bivariate interpolation.}
  \item{y}{The variable to be interpolated.}
  \item{xout}{Points on the x-axis where the function is to be evaluated. A single numeric variable in the case of univariate interpolation 
and a matrix of locations for bivariate interpolation.}
  \item{knum}{The number of target points used for bivariate interpolation.}
  \item{std}{If \emph{TRUE}, re-scales the columns of \emph{x} and \emph{xout} by dividing by the standard deviation of the columns in \emph{x}.  
Not applicable for univariate interpolation.}
}


\value{The values of y interpolated to the \emph{xout} locations.}
\details{
 The univariate version of the function is designed as a partial replacement for the \emph{aspline} function in the \emph{akima} package.  
 It produces a smooth function that closely resembles the interpolation that would be done by hand.
 Values of y are averaged across any ties for x.  The function does not allow for extrapolation beyond the min(x), max(x) range.


 The bivariate version of the function uses the modifed Shepard's method for interpolation.  The function uses the RANN package to find the 
nearest \emph{knum} target points to each location in \emph{xout}.  The following formula is used to interpolate from these target points to 
the locations given in \emph{xout}:

\deqn{\frac{\sum_{i=1}^{knum} w_i y_i}{\sum_{i=1}^{knum} w_i} }

where 

\eqn{w_{i} = ((maxd - d_{i})^2)/(maxd*d_{i})}

and

\eqn{maxd = max(d_{1}, ..., d_{knum})}.



}

\references{
Akima, Hiroshi, "A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures,"
\emph{Journal of the Association for Computing Machinery} 17 (1970), 589-602.

Franke, R. and G. Neilson, "Smooth Interpolation of Large Sets of Scatter data,"
\emph{International Journal of Numerical Methods in Engineering} 15 (1980), 1691 - 1704.
}

\examples{
set.seed(484849)
n = 1000
x <- runif(n,-2*pi,2*pi)
x <- sort(x)
y <- sin(x) + cos(x) - .5*sin(2*x) - .5*cos(2*x) + sin(3*x)/3 + cos(3*x)/3
x1 <- seq(-2*pi,2*pi,length=100)
y1 <- sin(x1) + cos(x1) - .5*sin(2*x1) - .5*cos(2*x1) + sin(3*x1)/3 + cos(3*x1)/3

yout <- smooth12(x1,y1,x)
plot(x,y,type="l")
lines(x,yout,col="red")

x <- seq(0,10)
xmat <- expand.grid(x,x)
y <- sqrt((xmat[,1]-5)^2 + (xmat[,2]-5)^2)
xout <- cbind(runif(n,0,10),runif(n,0,10))
y1 <- sqrt((xout[,1]-5)^2 + (xout[,2]-5)^2)
y2 <- smooth12(xmat,y,xout)
cor(y1,y2)
}

\seealso{
 \link{maketarget}
}


