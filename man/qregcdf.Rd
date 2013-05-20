\name{qregcdf}
\alias{qregcdf}
\title{
Nonparametric quantiles based on conditional CDF functions
}
\description{
Estimates conditional quantile functions based on nonparametric conditional CDF functions
}
\usage{ 
qregcdf(form,taumat=c(.10,.25,.50,.75,.90),hx=0,hy=0,nx=20,ny=100,
  targetx=0,targety=0,graph.target=FALSE,graph.yhat=FALSE, data=NULL)
}

\arguments{
  \item{form }{Model formula}
  \item{taumat }{Vector of quantiles.   Default:  taumat=c(.10, .25, .50, .75, .90) }
  \item{hx }{Bandwidth for \emph{x} in \eqn{K((X_i-x)/hx)}.      Default:
\emph{hx} = 1.06*min(sd(x), (quantile(x,.75)-quantile(x,.25))/1.349)*(n^{-.2})
 }
  \item{hy }{Bandwidth for \emph{y} in \eqn{\Phi((y-Y_i)/hy)}       Default:
\emph{hy} = 1.06*min(sd(y), (quantile(y,.75)-quantile(y,.25))/1.349)*(n^(-.2))
 }
  \item{nx }{Number of target points for \emph{x} if using an evenly spaced grid.  Default is \emph{nx} = 20.  
If \emph{nx}>0, then \emph{targetx} <- seq(min(x),max(x),length=nx) }
  \item{ny }{Number of target points for \emph{y} if using an evenly spaced grid.  Default is \emph{ny} = 200.  
If \emph{ny}>0, then \eqn{targety <- seq(min(y),max(y),length=ny)} }
  \item{targetx }{Vector of user-provided target values for \emph{x}.  An alternative to using \emph{nx} to specify a uniformly spaced grid.  Default:  not specified. }
  \item{targety }{Vector of user-provided target values for \emph{y}.  An alternative to using \emph{nxy} to specify a uniformly spaced grid.  Default:  not specified. }
  \item{graph.target }{If \emph{graph.target=T}, graph of results is produced based on target values.
Default:  \emph{graph.target=F}.
 }
  \item{graph.yhat }{If \emph{graph.yhat=T}, graph of results is produced based on interpolations to actual values of \emph{x} and \emph{y}.  Default:  \emph{graph.yhat=F}. }
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}

\value{
  \item{yhat }{Matrix of quantile predictions.  Dimension is \emph{n x length(taumat)} }
  \item{yhat.target }{Matrix of quantile predictions at target values of \emph{x}.  Dimension is \emph{length(targetx) x length(taumat)}.
 }
  \item{targetx }{Vector of target values for \emph{x}. }
  \item{targety }{Vector of target values for \emph{y}. }
  \item{taumat }{Vector of target quantile values. }
  \item{hx }{Bandwidth for \emph{x}. }
  \item{hy }{Bandwidth for \emph{y}. }
}
\details{ 
Following Li and Racine (2007), equation 6.3, a smoothed version of the conditional CDF of \emph{y} given a value of \emph{x} can be written: \cr

\deqn{F(y|x) = \frac{ \frac{1}{n} \sum_i \Phi (\frac{y-y_i}{h_y})h_x^{-1}K(\frac{X_i-x}{h_x})}{h_x^{-1}\sum_i K(\frac{X_i-x}{h_x})} }{F(y|x) = (1/n) \{ \sum_i \Phi ((y-y_i)/hy)hx^{-1}K((X_i-x)/hx)\}/\{hx^{-1}\sum_i K((X_i-x)/hx)\} }

The estimation procedure begins by evaluating this expression at each point in the grid determined by the values of \emph{target.x} and \emph{target.y}   

The result is an \emph{nx x ny} matrix of values for \emph{F(y|x)}.  
Let \eqn{f(x_j)}{fxj} represent the \emph{ny}-vector of values of \eqn{F(y|x_j)}{F(y|xj)}, and let \eqn{f_k}{fk} indicate the entry of \eqn{F(y|x_j)}{F(y|xj)} 
associated with \eqn{y_k}{yk}, k = 1,...,ny.  Finally, let \eqn{\tau}{tau} represent an entry of \emph{taumat}.  
Then the value of \emph{yhat.target} associated with quantile \eqn{\tau}{tau}
and \eqn{x_j}{xj} is the largest value of \eqn{f_k}{fk} such that \eqn{f_k < \tau < f_k+1}{fk < tau < fk+1}.  
The resulting \emph{nx x length(taumat)} matrix is available after estimation as \emph{yhat.target}.  
The \emph{smooth12} is used to interpolate each column of \emph{yhat.target} to span the full vector of original values of \emph{x}.  
The result is the \emph{n x length(taumat)} matrix \emph{yhat}.

Note:  The default bandwidth may prove too small if there are regions where \emph{x} is sparse.  
It may be necessary to experiment with larger bandwidths for \emph{hx} and \emph{hy}.  
The function \emph{qreglwr} is more flexible, allowing nearest neighbor approaches as well as fixed bandwidths.

}
\references{
Li, Oi and Jeffrey Scott Racine.  \emph{Nonparametric Econometrics:  Theory and Practice}.  Princeton, NJ:  Princeton University Press, 2007.  Chapter 6.
 }

\seealso{
  \link{condens}
}

\examples{
data(dupage99)
dupage99$ratio <- dupage99$av/dupage99$price
o <- order(dupage99$price)
dupage99 <- dupage99[o,]
attach(dupage99)
price <- price/1000

fit <- qregcdf(ratio~price)
ymin = min(fit$yhat)
ymax = max(fit$yhat)
plot(price, fit$yhat[,1],type="l",xlab="Sales Price (1000s)",ylab="Assessment Ratio",
  ylim=c(ymin,ymax),main="Nonparametric Conditional CDF Quantile Regression")
for (j in seq(2,5)) {
  lines(price,fit$yhat[,j])
}
fit$hx
fit$hy

}


\keyword{Quantile Regression}
\keyword{Nonparametric}
\keyword{Conditional Density}

