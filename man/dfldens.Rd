\name{dfldens }
\alias{dfldens }
\title{Counterfactual Kernel Density Functions}
\description{Uses the DiNardo, Fortin, and Lemieux approach to re-weight kernel density functions based on values of an explanatory variable from an earlier period.}
\usage{ 
dfldens(y,lgtform,window=0,bandwidth=0,kern="tcub",probit=FALSE,
  graph=TRUE,yname="y",alldata=FALSE,data=NULL)
}

\arguments{
  \item{y}{The dependent variable for which the counterfactual density is estimated.  The data frame must be specified if it has not been attached, e.g., 
\emph{y=mydata$depvar.}}
  \item{lgtform}{The formula for the logit or probit model for the time variable.  The dependent variable should be a 0-1 variable with 1's representing
the later time period.  Example:  \emph{lgtform=timevar~x1+x2.}}
  \item{window}{The window size for the kernel density function.  Default:  not used.}
  \item{bandwidth}{The bandwidth.  Default:  bandwidth = (.9*(quantile(y1,.75)-quantile(y1,.25))/1.34)*(n1^(-.20)), 
specified by setting \emph{bandwidth} = 0 and \emph{window} = 0.}
  \item{kern}{Kernel weighting function.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{probit}{If \emph{TRUE}, a probit model is used for the time variable rather than logit.  Default:  \emph{probit = FALSE}.}
  \item{graph}{If \emph{TRUE}, produces a graph showing the density function for time 1 and the counterfactual density.  Default:  \emph{graph=TRUE}.}
  \item{yname}{The name to be used for the variable whose density functions are drawn when \emph{graph=T}.  Default:  \emph{yname = "y"}.}
  \item{alldata}{If \emph{TRUE}, the density functions are calculated using each observation in turn as a target value.  
When \emph{alldata=F}, densities are calculated at a set of points chosen by the \emph{locfit} program using an adaptive decision tree approach, 
and the \emph{akima} program is used to interpolate to the full set of observations.}
  \item{data}{A data frame with the variables for the logit or probit model specified by \emph{lgtform}.  Note: the data frame for \emph{y} must be 
specified even if it is part of \emph{data}.}
}


\value{
  \item{target}{The vector of target values for \emph{y} for the density functions.}
  \item{dtarget1}{The vector of densities in period 1 at the target values of \emph{y}.}
  \item{dtarget10}{The counterfactual densities in period 1 at the target values of \emph{y}.}
  \item{dhat1}{The vector of densities in period 1 at the actual values of \emph{y}.}
  \item{dhat10}{The counterfactual densities in period 1 at the actual values of \emph{y}.}
}
\details{
The \emph{dfldens} command first calculates kernel density estimates for \emph{y} in time period \emph{timevar} = 1.
The density estimate at target point \emph{y} is \eqn{f(y_1) = (1/(hn_1)) \sum_i K((y_{1i} - y_1)/h)}.  
The following kernel weighting functions are available:

\tabular{lll}{

Kernel \tab Call abbreviation \tab Kernel function K(z) \cr
Rectangular \tab ``rect'' \tab \eqn{\frac{1}{2} I(|z| <1)}{1/2 * I(|z|<1)} \cr
Triangular \tab ``tria'' \tab \eqn{(1-|z|)I(|z|<1)}{(1-|z|) * I(|z|<1)}\cr
Epanechnikov \tab ``epan'' \tab \eqn{\frac{3}{4} (1-z^2) * I(|z| <1)}{3/4 * (1-z^2)*I(|z| < 1)} \cr
Bi-Square \tab ``bisq'' \tab \eqn{\frac{15}{16} (1-z^2)^2 * I(|z| <1)}{15/16 * (1-z^2)^2 * I(|z| < 1)} \cr
Tri-Cube \tab ``tcub'' \tab \eqn{\frac{70}{81} (1-|z|^3)^3 * I(|z| <1)}{70/81 * (1-|z|^3)^3 * I(|z| < 1)} \cr
Tri-Weight \tab ``trwt'' \tab \eqn{\frac{35}{32} (1-z^2)^3 * I(|z| <1)}{35/32 * (1-z^2)^3 * I(|z| < 1)} \cr
Gaussian \tab ``gauss'' \tab \eqn{(2\pi)^{-.5} e^{-z^2/2}}{2pi^{-.5} exp(-z^2/2)} \cr }

By default, \emph{dfldens} uses a tri-cube kernel with a fixed bandwidth of h = (.9*(quantile(y1,.75)-quantile(y1,.25))/1.34)*(n1^(-.20)).
The results are stored in \emph{dtarget1} and \emph{dhat1}.

The counterfactual density is an estimate of the density function for \emph{y} in time 1 if the explanatory variables
listed in \emph{lgtform} were equal to their time 0 values. 
DiNardo, Fortin, and Lemieux (1996) show that the the following re-weighting of  \eqn{f(y_1)} is an estimate of the counterfactual density:
\eqn{(1/(hn_1)) \sum_i \tau_i K((y_{1i} - y_1)/h)}.
The weights are given by \eqn{tau_i = (P(x_i)/(1-P(x_i)))/(p/(1-p)) }, where \eqn{p = n_0/(n_0 + n_1))} and 
\eqn{P(x_i))} is the estimated probability that \emph{timevar = 0} from the estimated logit or probit regression of \emph{timevar} on \emph{X}.

If \emph{X} includes a single variable \emph{x}, the counterfactual density shows how the \eqn{f(y_1)} would change if \eqn{x = x_0} rather than \eqn{x_1}.
Alternatively, \emph{X} can include multiple variables, in which case 
the counterfactual density shows how the \eqn{f(y_1)} would change if all of the variables in \emph{X} were equal to their \emph{timevar = 0} values.
 }

\references{
DiNardo, J., N. Fortin, and T. Lemieux, "Labor Market Institutions and the Distribution of Wages, 1973-1992:  A Semi-Parametric Approach,"
\emph{Econometrica} 64 (1996), 1001-1044.

Leibbrandt, Murray, James A. Levinsohn, and Justin McCrary, "Incomes in South Africa after the Fall of Apartheid,"
\emph{Journal of Globalization and Development} 1 (2010).
 }

\examples{
data(matchdata)
matchdata$year05 <- matchdata$year==2005
fit <- dfldens(matchdata$lnprice, year05~lnland+lnbldg, window=.2, 
  yname = "Log of Sale Price", data=matchdata)
matchdata$age <- matchdata$year - matchdata$yrbuilt
fit <- dfldens(matchdata$lnprice, year05~age, window=.2, 
  yname="Log of Sale Price", data=matchdata)

}


\seealso{

  \code{\link{qregsim2}}

}


\keyword{Nonparametric}
\keyword{DiNardo-Levinsohn-McCrary}

