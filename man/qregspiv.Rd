\name{qregspiv }
\alias{qregspiv }
\title{IV Estimator for the Spatial AR Quantile Model}
\description{Uses the Kim and Muller (2004) or Chernozhukov and Hansen (2006) method to estimate a quantile version of the spatial AR Model}
\usage{ 
qregspiv(form,wy=NULL,wmat=NULL,inst=NULL,winst=NULL,shpfile=NULL,
  tau=.5,rhomat=NULL,printsariv=FALSE,silent=FALSE,
  nboot=100,alpha=.05,data=NULL) 
 }

\arguments{
  \item{form}{Model formula}
  \item{wy}{The \emph{WY} variable.  Default:  not specified; program attempts to calculate \emph{WY} using \emph{wmat} or \emph{shpfile}.}
  \item{wmat }{Directly enter \emph{wmat} rather than creating it from a shape file.  Default:   not specified.  \emph{wmat} is needed 
unless \emph{WY} is provided and a full instrument list is specified using \emph{inst}. Default:  \emph{W} = NULL.}
  \item{inst}{List of instruments \emph{not} to be pre-multiplied by \emph{W}.  Entered as \emph{inst=~w1+w2 ...}
   Default:  \emph{inst=NULL}.  See \emph{details} for more information.}
  \item{winst}{List of instruments to be pre-multiplied by \emph{W} before use.  Entered as \emph{winst=~w1+w2 ...}  
Default:  \emph{inst=NULL}.  See \emph{details} for more information.}
  \item{shpfile}{Shape file used to construct \emph{W} based on first order contiguity using a queen criterion.  Needed if \emph{wmat} is not
 provided when the program requires it.  }
  \item{tau}{The quantile.  Default:  \emph{tau} = .5}
  \item{rhomat}{A vector of values for \eqn{\rho}.  If \emph{rhomat}=NULL, uses the Kim and Muller (2004) two-stage approach to estimate the model.
If \emph{rhomat} is a vector with two or more entries, uses the Chernozhukov and Hansen (2006) IV approach to estimate the model.
Default:  \emph{rhomat=NULL}.}
  \item{printsariv}{If \emph{TRUE}, also estimates a standard spatial AR model using an IV approach.  Instruments for WY are based on 
the \emph{inst} and \emph{winst} options.  Default:  \emph{printsariv=FALSE}.}
  \item{silent}{If \emph{TRUE}, no results are printed.  Useful for Monte Carlo.}
  \item{nboot}{The number of simulations for the bootstrap standard errors.  Needed for the Kim and Muller (2004) model.  Default:  \emph{nboot}=100.}
  \item{alpha}{Probability for the confidence intervals, calculated by the percentile method for the Kim and Muller (2004) model.  
Default:  \emph{alpha}=.05, i.e., a 95 percent confidence interval.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}


\value{
  A table showing the coefficient estimates, standard errors, and z-values. 
Also includes a 1-\emph{alpha} confidence interval based on the percentile method when the Kim and Muller is estimated, i.e., when \emph{rhomat = NULL}.
}
\details{
The procedure is intended for quantile estimation of the spatial AR model, \eqn{Y = \rho WY + X \beta + u}.  
It can also be for quantile IV estimation of any model with one endogenous explanatory variable.

\emph{Kim and Muller (2004):}

The Kim and Muller 2004 estimation procedure is the default.  The procedure has two stages.  In the first stage, an instrumental variable is constructed for \emph{WY} using the predicted values from a quantile
regression of \emph{WY} on a set of instruments, \emph{Z}. The second stage is a quantile regression of \emph{Y} on \emph{X} and the predicted values of \emph{WY}.
The same quantile, \emph{tau}, is used for both regressions.

Standard errors are calculated using a simple bootstrap estimator.  New samples are constructed by drawing with replacement from the rows of the
data frame holding \emph{y}, \emph{WY}, \emph{X}, and \emph{Z}.  Both stages are re-estimated \emph{nboot} times using the series of bootstrap samples.
The bootstrap standard errors are the standard deviations of the \emph{nboot} re-calculations of the coefficient estimates.
The confidence intervals are based on the percentile method:  for any coefficient \emph{b}, the 1-\emph{alpha} confidence interval is
  (\emph{quantile(b, alpha/2)}, \emph{quantile(b, 1-alpha/2)}).

\emph{Chernozhukov and Hansen (2006):}

The Chernozhukov and Hansen (2006) procedure is used when a vector of possible values for \eqn{\rho} is specified using the \emph{rhomat} option, 
e.g., \emph{rhomat} = seq(0,.9,.05).
The \emph{qregspiv} command implements a simple version of the Chernzhukov and Hansen estimator in which the explanatory variable \emph{WY} is replaced by the 
predicted values from an OLS regression of \emph{WY} on the instruments, \emph{Z}. This instrumental variable is then used as an explanatory variable for a series
of quantile regressions of \eqn{Y - \rho WY} on \emph{X} and \eqn{\hat{WY}} -- one regression for each value of \eqn{\rho} listed in \emph{rhomat}.
The estimated value of \eqn{\rho} is the value that leads the coefficient on \eqn{\hat{WY}} to be closest to zero. 
After finding \eqn{\hat{\rho}}, the estimated values of \eqn{\beta} are calculated by a quantile regression of \eqn{Y - \hat{\rho} WY} on \emph{X}.

Standard errors are based on equations 3.13 and 3.14 in Chernozhukov and Hansen (2006).  
Let \emph{e} represent the residuals from the quantile regression of \eqn{Y - \hat{\rho} WY} on \emph{X}, and define \eqn{f_i = I(|e_i|<h)/(2h)}, 
where \eqn{h = 1.06*sd(e)*n^{-.2}}.  Also, let \eqn{\Phi} represent the predicted value of \emph{WY} from an OLS regression of \emph{WY} on \emph{Z},
and let \emph{D} represent the actual values of \emph{WY}.  Finally, define \eqn{\Phi^*_i = f_i\Phi_i} and \eqn{X^*_i = f_iX_i}.
Then the covariance matrix for \eqn{\hat{\theta} = (\hat{\rho},  \hat{\beta})} is .

\deqn{V(\hat{\theta}) = J(\tau)^{-1}S(\tau)(J^t)^{-1}}

where \eqn{J(\tau) = \left(\begin{array}{ll}\Phi^{*t} D & \Phi^{*t} X  \\X^{*t} D & X^{*t} X \end{array}\right) }
and \eqn{S(\tau) = \tau(1-\tau)\left(\begin{array}{ll}\Phi^t D & \Phi^t X  \\X^t D & X^t X \end{array}\right) }


\emph{Instruments:}

By default, the instrument list includes \emph{X} and \emph{WX}, where \emph{X} is the original explanatory variable list and \emph{W} is the spatial weight matrix.  
It is also possible to directly specify the full instrument list or to include only a subset of the \emph{X} variables in the 
list that is to be pre-multiplied by \emph{W}. The results of both the quantile IV estimator and the standard IV estimator can be quite sensitive to the 
choice of instruments for the spatial AR model.

Let \emph{list1} and \emph{list2} be user-provided lists of the form \emph{list=~z1+z2}.  
The combinations of defaults (\emph{NULL}) and lists for \emph{inst} alter the final list of instruments as follows: \cr

\emph{inst = NULL},  \emph{winst = NULL}:   \emph{Z = (X, WX)} \cr
\emph{inst = list1}, \emph{winst = NULL}:   \emph{Z = list1} \cr
\emph{inst = NULL},  \emph{winst = list2}:  \emph{Z = (X, W*list2)} \cr
\emph{inst = list1}, \emph{winst = list2}:  \emph{Z = (list1, W*list2)} \cr

Note that when \emph{inst=list1} and \emph{winst=NULL} it is up to the user to specify at least one variable in \emph{list1} that is not also included in \emph{X}.


\emph{Non-Quantile IV Estimates:}

Standard, non-quantile IV estimates are presented if \emph{printsariv = T}.  The first stage is a regression of \emph{WY} on \emph{Z}.  
The second stage is a regression of \emph{Y} on \emph{X} and the predicted values of \emph{WY}.  
Let \eqn{ \hat{G} } be the matrix of explanatory variables in the second stage (i.e., \eqn{ \hat{G} = (Z,  \hat{WY} )} ), 
the covariance matrix is \eqn{ \hat{ \sigma }^2 (\hat{G}'\hat{G})^{-1} }, where \eqn{ \hat{ \sigma }^2 = e'e/n} and
\eqn{e = Y - \hat{\rho} WY - X \hat{\beta}}.  
Note that the variance calculation uses actual values of \emph{WY} while \eqn{(\hat{G}'\hat{G})^{-1}} uses predicted values.
 }

\references{
Chernozhukov, Victor and Christian Hansen, "Instrumental Quantile Regression Inference for Structural and Treatment Effect Models,"
  \emph{Journal of Econometrics} 132 (2006), 491-525.

Kim, Tae-Hwan and Christophe Muller, "Two-Stage Quantile Regression when the First Stage is Based on
Quantile Regression, \emph{Econometrics Journal} 7 (2004), 218-231.

Koenker, Roger.  \emph{Quantile Regression}.  New York:  Cambridge University Press, 2005.

Kostov, Philip, "A Spatial Quantile Regression Hedonic Model of Agricultural Land Prices," \emph{Spatial Economic Analysis} 4 (2009), 53-72.

Zietz, Joachim, Emily Norman Zietz, and G. Stacy Sirmans, "Determinants of House Prices:  A Quantile Regression Approach,"
\emph{Journal of Real Estate Finance and Economics} 37 (2008), 317-333.

 }


\examples{
data(matchdata)
set.seed(4849189)
mdata <- matchdata[matchdata$year==2005,]
obs <- sample(seq(1,nrow(mdata)),400)
mdata <- mdata[obs,]
mdata$age <- 2005 - mdata$yrbuilt
lmat <- cbind(mdata$longitude,mdata$latitude)

fit <- makew(coormat=lmat,method="ring",ringdist=.50)
wmat <- fit$wmat
form <- lnprice~lnland+lnbldg
fit <- qregspiv(form,wmat=wmat,data=mdata,tau=.5)
}


\seealso{
  \code{\link{sarml}}

  \code{\link{qregbmat}}

  \code{\link{qregsim1}}

  \code{\link{qregsim2}}

  \code{\link{qregcpar}}

  \code{\link{qreglwr}}

}


\keyword{Quantile Regression}
\keyword{Spatial AR Model}

