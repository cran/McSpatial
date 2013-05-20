\name{qregsim2}
\alias{qregsim2}
\title{Machado-Mata Decomposition of Changes in Distributions}
\description{Decomposes quantile regression estimates of changes in the distribution of a dependent variable into the 
components associated with changes in the distribution of the explanatory variables and the coefficient estimates.   }
\usage{ 
qregsim2(formall, formx, dataframe1, dataframe2, bmat1, bmat2,  
  graphx=TRUE, graphb=TRUE, graphy=TRUE, graphdy=TRUE, nbarplot=10, 
  yname=NULL, xnames=NULL, timenames=c("1","2"),
  leglocx="topright",leglocy="topright",leglocdy="topright",
  nsim=20000, bwadjx=1,bwadjy=1,bwadjdy=1)  
 }

\arguments{
  \item{formall}{Model formula.  Must match the model formula used for \emph{qregbmat}.}
  \item{formx}{Model formula for the variables used for the decompositions, e.g., \emph{formx=~x1+x2}.  The coefficients and variables for the other variables
are held at their time 2 values for the simulations.}
  \item{dataframe1}{The data frame for regime 1.  Should include all the variables listed in \emph{formall}.}
  \item{dataframe2}{The data frame for regime 2.  Should include all the variables listed in \emph{formall}.}
  \item{bmat1}{Matrix of values for regime 1 quantile coefficient matrices; the output from running \emph{qregbmat} using \emph{dataframe1}.}
  \item{bmat2}{Matrix of values for regime 2 quantile coefficient matrices; the output from running \emph{qregbmat} using \emph{dataframe2}.}

  \item{graphx}{If \emph{graphx=T}, presents kernel density estimates of each of the explanatory variables in \emph{formx}.}
  \item{graphb}{If \emph{graphb=T}, presents graphs of the quantile coefficient estimates for the variables in \emph{formx}.}
  \item{graphy}{If \emph{graphy=T}, presents of the predicted values of \emph{y} for time1, time2, and the counterfactual.}
  \item{graphdy}{If \emph{graphdy=T}, presents graphs of the changes in densities.}
  \item{nbarplot}{Specifies the maximum number of values taken by an explanatory variable before bar plots are replaced by smooth kernel
density functions.  Only relevant when \emph{graphx = T}. }

  \item{yname}{A label used for the dependent variable in the density graphs, e.g., \emph{yname = "Log of Sale Price"}.}
  \item{xnames}{Labels for graphs involving the explanatory variables, e.g., \emph{xnames = "x1"} for one explanatory variable,
or \emph{xnames = c("x1","x2")} for two variables.}
  \item{timenames}{A vector with labels for the two regimes.  Must be entered as a vector with character values.  Default:  c("1","2").}
  
  \item{leglocx}{Legend location for density plots of the explanatory variables, e.g., \emph{leglocx = "topright"} for one explanatory variable,
or \emph{leglocx = c("topright","topleft")} for two variables.}
  \item{leglocy}{Legend location for density plots of predicted values of the dependent variable.  Default: \emph{leglocy = "topright"}.}
  \item{leglocdy}{Legend location for plot of density changes.  Default: \emph{leglocdy = "topright"}.}

  \item{nsim}{Number of simulations for the decompositions.}
  \item{bwadjx}{Factor used to adjust bandwidths for kernel density plots of the explanatory variables.  Smoother functions are produced when \emph{bwadjust}>1.  
Passed directly to the \emph{density} function's \emph{adjust} option.  Default:  \emph{bwadjx=1}.}
  \item{bwadjy}{Factor used to adjust bandwidths for kernel density plots predicted values of the dependent variable.}
  \item{bwadjdy}{Factor used to adjust bandwidths for plots of the kernel density changes.}
}

\value{
  \item{ytarget}{The values for the x-axis for the density functions.}
  \item{yhat11}{The kernel density function for \eqn{X_1 \beta_1 + Z_1\gamma_1}.}
  \item{yhat22}{The kernel density function for \eqn{X_2 \beta_2 + Z_2\gamma_2}.}
  \item{yhat12}{The kernel density function for \eqn{X_1 \beta_2 + Z_2\gamma_2}.}
  \item{d2211}{The difference between the density functions for \eqn{X_2 \beta_2 + Z_2\gamma_2} and \eqn{X_1 \beta_1 + Z_1\gamma_1}.  
Will differ from yhat22 - yhat11 if \emph{bwadjy} and \emph{bwadjdy} are different.}
  \item{d2212}{The difference between the density functions for \eqn{X_2 \beta_2 + Z_2\gamma_2} and \eqn{X_1 \beta_2 + Z_2\gamma_2}.  
Will differ from yhat22 - yhat12 if \emph{bwadjy} and \emph{bwadjdy} are different.}
  \item{d1211}{The difference between the density functions for \eqn{X_1 \beta_2 + Z_2\gamma_2} and \eqn{X_1 \beta_1 + Z_1\gamma_1}.  
Will differ from yhat12 - yhat11 if \emph{bwadjy} and \emph{bwadjdy} are different.}

}
\details{
The base models are \eqn{y_1 = X_1\beta_1 + Z_1\gamma_1} for regime 1 and \eqn{y_2 = X_2\beta_2 + Z_2\gamma_2} for regime 2.  
The counterfactual model is \eqn{y_{12} = X_1\beta_2 + Z_2\gamma_2}.
The full list of variable (both X and Z) are provided by \emph{form}; this list must correspond exactly with the list provided to \emph{qregbmat}.
The subset of variables that are the subject of the decompositions are listed in \emph{formx}.  


The matrices \emph{bmat1} and \emph{bmat2} are intended to represent the output from \emph{qregbmat}. The models must include the same set of explanatory variables,
and the variables must be in the same order in both  \emph{bmat1} and \emph{bmat2}.
In contrast, the data frames \emph{dataframe1} and \emph{dataframe2} can have different numbers of observations and different sets of explanatory, 
as long as they include the dependent variable and the variables listed in  \emph{bmat1} and \emph{bmat2}.


The output from \emph{qregsim2} is a series of graphs.  If all options are specified, the graphs appear in the following order: 

1.	Kernel density estimates for each variable listed in \emph{formx}.  
Estimated using \emph{density} with default bandwidths and the specified value for \emph{bwadjx}.
Not shown if \emph{graphx=F}.  The \emph{xnames} and \emph{leglocx} options can be used to vary the names used to label the x-axis and the legend location.

2. 	Quantile coefficient estimates for the variables listed in \emph{formx}.  Not listed if \emph{graphb=F}.

3.	Kernel density estimates for the predicted values of \eqn{X_1\beta_1 + Z_1\gamma_1} and \eqn{X_2\beta_2 + Z_2\gamma_2}, 
and the counterfactual, \eqn{X_1\beta_2 + Z_2\gamma_2}.  
Estimated using \emph{density} with default bandwidths and the specified value for \emph{bwadjy}.  
Not shown if \emph{graphy=F}.  The label for the x-axis can be varied with the \emph{yname} option.  
The three estimated density functions are returned after estimation as \emph{yhat11}, \emph{yhat22}, and \emph{yhat12}.

4.	A graph showing the change in densities, d2211 = f22 - f11, along with the Machado-Mata decomposition showing:


(a) the change in densities due to the variables listed in \emph{formx}:  d2212 = f22 - f12.


(b) the change in densities due to the coefficients:  d1211 = f12 - f11.


These estimates are returned after estimation as \emph{d2211}, \emph{dd2212}, and \emph{d1211}.  
The density changes are not shown if \emph{graphdy=F}.  The label for the x-axis can be varied with the \emph{yname} option.
The bandwidth for the original density functions f11, f22, and f12 can be varied with \emph{bwadjdy}.  
It is generally desirable to set \emph{bwadjdy} > \emph{bwadjy} because additional smoothing is needed to make the change in densities appear smooth.



The distributions are simulated by drawing \emph{nsim} samples with replacement from \emph{xobs1 <- seq(1:n1)}, 
\emph{xobs2 <- seq(1:n2)}, and \emph{bobs <- seq(1:length(taumat))}.  The commands for the simulations are:


  xobs1 <- sample(seq(1:n1),nsim,replace=TRUE)

  xobs2 <- sample(seq(1:n2),nsim,replace=TRUE)

  bobs <- sample(seq(1:ntau),nsim,replace=TRUE)

  xhat1 <- allmat1[xobs1,]

  xhat2 <- allmat2[xobs2,]

  znames <- setdiff(colnames(allmat1),colnames(xmat1))

  if (identical(znames,"(Intercept)")) {xhat12 <- xhat1}

  if (!identical(znames,"(Intercept)")) {

    xhat12 <- cbind(xhat2[,znames],xhat1[,colnames(xmat1)])

    xhat12 <- xhat12[,colnames(allmat1)]

  }

  bhat1 <- bmat1[bobs,]

  bhat2 <- bmat2[bobs,]


where \emph{allmat} and \emph{xmat} denote the matrices defined by explanatory variables listed in \emph{formall} (including the intercet) and \emph{formx}.
Since the bandwidths are simply the defaults from the \emph{density} function, 
they are likely to be different across regimes as the number of observations and the standard deviations may vary across times.  
Thus, the densities are re-estimated using the average across regimes of the original bandwidths. 

}

\references{
Koenker, Roger.  \emph{Quantile Regression}.  New York:  Cambridge University Press, 2005.

Machado, J.A.F. and Mata, J., "Counterfactual Decomposition of Changes in Wage Distributions using Quantile Regression," \emph{Journal of Applied Econometrics} 20 (2005), 445-465. 

McMillen, Daniel P., "Changes in the Distribution of House Prices over Time:  Structural Characteristics, Neighborhood or Coefficients?" \emph{Journal of Urban Economics} 64 (2008), 573-589. 

 }

\examples{
par(ask=TRUE)

n = 5000
set.seed(484913)
x1 <- rnorm(n,0,1)
u1 <- rnorm(n,0,.5)
y1 <- x1 + u1

# no change in x.  Coefficients show quantile effects
tau <- runif(n,0,.5)
x2 <- x1
y2 <- (1 + (tau-.5))*x2 + .5*qnorm(tau)

dat <- data.frame(rbind(cbind(y1,x1,1), cbind(y2,x2,2)))
names(dat) <- c("y","x","year")
bmat1 <- qregbmat(y~x,data=dat[dat$year==1,],graphb=FALSE)
bmat2 <- qregbmat(y~x,data=dat[dat$year==2,],graphb=FALSE)
fit1 <- qregsim2(y~x,~x,dat[dat$year==1,],dat[dat$year==2,],
  bmat1,bmat2,bwadjdy=2)

# Distribution of x changes.  Coefficients and u stay the same
x2 <- rnorm(n,0,2)
y2 <- x2 + u1
dat <- data.frame(rbind(cbind(y1,x1,1), cbind(y2,x2,2)))
names(dat) <- c("y","x","year")
bmat1 <- qregbmat(y~x,data=dat[dat$year==1,],graphb=FALSE)
bmat2 <- qregbmat(y~x,data=dat[dat$year==2,],graphb=FALSE)
fit1 <- qregsim2(y~x,~x,dat[dat$year==1,],dat[dat$year==2,],
  bmat1,bmat2,bwadjdy=2)

}


\seealso{
  \code{\link{dfldens}}

  \code{\link{qregbmat}}

  \code{\link{qregsim1}}

  \code{\link{qregcpar}}

  \code{\link{qreglwr}}

}


\keyword{Quantile Regression}
\keyword{Machado-Mata}

