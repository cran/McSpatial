\name{semip}
\alias{semip}
\title{Semi-Parametric Regression}
\description{
Estimates a semi-parametric model with the form \eqn{y = X \beta + f(z) + u}, 
where f(z) is either fully nonparametric with \eqn{f(z) = f(z_1)} or conditionally parametric with \eqn{f(z) = z_2 \lambda (z_1)}. }

\usage{ 
semip(form,nonpar,conpar,window1=.25,window2=.25,bandwidth1=0,bandwidth2=0, 
 kern="tcub",distance="Mahal",alldata=FALSE, print.summary=TRUE, data=NULL)
}

\arguments{
  \item{form }{Model formula.  Specifies the base parametric form of the model, \eqn{y = X \beta}.  
  Any number of variables can be included in X.  Format:  semip(y~x1+x2..., ...). }
  \item{nonpar }{List of variables in \eqn{z_1}.  Formats: semip(..., nonpar=~z1a, ...) or semip(..., nonpar=~z1a+zb, ...).  
  Important:  note the "~" before the first z1 variable.  At most two variables can be included in \eqn{z_1}.
  }
  \item{conpar}{List of variables in \eqn{z_2}.  By default, \emph{conpar} = NULL and f(z) has the fully nonparametric form \eqn{f(z) = f(z_1)};
  in this case the variables in \eqn{z_1} are taken from the list provided by \emph{nonpar}.  If a list of variables is provided for \emph{nonpar},
  the conditionally parametric form  \eqn{f(z) = z_2 \lambda (z_1)} is assumed for f(z), and the variables for \eqn{z_2} are provided by \emph{conpar}.
  Any number of variables can be included in \emph{conpar}.  Format:  semip(..., conpar=~z2a+z2b+z2c+..., ...).  
   Important:  note the "~" before the first \emph{z2} variable. 
  }
  \item{window1 }{Window size for the LWR or CPAR regressions of \emph{y} and \emph{x} on \emph{z}.  Default = .25. }
  \item{window2 }{Window size for the LWR or CPAR regression of \eqn{y-X \hat{\beta}}{y-X \beta}  on \emph{z}.  Default = .25. }
  \item{bandwidth1 }{Bandwidth for the LWR or CPAR regressions of \emph{y} and \emph{x} on \emph{z}.  Default:  not specified.}
  \item{bandwidth2 }{Bandwidth for the LWR or CPAR regression of \eqn{y-X \hat{\beta}}{y-X \beta}  on \emph{z}.  Default:  not specified.}
  \item{kern }{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{distance }{Options:  "Euclid", "Mahal", or "Latlong" for Euclidean, Mahalanobis, or "great-circle" geographic distance.   
May be abbreviated to the first letter but must be capitalized.  
Note:  \emph{semip} looks for the first two letters to determine which variable is latitude and which is longitude, 
so data set must be attached first or specified using the data option; options like data$latitude will not work.  Default:  Mahal. }
  \item{alldata }{If \emph{alldata=T}, each observation is used as a target value for \emph{x}.  
When \emph{alldata=F}, the function is estimated at a set of points chosen by the \emph{locfit} program using an adaptive decision tree approach, 
and the \emph{akima} program is used to interpolate to the full set of observations.  Specifying \emph{alldata=T} can lead to long estimation times.}
  \item{print.summary }{If \emph{print.summary=T}, prints a summary of the regression results for \emph{ey} on \emph{ex}, i.e., 
 the parametric portion of the model.  Default:  \emph{print.summary=T}.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}
}


\value{
  \item{xcoef}{The estimated coefficients for the parametric part of the model, \eqn{\beta}.}
  \item{vmat}{The covariance matrix for the estimates of \eqn{\beta}.}
  \item{xbhat}{The predicted values of \emph{y} for the full data set.}
  \item{nphat}{The predicted values of \emph{f(z)} for the full data set.  
\emph{mean(xbhat)+mean(nphat)} will be close but not necessarily identical to \emph{mean(y)}.}
  \item{nphat.se}{Standard errors for the predicted values of \emph{y} for the full data set.}
  \item{npfit}{The complete set of \emph{lwr} or \emph{cparlwr} results from the nonparametric regression of \eqn{y - X \beta} on Z.}
  \item{df1}{\emph{k + tr(L)}, where \emph{k} is the number of explanatory variables in \eqn{X \beta} (including the constant) and \emph{L} 
  is the \emph{nxn} matrix used to calculate the final-stage nonparametric or conditionally parametric regression of \eqn{Y - X \hat{\beta}}{Y - X \beta}  on \emph{Z}.   
  \emph{df1} is one measure of the degrees of freedom used in estimation.  }
  \item{df2}{An alternative measure of the degrees of freedom used in estimation, \eqn{df2 = k + tr(L'L)}. }
  \item{sig2}{Estimated residual variance, \eqn{sig2 = rss/(n-2*df1+df2)}.}
}
\details{
If \emph{conpar = NULL}, the function implements Robinson's (1988) semi-parametric estimator for the model \eqn{y = X \beta + f(z) + u}.  
In this case, the list of variables in \emph{z} is taken from \emph{nonpar} and \emph{z} can have at most two variables. 
If a list of variables is provided for \emph{conpar}, the function implements the semi-parametric estimator for the model 
\eqn{f(z) = z_2 \lambda (z_1)}.  In this case, the list of variables in \eqn{z1} is taken from \emph{nonpar} and 
the list of variables in \eqn{z_2} is taken from \emph{conpar}.  \eqn{z_1} can have at most two variables.  
There is no limit on the number of variables in \eqn{z_2}.

The estimation procedure has the following three steps under either specification:

1.  Nonparametric regressions of \emph{y} on \emph{z} and each \emph{X} on \emph{z} using the \emph{lwr} function when \emph{conpar=NULL} 
and the \emph{cparlwr} function when a list of variables is provided for \emph{cparlwr}.
The window or bandwidth for these regressions is set by \emph{window1} or \emph{bandwidth1}. 

2.  OLS regression of \eqn{y-\hat{y}}{y-fitted(y)} on the \emph{k-1} variables in \eqn{X - \hat{X}}{X - fitted(X)}, omitting the intercept.  
The coefficients from this regression are the estimated values of \eqn{\beta}.

3. Nonparametric regression of \eqn{y- X \hat{\beta}}{y-X \beta}  on \emph{z} using the \emph{lwr} function when \emph{conpar=NULL} 
and the \emph{cparlwr} function when a list of variables is provided for \emph{cparlwr}. 
The window or bandwidth for these regressions is set by \emph{window2} or \emph{bandwidth2}.

The stage-two OLS regressions use \emph{k} degrees of freedom.  The stage-three nonparametric regression uses \emph{2*df1-df2} degrees of freedom, 
where \eqn{df1 = tr(L)} and \eqn{df2 = tr(L'L)} and \emph{L} is the \emph{nxn} matrix for the \emph{lwr} or \emph{cparlwr} regression 
\eqn{L(Y - X \hat{\beta})}{L(Y - X \beta)}.  
The estimated variance is \eqn{\hat{\sigma^2} = rss/(n-2*df1+df2)}{sig2 = rss/(n-2*df1+df2)}, 
where \eqn{rss = \sum_i(y_i-X_i\beta-f(z_i))^2}{rss = sum(y-XB-f(z))^2} .  The covariance matrix estimate for \eqn{\beta} is 
\eqn{\hat{\sigma^2}((X-\hat{X})'(X-\hat{X}))^{-1}.}{sig2*((X-fitted(X))'(X-fitted(X)))^(-1).}  
The covariance matrix is stored as \emph{vmat}. 

Estimation can be very slow when \emph{alldata=T}.  When \emph{alldata=F}, the package \emph{locfit} is used to find a good set of target points at which to 
evaluate the function.  See Loader (1999, section 12.2) for a description of the algorithm used to determine the target points.  
The \emph{akima} package is then used to interpolate the coefficient estimates, the standard errors, and the values used to form \emph{df1} and \emph{df2}.  

The nonparametric regressions are estimated using either the \emph{lwr} or \emph{cparlwr} function.  See their descriptions for more information.
}

\references{
Cleveland, William S. and Susan J. Devlin, "Locally Weighted Regression:  An Approach to Regression Analysis by Local Fitting," \emph{Journal of the American Statistical Association} 83 (1988), 596-610. 

Loader, Clive.  \emph{Local Regression and Likelihood}.  New York:  Springer, 1999. 

McMillen, Daniel P., "Issues in Spatial Data Analysis," \emph{Journal of Regional Science} 50 (2010), 119-141. 

McMillen, Daniel P., "Employment Densities, Spatial Autocorrelation, and Subcenters in Large Metropolitan Areas," \emph{Journal of Regional Science} 44 (2004), 225-243. 

McMillen, Daniel P. and Christian Redfearn, "Estimation and Hypothesis Testing for Nonparametric Hedonic House Price Functions," \emph{Journal of Regional Science} 50 (2010), 712-733.

Pagan, Adrian and Aman Ullah.  \emph{Nonparametric Econometrics}.  New York:  Cambridge University Press, 1999. 

Robinson, Paul M. 1988.   "Root-N-Consistent Semiparametric Regression," \emph{Econometrica}, 56, 931-954. \cr

}

\examples{

# Single variable in f(z)
par(ask=TRUE)
n = 1000
x <- runif(n,0,2*pi)
x <- sort(x)
z <- runif(n,0,2*pi)
xsq <- x^2
sinx <- sin(x)
cosx <- cos(x)
sin2x <- sin(2*x)
cos2x <- cos(2*x)
ybase1 <-  x - .1*xsq + sinx - cosx - .5*sin2x + .5*cos2x
ybase2 <- -z + .1*(z^2) - sin(z) + cos(z) + .5*sin(2*z) - .5*cos(2*z)
ybase <- ybase1+ybase2
sig = sd(ybase)/2
y <- ybase + rnorm(n,0,sig)

# Correct specification for x; z in f(z)
fit <- semip(y~x+xsq+sinx+cosx+sin2x+cos2x,nonpar=~z,window1=.20,window2=.20)
2*fit$df1 - fit$df2
yvect <- c(min(ybase1,fit$xbhat), max(ybase1, fit$xbhat))
xbhat  <- fit$xbhat - mean(fit$xbhat) + mean(ybase1)
plot(x,ybase1,type="l",xlab="x",ylab="ybase1",ylim=yvect, main="Predictions for XB")
lines(x, xbhat, col="red")

predse <- sqrt(fit$sig2 + fit$nphat.se^2)
nphat <- fit$nphat - mean(fit$nphat) + mean(ybase2)
lower <- nphat + qnorm(.025)*fit$nphat.se
upper <- nphat + qnorm(.975)*fit$nphat.se
o <- order(z)
yvect <- c(min(lower), max(upper))
plot(z[o], ybase2[o], type="l", xlab="z", ylab="f(z) ",
   main="Predictions for f(z) ", ylim=yvect)
lines(z[o], nphat[o], col="red")
lines(z[o], lower[o], col="red", lty="dashed")
lines(z[o], upper[o], col="red", lty="dashed")

rm(list=ls())

# Chicago Housing Sales
data(matchdata)
match05 <- data.frame(matchdata[matchdata$year==2005,])
match05$age <- 2005-match05$yrbuilt

library(maptools)
library(RColorBrewer)
library(akima)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1,]
cmap <- cmap[
cmap$CAREA=="Albany Park"|
cmap$CAREA=="Edgewater"|
cmap$CAREA=="Edison Park"|
cmap$CAREA=="Forest Glen"|
cmap$CAREA=="Jefferson Park"|
cmap$CAREA=="Lincoln Square"|
cmap$CAREA=="North Park"|
cmap$CAREA=="Norwood Park"|
cmap$CAREA=="Rogers Park"|
cmap$CAREA=="Uptown"|
cmap$CAREA=="West Ridge", ]

lmat <- coordinates(cmap)

#fixed effects for community areas
fit1 <- lm(lnprice~lnland+lnbldg+rooms+bedrooms+bathrooms+centair+fireplace+brick+
garage1+garage2+age+rr+dcbd +factor(carea), data=match05)
length(fit1$coef)
attach(match05)
xbhat <- fit1$coef[1] + fit1$coef[2]*lnland + fit1$coef[3]*lnbldg + 
 fit1$coef[4]*rooms +  fit1$coef[5]*bedrooms + fit1$coef[6]*bathrooms + 
 fit1$coef[7]*centair + fit1$coef[8]*fireplace +  fit1$coef[9]*brick + 
 fit1$coef[10]*garage1 + fit1$coef[11]*garage2 + fit1$coef[12]*age +
 fit1$coef[13]*rr + fit1$coef[14]*dcbd
lochat <- fitted(fit1) - xbhat
lmean = mean(lochat)
lochat <- interpp(longitude, latitude, lochat, lmat[,1], lmat[,2], extrap=FALSE, 
	duplicate="mean")$z
cmap$lochat1 <- lochat - lmean + mean(lnprice)
mapplot(cmap,"lochat1",title="Fixed  Effects & dcbd")

# nonparametric control for dcbd
fit2 <- semip(lnprice~lnland+lnbldg+rooms+bedrooms+bathrooms+centair+fireplace+brick+
garage1+garage2+ age+rr, nonpar=~dcbd, data=match05)
2*fit2$df1 - fit2$df2
lochat <- fit2$nphat
lmean = mean(lochat)
lochat <- interpp(longitude, latitude, lochat, lmat[,1], lmat[,2], extrap=FALSE, 
	duplicate="mean")$z
cmap$lochat2 <- lochat - lmean + mean(lnprice)
mapplot(cmap,"lochat2",title="f(dcbd)")

# nonparametric controls for longitude and latitude
fit3 <- semip(lnprice~lnland+lnbldg+rooms+bedrooms+bathrooms+centair+fireplace+brick+
garage1+garage2+ age+rr+dcbd, nonpar=~longitude+latitude, data=match05,
distance="Latlong")
2*fit3$df1 - fit3$df2
lochat <- dcbd*fit3$xcoef["dcbd"] + fit3$nphat
lmean = mean(lochat)
lochat <- interpp(longitude, latitude, lochat, lmat[,1], lmat[,2], extrap=FALSE, 
	duplicate="mean")$z
cmap$lochat3 <- lochat - lmean + mean(lnprice)
mapplot(cmap,"lochat3",title="f(longitude, latitude) & dcbd")

# Conditional parametric model:  y = XB + dcbd*lambda(longitude,latitude) + u
fit4 <- semip(lnprice~lnland+lnbldg+rooms+bedrooms+bathrooms+centair+fireplace+
 brick+garage1+garage2+age+rr, nonpar=~longitude+latitude, conpar=~dcbd, 
 data=match05, distance="Latlong")
2*fit4$df1 - fit4$df2
lochat <- fit4$npfit$yhat
lmean = mean(lochat)
lochat <- interpp(longitude, latitude, lochat, lmat[,1], lmat[,2], extrap=FALSE, 
	duplicate="mean")$z
cmap$lochat4 <- lochat - lmean + mean(lnprice)
mapplot(cmap,"lochat4",title="dcbd*lambda(longitude, latitude)")
cmap$dcbdhat <- interpp(fit4$npfit$target[,1],fit4$npfit$target[,2],
  fit4$npfit$xcoef.target[,"dcbd"],lmat[,1],lmat[,2],duplicate="mean")$z
mapplot(cmap,"dcbdhat",title="dcbd coefficients")


# Conditional parametric model:  y = XB + Z*lambda(longitude,latitude) + u
# Z = (dcbd,lnland,lnbldg,age)
fit4 <- semip(lnprice~rooms+bedrooms+bathrooms+centair+fireplace+brick+
garage1+garage2+rr, nonpar=~longitude+latitude, conpar=~dcbd+lnland+lnbldg+age, 
 data=match05, distance="Latlong")
2*fit4$df1 - fit4$df2
lochat <- fit4$npfit$yhat
lmean = mean(lochat)
lochat <- interpp(longitude, latitude, lochat, lmat[,1], lmat[,2], extrap=FALSE, 
	duplicate="mean")$z
cmap$lochat4 <- lochat - lmean + mean(lnprice)
mapplot(cmap,"lochat4",title="Nonparametric Part of Model")
cmap$dcbdhat <- interpp(fit4$npfit$target[,1],fit4$npfit$target[,2],
 fit4$npfit$xcoef.target[,"dcbd"],lmat[,1],lmat[,2],duplicate="mean")$z
mapplot(cmap,"dcbdhat",title="dcbd coefficients")
cmap$lnlandhat <- interpp(fit4$npfit$target[,1],fit4$npfit$target[,2],
 fit4$npfit$xcoef.target[,"lnland"], lmat[,1],lmat[,2],duplicate="mean")$z
mapplot(cmap,"lnlandhat",title="lnland coefficients")
cmap$lnbldghat <- interpp(fit4$npfit$target[,1],fit4$npfit$target[,2],
 fit4$npfit$xcoef.target[,"lnbldg"], lmat[,1],lmat[,2],duplicate="mean")$z
mapplot(cmap,"lnbldghat",title="lnbldg coefficients")
cmap$agehat <- interpp(fit4$npfit$target[,1],fit4$npfit$target[,2],
 fit4$npfit$xcoef.target[,"age"], lmat[,1],lmat[,2],duplicate="mean")$z
mapplot(cmap,"agehat",title="age coefficients")



}



\seealso{
  \link{cparlwr}

  \link{lwr}
}

\keyword{Semi-Parametric Models}
