\name{matchdata}
\alias{matchdata}
\docType{data}
\title{
Matched samples of house sales in Chicago for 1995 and 2005
}
\description{
Sales prices, structural characteristics, and location variables for 1602 single-family homes in the City of Chicago in 1995 and a matched sample of 1602 homes in 2005.  
}
\usage{data(matchdata)}
\format{
  A data frame with 3204 observations on the following 18 variables.
  \describe{
    \item{\code{year}}{Year of sale, 1995 or 2005}
    \item{\code{lnland}}{Log of land area in square feet}
    \item{\code{lnbldg}}{Log of building area in square feet}
    \item{\code{rooms}}{Number of rooms}
    \item{\code{bedrooms}}{Number of bedrooms}
    \item{\code{bathrooms}}{Number of bathrooms}
    \item{\code{centair}}{Home has central air conditioning}
    \item{\code{fireplace}}{Home has one or more fireplaces}
    \item{\code{brick}}{Brick or brick/frame construction}
    \item{\code{garage1}}{Garage, 1 car}
    \item{\code{garage2}}{Garage, 2+ cars}
    \item{\code{dcbd}}{Distance from the central business district or "CBD" - the traditional center of Chicago at the intersection of State and Madison Streets, at approximately -87.627800 longitude and 41.881998 latitude}
    \item{\code{rr}}{Within .25 miles of a rail line}
    \item{\code{yrbuilt}}{Year the home was built}
    \item{\code{carea}}{a factor with levels. Community area, a traditional definition of neighborhood in Chicago.}
    \item{\code{latitude}}{Latitude in degrees}
    \item{\code{longitude}}{Longitude in degrees}
    \item{\code{lnprice}}{Log of sales price}
  }
}
\details{
Includes all sales of single-family homes on the Far North Side of Chicago listed in the cleaned Illinois Department of Revenue file for 1995.  
A matched sample is created from comparable 2005 sales using the MatchIt package.  
Matches are created based on propensity scores estimated using a logit model for the probability that a home sold in 2005 rather than 1995.  
The commands used to create the matched sample are the following:

hedonic$carea <- as.factor(hedonic$cname) \cr
m.out <- matchit(y~lnland + lnbldg + rooms + bedrooms + bathrooms + centair + fireplace + brick + garage1 + garage2 + dcbd + elstop + lake + rr + yrbuilt + carea + latitude + longitude,
  data=hedonic,method="nearest",discard="both") \cr
mdata <- match.data(m.out) \cr
attach(mdata) \cr
matchdata <- data.frame(year, lnland, lnbldg, rooms, bedrooms, bathrooms, centair, fireplace, brick, garage1, garage2, 
  dcbd, rr, yrbuilt, carea, latitude, longitude, lnprice) 
\cr



The elstop and lake variables, which are not included here, indicate whether a home is within .25 miles of and EL stop and within .5 miles of Lake Michigan.

}
\source{
Daniel McMillen.  
Sales data were provided originally by the Illinois Department of Revenue.  
Structural characteristics are drawn from the 1997 assessment file from the Cook County Assessor's Office.
}
\references{
Deng, Yongheng, Sing Tien Foo, and Daniel P. McMillen, "Private Residential Price Indices in Singapore," \emph{Regional Science and Urban Economics}, 
42 (2012), 485-494.

Ho, D., Imai, K., King, G, Stuart, E., "Matching as Nonparametric Preprocessing for Reducing Model Dependence in Parametric Causal Inference," \emph{Political Analysis} 15 (2007), 199-236.

Ho, D., Imai, K., King, G, Stuart, E., "MatchIt:  Nonparametric preprocessing for parametric causal inference," \emph{Journal of Statistical Software} 42 (2011), 1-28..

McMillen, Daniel P., "Repeat Sales as a Matching Estimator," \emph{Real Estate Economics} 40 (2012), 743-771.
}

\examples{
data(matchdata)
matchdata$year05 <- matchdata$year==2005
matchdata$age <- matchdata$year - matchdata$yrbuilt 
fit <- lm(lnprice~lnland+lnbldg+rooms+bedrooms+bathrooms+centair+fireplace+brick+
  garage1+garage2+dcbd+rr+age+year05+factor(carea), data=matchdata)
summary(fit)
}

\keyword{datasets}
