\name{matchmahal}
\alias{matchmahal}
\title{Matched sample data frame based on mahalanobis distances}
\description{Creates a matched sample data frame based on mahalanobis distances}
\usage{
matchmahal(form,data=NULL,discard="none",
  distance="logit",m.order="none",nclose=0,ytreat=1)
}
\arguments{
  \item{form}{Model formula}

  \item{data}{A data frame containing the data. Default:  use data in the current working directory}

  \item{discard}{Observations to be discarded based on the propensity score.
               If \emph{discard} = "control", only control observations are discarded.  
               If \emph{discard} = "treat", only treatment observations are discarded.
               If \emph{discard} = "both", both control and treatment observations are deleted. 
Default:  \emph{discard} = "none"; no options are discarded and propensity scores are not estimated.
  }

  \item{distance}{The link formula to be passed on to the glm command if discard = "control", "treat", or "both"; default = "logit"}

  \item{m.order}{Order by which estimated distances are sorted before starting the matching process.  
  Options:  "decreasing", "increasing", "random", and "none". 
  As the "decreasing" and "increasing" options are based on propensity scores, they are only applicable when discard = "control", "treat", or "both".}

  \item{nclose}{If \emph{nclose}>0, sorts the matched observations by the distance measure and chooses the \emph{nclose} matches with the smallest distances.}

  \item{ytreat}{The value of the dependent variable for the treatment group.  Default:  \emph{ytreat} = 1.  
Constructs matched samples for all other values of the dependent variable.  
If \emph{discard}="treat" or \emph{discard}="both", only treatment observations that were discarded 
for \emph{every} control value of the dependent variable are omitted from the final data set. }
}

\value{

  Returns the matched sample data frame.  Adds the following variables to the data set:

  \emph{origobs}:  The observation number in the original data set

  \emph{matchobs}: The observation number in the matched data set to which the observation is matched.  
\emph{matchobs} refers to the observation's number in the original data set, i.e., to the variable \emph{origobs}.

\emph{Note:}  If the original data set includes variables named \emph{origobs} and \emph{matchobs}, 
they will be overwritten by the variables produced by \emph{matchmahal}.

}

\examples{

set.seed(189)
n = 1000
x <- rnorm(n)
x <- sort(x)
y <- x*1 + rnorm(n, 0, sd(x)/2)
y <- ifelse(y>0,1,0)
table(y)
fit <- matchmahal(y~x,ytreat=1)
table(fit$y)

}


\details{

Creates a matched sample data set by matching each treatment variable to the closest control variable based on mahalanobis distances.
Like \emph{matchprop},  \emph{matchmahal} is particularly useful for creating a series of matched sample data sets over time 
relative to a base time period. 

Let \emph{X1} be the matrix of explanatory variables for the treatment observations and let \emph{X2} be the comparable matrix for the control observations.
The mahalanobis measure of distance between the \emph{i}th row of \emph{X1} and all control observations is \eqn{d_i = mahalanobis(X2, X1[i,], cov(rbind(X2,X1)))}.  
The first observation of \emph{X1} is matched with the closest observation in \emph{X2} based on this distance measure.  
The row is then removed from \emph{X2} and the second observation of \emph{X1} is matched with the closest of the remaining control observations.  
The process is repeated until there are no more observations left in one of the matrices.

By default, \emph{matchprop} matches every treatment observation with a control observation.  If the number of treatment observations (n1) is less than the 
number of control observations (n2), then the first n2 treatment observations will be in the final matched sample data set.  
By default, the observations are matched in the order in which they appear in the original data set. 
Alternatively, the observations can be matched in random order by specifying \emph{m.order} = "random".

The \emph{distance} option allows the user to specify a metric by which observations are determined to be outside the probability support.
The same options are available as in the \emph{matchprop} command.
The natural one is \emph{distance} = "mahal" combined with \emph{discard} = "control", "treat" or "both" and \emph{m.order} = "increasing", "decreasing", or "random".
Other options are listed in the documentation for the \emph{matchprop} command, e.g., \emph{distance} = "logit" or "probit".
Any of the these \emph{distance} options produces a propensity score, \emph{p}.  When \emph{distance} = "mahal", the propensity score is the 
mahalanobis distance of each observation from the vector of means. 
The \emph{discard} option determines how observations are handled that are outside the probability support.  
For example, if the treatment is set to \emph{ytreat} = 1 and the alternative value of the dependent variable is \emph{y} = 2, then:

\emph{discard} = "control":	observations with p[y==2]<min(p[y==1]) are discarded from the y==2 sample

\emph{discard} = "treat":	observations with p[y==1]>max(p[y==2]) are discarded from the y==1 sample

\emph{discard} = "both":	both sets of observations are deleted

If \emph{discard} = "treat" or "both" and the dependent variable has more than two values, 
a different set of treatment observations may be discarded as being outside the support of the two propensity measures.
Only treatment observations that are rejected by \emph{both} models will end up being omitted from the final data set.

}

\references{
Deng, Yongheng, Sing Tien Foo, and Daniel P. McMillen, "Private Residential Price Indices in Singapore," \emph{Regional Science and Urban Economics}, 
42 (2012), 485-494.

Ho, D., Imai, K., King, G, Stuart, E., "Matching as Nonparametric Preprocessing for Reducing Model Dependence in Parametric Causal Inference," \emph{Political Analysis} 15 (2007), 199-236.

Ho, D., Imai, K., King, G, Stuart, E., "MatchIt:  Nonparametric preprocessing for parametric causal inference," \emph{Journal of Statistical Software} 42 (2011), 1-28..

McMillen, Daniel P., "Repeat Sales as a Matching Estimator," \emph{Real Estate Economics} 40 (2012), 743-771.
}

\seealso{
  \code{\link{matchprop}}

  \code{\link{matchqreg}}
}

\keyword{Matching}
