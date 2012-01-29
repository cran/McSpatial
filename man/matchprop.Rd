\name{matchprop}
\alias{matchprop}
\title{Matched sample data frame based on propensity scores}
\description{Creates a matched sample data frame based on propensity scores}
\usage{
matchprop(form,data=NULL,distance="logit",discard="both",
  reestimate="FALSE",m.order="none",nclose=0,ytreat=1)
}
\arguments{
  \item{form}{Model formula}

  \item{data}{A data frame containing the data. Default:  use data in the current working directory}

  \item{distance}{The link formula to be passed on to the glm command -- usually "probit" or "logit", but other standard options also work.}

  \item{discard}{Observations to be discarded based on the propensity score or the value of the mahalanobis distance measure if \emph{distance}="mahal".
               If \emph{discard} = "control", only control observations are discarded.  
               If \emph{discard} = "treat", only treatment observations are discarded.
               If \emph{discard} = "both", both control and treatment observations are deleted. 
  }
  \item{reestimate}{If \emph{reestimate}=TRUE, the propensity score is reestimated after observations are discarded}

  \item{m.order}{Order by which estimated distances are sorted before starting the matching process.  
  Options:  "decreasing", "increasing", "random", and "none".   }

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
they will be overwritten by the variables produced by \emph{matchprop}.

}

\seealso{
  \code{\link{matchmahal}}

  \code{\link{matchqreg}}
}

\examples{

set.seed(189)
n = 1000
x <- rnorm(n)
x <- sort(x)
y <- x*1 + rnorm(n, 0, sd(x)/2)
y <- ifelse(y>0,1,0)
table(y)
fit <- matchprop(y~x,m.order="largest",ytreat=1)
table(fit$y)

}


\details{
  Creates a matched sample data set using procedures based on the \emph{MatchIt} program.
\emph{MatchIt}'s routines are generally preferable for creating a single matched data set, although by doing less the \emph{matchprop} command is 
somewhat faster than \emph{MatchIt}.  \emph{matchprop} is particularly useful for creating a series of matched sample data sets over time 
relative to a base time period. 

Unless \emph{distance} = "mahal", the \emph{glm} command is used to estimate the propensity scores using
a series of discrete choice models for the probability, \emph{p}, 
that the dependent variable equals \emph{ytreat} rather than each alternative value of the dependent variable.  
The default link function is \emph{distance = "logit"}.  
Alternative link functions are specified using the \emph{distance} option. 
Links include the standard ones for a \emph{glm} model with \emph{family = binomial}, e.g., "probit", "cauchit", "log", and "cloglog".  

If \emph{mahal}= T, \emph{matchprop} implements \emph{MatchIt's} version of mahalanobis matching.  
Letting \emph{X} be the matrix of explanatory variables specified in \emph{form}, the mahalanobis measure of distance from the vector of mean values is
\eqn{p = mahalanobis(X, colMeans(X), cov(X)) }.  Although this version of mahalanobis matching is fast, it may not be the best way to construct matches 
because it treats observations that are above and below the mean symmetrically.  
For example, if \emph{X} is a single variable with \eqn{mean(X)} = .5 and \eqn{var(X)} = 1,
mahalanobis matching treats \emph{X} =.3 and \emph{X} = .7 the same:  mahalanobis(.3,.5,1) = .04 and mahalanobis(.7,.5,.1) = .04.  
The function \emph{matchmahal} is slower but generally preferable for mahalanobis matching because it pairs each treatment observation with 
the closest control observation, i.e., \eqn{min(mahalanobis(X0, X1[i,], cov(X)))}, where \emph{X0} is the matrix of explanatory variables for the control observations,
\emph{X1} is the matrix for the treatment observations, \emph{X} is the pooled explanatory variable matrix, and \emph{i} is the target treatment observation.

To illustrate how \emph{matchprop} constructs matched samples, suppose that the dependent variable takes on three values, 
\emph{y} = 1, 2, 3, and assume that \emph{y} = 1 is the treatment group.  First, the \emph{y} = 1 and \emph{y} = 2 observations are pooled 
and a propensity score \emph{p} is constructed by, e.g., estimating a logit model for the probability that \emph{y} = 1 rather than 2.  
Unless \emph{m.order} = "none", the data frame is then sorted by \emph{p} -- from largest to smallest if \emph{m.order} = "largest", 
from smallest to largest if \emph{m.order} = "smallest", and randomly if \emph{m.order} = "random".  The first treatment observation is then paired with 
the closest control observation, the second treatment observation is paired with the closest of the remaining control observations, and so on 
until the last observation is reached for one of the groups.  No control observation is matched to more than one treatment observation,
and only pairwise matching is supported using \emph{matchprop}. The process is then repeated using the \emph{y} = 1 and \emph{y} = 3 observations.
If the number of treatment observations is n1, the final data set will have roughly 3*n1 observations -- the n1 treatment observations and n1 observations
each from the \emph{y} = 2 and \emph{y} = 3 observations.  The exact number of observations will differ depending on how observations are treated by the
\emph{discard} option, and there will be fewer than n1 observations for, e.g., group 2 if n2<n1.


The \emph{discard} option determines how observations are handled that are outside the probability support.  
In the above example, let p be the propensity score for the logit model for the probability that \emph{y} = 1 rather than 2.  
If \emph{discard} = "control", observations with p[y==2]<min(p[y==1]) are discarded from the y==2 sample.  
If \emph{discard} = "treat", observations with p[y==1]>max(p[y==2]) are discarded from the y==1 sample.
If \emph{discard} = "both", both sets of observations are deleted. 
The process is then repeated for the \emph{y} = 1 and \emph{y} = 3 observations.  
If \emph{discard} = "treat" or "both", a different set of treatment observations may be discarded as being outside the support of the two propensity measures.
Only treatment observations that are rejected by \emph{both} models will end up being omitted from the final data set.

If \emph{reestimate} = T, the propensity scores are reestimated after any observations are discarded.  Otherwise, matches are based on the original propensity scores.

}

\references{
Deng, Yongheng, Sing Tien Foo, and Daniel P. McMillen, "Private Residential Price Indices in Singapore," \emph{Regional Science and Urban Economics}, (forthcoming).

Ho, D., Imai, K., King, G, Stuart, E., "Matching as Nonparametric Preprocessing for Reducing Model Dependence in Parametric Causal Inference," \emph{Political Analysis} 15 (2007), 199-236.

Ho, D., Imai, K., King, G, Stuart, E., "MatchIt:  Nonparametric preprocessing for parametric causal inference," \emph{Journal of Statistical Software} 42 (2011), 1-28..

McMillen, Daniel P., "Repeat Sales as a Matching Estimator," \emph{Real Estate Economics} (forthcoming).
}

\keyword{Matching}
