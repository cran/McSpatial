\name{misclass}
\alias{misclass}
\title{
Probit Model with Misclassification of the Dependent Variable
}
\description{
Implements the Hausman, Arbrevaya, and Scott-Morton (1998) maximum likelihood estimator for probit models with potential misclassification of the dependent variable.
}
\usage{ 
misclass(form,a0=0,a1=0,bmat=0,print.summary=TRUE,data=NULL) 

}

\arguments{
  \item{form }{Model formula}
  \item{a0 }{Starting value for \eqn{\alpha_0}{a0}.  Default:  \emph{a0} = 0.}
  \item{a1 }{Starting value for \eqn{\alpha_1}{a1}.  Default:  \emph{a1} = 0.}
  \item{bmat }{Starting values for \eqn{\beta}.  Default:  \emph{bmat} = 0, uses standard probit values.}
  \item{print.summary }{If \emph{print.summary=T}, prints a summary of the final \emph{nlm} estimates.  Default:  \emph{print.summary=T}.}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory}

}

\value{
  \item{a0 }{Estimate of \eqn{\alpha_0}{a0}, the probability that a true value of 0 is misclassified as a 1.}
  \item{a1 }{Estimate of \eqn{\alpha_1}{a1}, the probability that a true value of 1 is misclassified as a 0.}
  \item{estimate }{Coefficient estimates.}
  \item{stderr }{Standard errors for \emph{estimate} }
  \item{vmat }{Full covariance matrix.}
  \item{iterations }{The number of iterations taken to convergence.}
  \item{minimum }{The value of the log-likelihood function.}
  \item{gradient }{The gradient vector.}
}
\details{ 
Let y be the observed value of the 0-1 dependent variable and let \eqn{\tilde{y}}{yhat}  be the true value.  
The probability that a 0 is incorrectly classified as a 1 is \eqn{\alpha_0 = Pr(y=1|\tilde{y} = 0)}{a0 = Pr(y=1|yhat=0)}  
and the probability that a 1 is incorrectly classified as a 0 is \eqn{\alpha_1 = Pr(y=0|\tilde{y} = 1)}{a1 = Pr(y=0|yhat=1)} .  
Under the assumption that the errors in the underlying latent variable \eqn{X \beta + u}  are normally distributed, 
the probabilities of observing the correctly classified values of the dependent variable are 
\eqn{Pr(\tilde{y} = 1|X) = \Phi(X \beta)}{Pr(yhat=1|X) = \Phi(X \beta)}  
and \eqn{Pr(\tilde{y} = 0|X) = 1 - \Phi(X \beta)}{Pr(yhat=0|X) = 1 - \Phi(X \beta)}.   
The probability that an observation is classified as a 1 is 
\eqn{Pr(y = 1|X) = (1-\alpha_1)\Phi(X \beta) + \alpha_0 (1-\Phi(X \beta)) =\alpha_0 + (1-\alpha_0-\alpha_1)\Phi(X \beta)}{Pr(y = 1|X) = (1-a1)\Phi(X \beta) + a0 (1-\Phi(X \beta)) =a0 + (1-a0-a1)\Phi(X \beta)}. 
The probability that an observation is classified as a 0 is 
\eqn{Pr(y = 0|X) = \alpha_1\Phi(X \beta) + (1-\alpha_0)(1-\Phi(X \beta)) = 1-\alpha_0-(1-\alpha_0-\alpha_1)\Phi(X \beta)}{Pr(y = 0|X) = a1\Phi(X \beta) + (1-a0)(1-\Phi(X \beta)) = 1-a0-(1-a0-a1)\Phi(X \beta)}.
The log-likelihood function for the probit model with misclassification is given by 
\deqn{lnL = \sum_i \{y_i ln(Pr(y_i=1|X_i)) + (1-y_i)ln(Pr(y_i=0)|X_i)) \} = }{lnL = \sum {y ln(Pr(y=1|X)) + (1-y)ln(Pr(y=0)|X)) } =}
\deqn{\sum_i  \{y_i ln(\alpha_0+(1-\alpha_0-\alpha_1)\Phi(X \beta)) + (1-y_i)ln(1-\alpha_0-(1-\alpha_0-\alpha_1)\Phi(X \beta)) \}. }{\sum  {y ln(a0+(1-a0-a1)\Phi(X \beta)) + (1-y)ln(1-a0-(1-a0-a1)\Phi(X \beta)) }. }


The log-likelihood function is maximized using the \emph{nlm} function.  
In practice, the model sometimes has difficulties converging because the maximization procedure attempts to set the 
misclassification probabilities outside the (0,1) interval.  
To avert this problem, the \emph{misclass} function estimates
\eqn{\alpha_0 = \Phi(\alpha_0^*) }{a0 = \Phi(a0*) }   and  \eqn{\alpha_1 = \Phi(\alpha_1^*) }{a1 = \Phi(a1*) }.  
The covariance matrix estimate is calculated using the \emph{hessian} option in \emph{nlm}.  
The vector \emph{estimate} contains the estimated values of \eqn{\beta},  
\eqn{\alpha_0^* = \Phi^{-1}(\alpha_0) = qnorm(\alpha_0)}{a0* = qnorm(a0)}, 
and \eqn{\alpha_1^* = \Phi^{-1}(\alpha_1) = qnorm(\alpha_1)}{a1* = qnorm(a1)}.
Similarly, \emph{stderr} and \emph{vmat} report the standard error estimates and the full covariance matrix 
for \eqn{(\beta \,\,\, \alpha_0^* \,\,\, \alpha_1^*)}{\beta,  a0*,  a1*}.
The estimated probabilities are reported in \emph{a0} and \emph{a1}.
	

By default, the starting values are obtained using a standard probit model 
with \eqn{\alpha_0}{a0} = 0 and \eqn{\alpha_1}{a1} = 0.  
The standard probit models are presented also.  
The starting values can be changed using the \emph{a0}, \emph{a1}, and \emph{bmat} options in \emph{misclass}.


}
\references{
Dye, Richard F. and Daniel P. McMillen, "Teardowns and Land Values in the Chicago Metropolitan Area," \emph{Journal of Urban Economics} 61 (2007), 45-64.

Hausman, J.A., Jason Arbrevaya, and F.M. Scott-Morton, "Misclassification of the Dependent Variable in a Discrete-Response Setting," \emph{Journal of Econometrics} 87 (1998), 239-269.

}

\examples{
set.seed(189)
n = 1000
x <- rnorm(n)
x <- sort(x)
y <- x*1 + rnorm(n, 0, sd(x)/2)
y <- ifelse(y>0,1,0)
e <- runif(n)
misy <- y
misy <- ifelse(e<.10&y==0,1,y)
misy <- ifelse(e>.90&y==1,0,misy)
table(y,misy)
fit <- misclass(misy~x)
}


\keyword{Misclassification}
\keyword{Discrete Choice Models}

