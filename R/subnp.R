subnp <- function(ydens,long,lat,window=.5,pval=.10) {
  library(locfit)
  data <- data.frame(ydens,long,lat)
  names(data) <- c("ydens","long","lat")
  fit <- locfit(ydens~lp(long,lat,nn=window,deg=1),kern="tcub",ev=dat(cv=FALSE),data=data)
  mat <- predict(fit,se.fit=TRUE,band="pred")
  yhat <- mat$fit
  sehat <- mat$se.fit
  upper <- yhat - qnorm(pval/2)*sehat
  subobs <- ifelse(ydens>upper,1,0)
  
  cat("Number of tracts identified as part of subcenters:  ",sum(subobs),"\n")
  out <- list(subobs)
  names(out) <- c("subobs")
  return(out)
}

