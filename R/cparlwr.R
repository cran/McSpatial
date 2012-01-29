cparlwr <- function(form,nonpar,window=.25,bandwidth=0,kern="tcub",distance="Mahal",alldata=FALSE,data=NULL) {
  library(locfit)
  library(akima)

  mat <- model.frame(form,data=data)
  y <- mat[,1]
  n = length(y)
  xmat <- model.matrix(form,data=data)
  zmat <- model.frame(nonpar,data=data)
  nz = ncol(zmat)
  nk = ncol(xmat)

  xname <- colnames(xmat)
  nocons = xname[1]!="(Intercept)"

  if (nz==1) {vzmat <- var(zmat) }
  if (nz==2) {
    vzmat <- cov(zmat) 
    if (distance=="Euclid"|distance=="E") {vzmat <- diag(diag(vzmat)) }
  }

  if (kern=="rect")  { wgt <- function(psi) {ifelse(abs(psi)>=0,1,0) } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { 1-psi^2 } }
  if (kern=="bisq")  { wgt <- function(psi) { (1-psi^2)^2 } }
  if (kern=="tcub")  { wgt <- function(psi) { (1 - abs(psi)^3)^3 } }
  if (kern=="trwt")  { wgt <- function(psi) { (1 - psi^2)^3 } }
  if (kern=="gauss") { wgt <- function(psi) { exp(-((2.5*psi)^2)/2) } }

  if (bandwidth>0) {window = 0}

  if (nz==1&window>0)    {fit <- locfit(~lp(zmat[,1],nn=window,deg=1),kern=kern) }
  if (nz==2&window>0)    {fit <- locfit(~lp(zmat[,1],zmat[,2],nn=window,deg=1),kern=kern) }
  if (nz==1&bandwidth>0) {fit <- locfit(~lp(zmat[,1],h=2*bandwidth,deg=1),kern=kern) }
  if (nz==2&bandwidth>0) {fit <- locfit(~lp(zmat[,1],zmat[,2],h=2*bandwidth,deg=1),kern=kern) }
 
  if (alldata==FALSE) {
    zev <- lfeval(fit)$xev
    nt = length(zev)/nz
    target <- t(array(zev,dim=c(nz,nt)))
    obs <- array(0,dim=nt)
    for (i in seq(1:nt)) {
      dist <- sqrt(mahalanobis(zmat, target[i,], vzmat))
      obs[i] <- which.min(dist)
    }
    colnames(target) <- colnames(zmat)
    obs <- sort(unique(c(obs,chull(zmat))))
    nt = length(obs)
    xvect <- as.matrix(xmat[obs,])
    target <- as.matrix(zmat[obs,])
  }
  if (alldata==TRUE) {
    target <- as.matrix(zmat)
    obs <- seq(1:n)
    nt = n
    xvect <- xmat
  }

  if (distance=="Latlong"|distance=="L") {
    tvect <- attr(terms(nonpar),"term.labels")
    if (substr(tvect[1],1,2)=="la"|substr(tvect[1],1,2)=="La"|substr(tvect[1],1,2)=="LA") {
      la  <- 2*pi*zmat[,1]/360 
      la1 <- 2*pi*target[,1]/360
    }
    if (substr(tvect[2],1,2)=="la"|substr(tvect[2],1,2)=="La"|substr(tvect[2],1,2)=="LA") {
      la  <- 2*pi*zmat[,2]/360 
      la1 <- 2*pi*target[,2]/360 
    }
    if (substr(tvect[1],1,2)=="lo"|substr(tvect[1],1,2)=="Lo"|substr(tvect[1],1,2)=="LO") {
      lo  <- 2*pi*zmat[,1]/360 
      lo1 <- 2*pi*target[,1]/360
    }
    if (substr(tvect[2],1,2)=="lo"|substr(tvect[2],1,2)=="Lo"|substr(tvect[2],1,2)=="LO") {
      lo  <- 2*pi*zmat[,2]/360 
      lo1 <- 2*pi*target[,2]/360
    }
  }

  xcoef.target <- array(0,dim=c(nt,nk))
  df1target <- array(0,dim=nt)
  df2target <- array(0,dim=nt)
  xcoef.target.se <- array(0,dim=c(nt,nk))

  for (i in seq(1:nt)) {
    if (distance!="Latlong"&distance!="L") {dist <- sqrt(mahalanobis(zmat, target[i,], vzmat)) }
    if (distance=="Latlong"|distance=="L") {
      dist <- pmin(sin(la)*sin(la1[i]) + cos(la)*cos(la1[i])*cos(lo1[i]-lo),  1)
      dist <- acos(dist)*3958
    }
    if (window>0) {h = quantile(dist,window) }
    if (bandwidth>0) {h = bandwidth}
    samp <- dist<=h
    if (kern=="gauss") {samp <- dist<=max(dist)}
    
    xmat1 <- as.matrix(xmat[samp,]) 

    k <- wgt(dist[samp]/h)
    xmat2 <- k*xmat1
    xx <- solve(crossprod(xmat1,xmat2))
    xmat1 <- xx%*%t(xmat2)
    xcoef.target[i,] <- xmat1%*%y[samp]
    vmat <- tcrossprod(xmat1)

    xmat1 <- xvect[i,]%*%xmat1
    df2target[i] = sum(xmat1^2)
    dist[samp] <- xmat1
    dist[!samp] <- 0
    df1target[i] = dist[obs[i]]

    xcoef.target.se[i,] <- sqrt(diag(vmat))
  }

  xcoef <- array(0,dim=c(n,nk))
  xcoef.se <- array(0,dim=c(n,nk))

  if (alldata==FALSE) {
    if (nz==1) {
      for (j in seq(1:nk)) {
        hat <- aspline(target,xcoef.target[,j],zmat[,1])
        xcoef[,j] <- hat$y
        hat <- aspline(target,xcoef.target.se[,j],zmat[,1])
        xcoef.se[,j] <- hat$y
       }
      hat <- aspline(target,df1target,zmat[,1])
      infl <- hat$y
      df1 = sum(hat$y)
      hat <- aspline(target,df2target,zmat[,1])
      df2 = sum(hat$y)
    }
    if (nz==2) {
      for (j in seq(1:nk)) {
        hat <- interpp(target[,1],target[,2],xcoef.target[,j], zmat[,1],zmat[,2],duplicate="mean")
        xcoef[,j] <- hat$z
        hat <- interpp(target[,1],target[,2],xcoef.target.se[,j],zmat[,1],zmat[,2],duplicate="mean")
        xcoef.se[,j] <- hat$z
       }
      hat <- interpp(target[,1],target[,2],df1target, zmat[,1],zmat[,2],duplicate="mean")
      infl <- hat$z
      df1 = sum(hat$z)
      hat <- interpp(target[,1],target[,2],df2target, zmat[,1],zmat[,2],duplicate="mean")
      infl <- hat$z
      df2 = sum(hat$z)
    }
  }

  if (alldata==TRUE) {
    xcoef <- xcoef.target
    xcoef.se <- xcoef.target.se
    infl <- df1target
    df1 = sum(infl)
    df2 = sum(df2target)
  }

#  yhat <- diag(tcrossprod(xmat,xcoef))   
  yhat <- xmat[,1]*xcoef[,1]
  if (nk>1) {
    for (j in seq(2,nk)) {
      yhat <- yhat + xmat[,j]*xcoef[,j]
    }
  }

  rss = sum((y-yhat)^2)
  sig2 = rss/(n-2*df1 + df2)
  cv = mean(((y-yhat)/(1-infl))^2)
  gcv = n*rss/((n-df1)^2)

  xcoef.target.se <- sqrt(sig2)*xcoef.target.se
  xcoef.se <- sqrt(sig2)*xcoef.se

  colnames(xcoef.target) <- xname
  colnames(xcoef.target.se) <- xname
  colnames(xcoef) <- xname
  colnames(xcoef.se) <- xname

  out <- list(target,xcoef.target,xcoef.target.se,yhat,xcoef,xcoef.se,df1,df2,sig2,cv,gcv,infl)
  names(out) <- c("target","xcoef.target","xcoef.target.se","yhat","xcoef","xcoef.se","df1","df2","sig2","cv","gcv","infl")
  return(out)   
  detach(data)
}


