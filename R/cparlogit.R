cparlogit <- function(form,nonpar,window=.25,bandwidth=0,kern="tcub",distance="Mahal",alldata=FALSE,data=NULL) {

  library(locfit)
  library(akima)

  mat <- model.frame(form,data=data)
  y <- as.numeric(mat[,1])
  n = length(y)
  xmat <- model.matrix(form,data=data)
  zmat <- model.frame(nonpar,data=data)
  nz = ncol(zmat)
  nk = ncol(xmat)

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
  xcoef.target.se <- array(0,dim=c(nt,nk))
  fit <- glm(y~xmat+0,family=binomial(link="logit"))
  bstart <- fit$coef

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
    y1 <- y[samp]
    w <- wgt(dist[samp]/h)

    nlfit <- glm(y1~xmat1+0, family=binomial(link="logit"), weights=w, start=bstart)
    xcoef.target[i,] <- nlfit$coef

    xb <- as.numeric(xmat1%*%xcoef.target[i,])
    p <- exp(xb)/(1+exp(xb))
    u <- y1-p
    gmat1 <- crossprod(as.matrix(w*u*(as.data.frame(xmat1))))
    d <- sqrt( p*(1-p) )
    gmat2 <- solve( crossprod( as.matrix(d*(as.data.frame(xmat1))) ))
    vmat <- gmat2%*%gmat1%*%gmat2
    xcoef.target.se[i,] <- sqrt(diag(vmat))
 }

  hat1 <- function(x) {
    if (nz==1) {hat <- aspline(target,x,zmat[,1])$y}
    if (nz==2) {hat <- interpp(target[,1],target[,2], x, zmat[,1],zmat[,2],duplicate="mean")$z }
    return(hat)
  }
  
  xcoef <- array(0,dim=c(n,nk))
  xcoef.se <- array(0,dim=c(n,nk))
  for (j in seq(1:nk)) {
    xcoef[,j] <- hat1(xcoef.target[,j])
    xcoef.se[,j] <- hat1(xcoef.target.se[,j])
  }
  yhat <- array(0,dim=n)
  xb <- rowSums(as.data.frame(xmat)*xcoef)
  p <- exp(xb)/(1+exp(xb))
  lnl = sum(ifelse(y==1,log(p),log(1-p)))
    
  out <- list(target,xcoef.target,xcoef.target.se,xcoef,xcoef.se,p,lnl)
  names(out) <- c("target","xcoef.target","xcoef.target.se","xcoef","xcoef.se","p","lnl")
  return(out)
}
