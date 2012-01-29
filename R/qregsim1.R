qregsim1 <- function(formall, formx, bmat, xmin=NULL, xmax=NULL, graphx=TRUE, graphb=TRUE, graphy=TRUE, graphxb=TRUE, graphols=FALSE,
  histogram=FALSE, histfreq=FALSE, yname=NULL, xname=NULL, nsim=20000, bwadjust=1, data=NULL)  {

  xlabel = xname
  ylabel = yname

  xmat <- model.matrix(formx,data=data)
  xname <- colnames(model.frame(formx,data=data))
  x <- as.numeric(data[,xname])

  xmat <- model.frame(formall,data=data)
  y <- xmat[,1]
  yname <- colnames(xmat)[1]

  xmat <- model.matrix(formall,data=data)
  taumat <- as.numeric(rownames(bmat))
  ntau = length(taumat)
  n = nrow(xmat)

  if (identical(xlabel,NULL)) {xlabel = xname}
  if (identical(ylabel,NULL)) {ylabel = yname}

  if (graphx==TRUE) {
    if (histogram==FALSE) {
      fit <- density(x,adjust=bwadjust)
      plot(fit$x,fit$y,xlab=xlabel,ylab="Density",type="l")
    }
    if (histogram==TRUE) {
      hist(x,freq=histfreq,xlab=xlabel,main="")
    }
  }

  if (graphb==TRUE) {
    plot(taumat,bmat[,xname],xlab="Quantile",ylab="Coefficients",type="l",main="Quantile Regression Coefficients")
  }

  xobs <- sample(seq(1:n), nsim, replace=TRUE)
  bobs <- sample(seq(1:ntau), nsim, replace=TRUE)
  xhat <- xmat[xobs,]
  bhat <- bmat[bobs,]
  xbhat <- rowSums(as.matrix(xhat*bhat))

  if (graphols==TRUE) {
    fit <- lm(y~xmat[,-1])
    olshat <- fitted(fit)
    betaols <- fit$coef
    names(betaols) <- colnames(xmat)
  }

  ytarget <- NULL
  densy1 <- NULL
  densy2 <- NULL

  if (graphy==TRUE) {
    fit <- density(y,adjust=bwadjust)
    ytarget <- fit$x
    densy1 <- fit$y
    densy2 <- density(xbhat,adjust=bwadjust,from=min(ytarget),to=max(ytarget))$y
    ymin = min(densy1,densy2)
    ymax = max(densy1,densy2)
    if (graphols==TRUE){
      densy3 <- density(olshat,adjust=bwadjust,from=min(ytarget),to=max(ytarget))$y
      ymin = min(ymin,densy3)
      ymax = max(ymax,densy3)
    }
    plot(ytarget,densy1,type="l",xlab=ylabel,ylab="Density",main="Density Functions for Actual and Predicted Values",ylim=c(ymin,ymax))
    lines(ytarget,densy2,col="red")
    if (graphols==FALSE){legend("topright",c("Actual","Predicted"),col=c("black","red"),lwd=1) }
    if (graphols==TRUE){
      lines(ytarget,densy3,col="blue")
      legend("topright",c("Actual","Quantile","OLS"),col=c("black","red","blue"),lwd=1)
    }
  }

  if (identical(xmin,NULL)&identical(xmax,NULL)) {
    xbar = mean(x)
    xhat <- xmat[xobs,xname]
    bhat <- bmat[bobs,xname]
    xbmean <- xbhat - xhat*bhat + xbar*bhat
    dmin = min(xbmean,xbhat)
    dmax = max(xbmean,xbhat)
    if (graphols==TRUE){
      xbols <- olshat - betaols[xname]*x + xbar*betaols[xname]
      dmin = min(dmin,xbols)
      dmax = max(dmax,xbols)
      densx3 <- density(xbols,adjust=bwadjust,from=dmin,to=dmax)$y
    }
    densx1 <- density(xbhat,adjust=bwadjust,from=dmin,to=dmax)$y
    densx2 <- density(xbmean,adjust=bwadjust,from=dmin,to=dmax)$y
    xtarget <- seq(dmin,dmax,length=512)
    ymin = ifelse(graphols==TRUE,min(densx1,densx2,densx3),min(densx1,densx2))
    ymax = ifelse(graphols==TRUE,max(densx1,densx2,densx3),max(densx1,densx2))
    if (graphxb==TRUE){
      plot(xtarget,densx1,type="l",xlab=ylabel,ylab="Density",ylim=c(ymin,ymax),
        main="Counterfactual Density Functions")
      lines(xtarget,densx2,col="red")
      if (graphols==FALSE){legend("topright",c("Actual","Quantile"),col=c("black","red"),lwd=1) }
      if (graphols==TRUE){
        lines(xtarget,densx3,col="blue")
        legend("topright",c("Actual","Quantile","OLS"),col=c("black","red","blue"),lwd=1)
      }
    }
  }
  
  if (!identical(xmin,NULL)|!identical(xmax,NULL)) {
    xhat <- xmat[xobs,xname]
    bhat <- bmat[bobs,xname]
    if (!identical(xmin,NULL)&!identical(xmax,NULL)) {
      x1 <- array(xmin,dim=nsim)
      x2 <- array(xmax,dim=nsim)
      x1ols <- array(xmin,dim=n)
      x2ols <- array(xmax,dim=n)
    }
    if (identical(xmin,NULL)|identical(xmax,NULL)) {
      xdelta <- ifelse(identical(xmin,NULL),xmax,xmin)
      x1 <- xhat
      x2 <- xhat + xdelta
      x1ols <- x
      x2ols <- x + xdelta
    }

    xb1 <- xbhat - xhat*bhat + x1*bhat
    xb2 <- xbhat - xhat*bhat + x2*bhat
    dmin = min(xb1,xb2)
    dmax = max(xb1,xb2)
    if (graphols==TRUE){
      xb3 <- olshat - betaols[xname]*x + x1ols*betaols[xname]
      xb4 <- olshat - betaols[xname]*x + x2ols*betaols[xname]
      dmin = min(xb1,xb2,xb3,xb4)
      dmax = max(xb1,xb2,xb3,xb4)
      densx3 <- density(xb3,adjust=bwadjust,from=dmin,to=dmax)$y
      densx4 <- density(xb4,adjust=bwadjust,from=dmin,to=dmax)$y
    }
    densx1 <- density(xb1,adjust=bwadjust,from=dmin,to=dmax)$y
    densx2 <- density(xb2,adjust=bwadjust,from=dmin,to=dmax)$y
    xtarget <- seq(dmin,dmax,length=512)

    minlab = as.character(xmin)
    maxlab = as.character(xmax)
    if (identical(xmin,NULL)|identical(xmax,NULL)) {
      minlab = "Actual"
      maxlab = paste("Actual + ",as.character(xdelta))
    }
    if (graphxb==TRUE){
      if (graphols==FALSE){
        ymin = min(densx1,densx2)
        ymax = max(densx1,densx2)
        plot(xtarget,densx1,type="l",xlab=ylabel,ylab="Density",ylim=c(ymin,ymax))
        lines(xtarget,densx2,col="red")
        legend("topright",c(minlab,maxlab),col=c("black","red"),lwd=1)
      }
      if (graphols==TRUE){
        ymin = min(densx1,densx2,densx3,densx4)
        ymax = max(densx1,densx2,densx3,densx4)
        plot(xtarget,densx1,type="l",xlab=ylabel,ylab="Density",ylim=c(ymin,ymax))
        lines(xtarget,densx2,lty="dashed")
        lines(xtarget,densx3,col="red")
        lines(xtarget,densx4,col="red",lty="dashed")
        legend("topright",c(minlab,maxlab),lty=c("solid","dashed"),lwd=1)
        legend("topleft",c("Quantile","OLS"),col=c("black","red"),lwd=1)
      }
    }

  }
  
  out <- list(xtarget,densx1,densx2,ytarget,densy1,densy2)
  names(out) <- c("xtarget","densx1","densx2","ytarget","densy1","densy2")
  return(out)

}

