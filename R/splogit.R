splogit <- function(form,inst=NULL,winst=NULL,wmat=NULL,shpfile=NULL, data=NULL, silent=FALSE) {
   
  library(car)
  if (length(wmat)==0) {
    library(spdep)
    neighbors <- poly2nb(shpfile,queen=TRUE)
    wmat <- nb2mat(neighbors,zero.policy=TRUE)
  }
  xmat <- model.frame(form,data=data)
  xnames <- colnames(xmat)

  y <- xmat[,1]
  xmat <- model.matrix(form,data=data)
  if (identical(inst,NULL)&identical(winst,NULL)) {zmat <- cbind(xmat, wmat%*%xmat[,-1])}
  if (identical(inst,NULL)&!identical(winst,NULL)) {zmat <- cbind(xmat, wmat%*%(model.matrix(winst,data=data)[,-1])) }
  if (!identical(inst,NULL)&identical(winst,NULL)) {zmat <- model.matrix(inst,data=data)}
  if (!identical(inst,NULL)&!identical(winst,NULL)) {zmat <- cbind(model.matrix(inst,data=data), wmat%*%(model.matrix(winst,data=data)[,-1])) }
  if (!identical(inst,NULL)&identical(winst,NULL)) {
   cat("Warning:  list provided for inst but not winst", "\n")
   cat("inst list should include variables that are omitted from orginal explanatory variable list", "\n")
   cat("\n")
  }

  if (silent==FALSE) {cat("STANDARD LOGIT ESTIMATES","\n")}
  logit <- glm(form,family=binomial(link="logit"),data=data)
  if (silent==FALSE) {print(summary(logit))}
  xb <- xmat%*%logit$coef
  p <- exp(xb)/(1+exp(xb))
  grad <- as.vector(p*(1-p))
  gmat <- grad*xmat
  u <- y-p + gmat%*%logit$coef
  wxb <- grad*(wmat%*%xb)
  gmat <- cbind(gmat,wxb)

  gmat <- zmat%*%(solve(crossprod(zmat))%*%(t(zmat)%*%gmat))
  fit <- lm(u~gmat+0)
  v <- diag(hccm(fit))
  summat <- cbind(fit$coef,sqrt(v),fit$coef/sqrt(v),2*(1-pnorm(abs(fit$coef)/sqrt(v) )) )
  rownames(summat) <- c(xnames, "WXB") 
  colnames(summat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")

  if (silent==FALSE) {
    cat("LINEARIZED GMM LOGIT ESTIMATES","\n")
    print(round(summat,5))
  }

  out <- list(fit$coef,sqrt(v))
  names(out) <- c("coef","se")
  return(out)

}

