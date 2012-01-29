spprobitml <- function(form,wmat,blockid=NULL,stdprobit=TRUE,data=NULL) {
  library(spdep)
  xmat <- model.frame(form,data=data)
  y <- xmat[,1]
  xmat <- model.matrix(form,data=data)
  n = length(y)

  if (identical(blockid,NULL)) {blockid <- array(1,dim=n)}
  blockmat <- table(blockid)
  totblock <- as.numeric(blockmat)
  blocknum <- as.numeric(names(blockmat))
  nblock = length(totblock)

  if (stdprobit==TRUE) {
    fit <- glm(form,family=binomial(link="probit"),data=data)
    cat("Standard Probit Estimates","\n")
    print(summary(fit))
  }

  makevar <- function(rho) {
    xstar <- xmat
    v <- array(1,dim=n)
    for (i in seq(1,nblock)) {
      sampvar <- blockid==blocknum[i]
      vmat <- solve(diag(totblock[i]) - rho*wmat[sampvar==TRUE, sampvar==TRUE])
      xstar[sampvar==TRUE,] <- vmat%*%xmat[sampvar==TRUE,]
      vmat <- tcrossprod(vmat)
      v[sampvar==TRUE] <- sqrt(diag(vmat))
    }
    out <- list(xstar,v) 
    names(out) <- c("xstar","v")
    return(out)
  }

  probitrho <- function(rho) {
    fit <- makevar(rho)
    xstar <- as.matrix(as.data.frame(fit$xstar)/fit$v)
    fit <- glm(y~xstar+0,family=binomial(link="probit"))
    xb <- as.numeric(as.matrix(xstar)%*%fit$coef)
    lvar <- sum(ifelse(y==1, log(pnorm(xb)), log(1-pnorm(xb))))
    out <- list(fit$coef,lvar)
    names(out) <- c("coef","logl")
    return(out)
  }

  logl <- function(rho) {-probitrho(rho)$logl}  
  rho = optimize(logl,lower=-.99,upper=.99)$minimum
  fit <- probitrho(rho)
  bvect <- c(fit$coef, rho)
  names(bvect) <- c(colnames(xmat),"rho")
   
  nk = length(bvect)-1
  rho = bvect[nk+1]
  b <- bvect[1:nk]
  fit <- makevar(rho)
  xb <- as.numeric(fit$xstar%*%b)/fit$v
  p <- pnorm(xb)
  g <- (dnorm(xb)^2)/(p*(1-p))
  gmat <- as.matrix((sqrt(g)/fit$v)*data.frame(fit$xstar))
  vmat1 <- solve(crossprod(gmat))

# Use numeric derivatives to calculate unconditional standard errors
  logl <- function(rho) {
    fit <- makevar(rho)
    sv <- fit$v
    xb <- as.numeric(as.matrix(fit$xstar)%*%b)/sv
    lvar <- ifelse(y==1, log(pnorm(xb)), log(pnorm(-xb)))
    return(lvar)
  }
  lvar <- sum(logl(rho))
  g <- dnorm(xb)*(y-p)/(p*(1-p))
  gmat <- as.matrix(data.frame(fit$xstar)*(g/fit$v))
  g1 <- logl(rho+.001)
  g0 <- logl(rho-.001)
  g <- (g1-g0)/.002
  vmat2 <- solve(crossprod(cbind(gmat,g)))

  semat1 <- sqrt(diag(vmat1))
  semat2 <- sqrt(diag(vmat2))
  cat("Conditional on rho","\n")
  cat("rho = ", rho, "\n")
  outmat <- cbind(bvect[1:nk], semat1, bvect[1:nk]/semat1, 2*(1-pnorm(abs(bvect[1:nk])/semat1)) )
  colnames(outmat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  rownames(outmat) <- colnames(xmat)
  print(outmat) 

  cat("Unconditional Standard Errors","\n")
  outmat <- cbind(bvect, semat2, bvect/semat2, 2*(1-pnorm(abs(bvect)/semat2)) )
  colnames(outmat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  rownames(outmat) <- c(colnames(xmat),"rho")
  print(outmat)

  out <- list(bvect,lvar,vmat1,vmat2)
  names(out) <- c("coef","logl","vmat1","vmat2")
  return(out)

}



