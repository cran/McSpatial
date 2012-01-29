misclass <- function(form,a0=0,a1=0,bmat=0,print.summary=TRUE,data=NULL) {

  probit <- glm(form,family=binomial(link="probit"),data=data)
  cat("Standard Probit Estimates","\n")
  print(summary(probit))
  xmat <- model.matrix(form,data=data)
  if (identical(bmat,0)) {bmat <- probit$coef}
  nk = length(bmat)
  y <- model.frame(form,data=data)[,1]

# ey <- a0 + (1-a0-a1)*pnorm(xb)

  a0 = ifelse(identical(a0,0), .001, a0)
  a1 = ifelse(identical(a1,0), .001, a1)

  if (!identical(a0,FALSE)&!identical(a1,FALSE)) {
    bstart <- c(bmat, qnorm(a0), qnorm(a1))
    logl <- function(b) {
      ey <- pnorm(b[nk+1]) + (1-pnorm(b[nk+1])-pnorm(b[nk+2]))*pnorm(xmat%*%b[1:nk])
      out <- -sum(ifelse(y==1,log(ey),log(1-ey)))
      return(out)
    }
  }
  if (identical(a0,FALSE)) {
    bstart <- c(bmat, qnorm(a1))
    logl <- function(b) {
      ey <- (1-pnorm(b[nk+1]))*pnorm(xmat%*%b[1:nk])
      out <- -sum(ifelse(y==1,log(ey),log(1-ey)))
      return(out)
    }
  }
  if (identical(a1,FALSE)) {
    bstart <- c(bmat, qnorm(a0))
    logl <- function(b) {
      ey <- pnorm(b[nk+1]) + (1-pnorm(b[nk+1]))*pnorm(xmat%*%b[1:nk])
      out <- -sum(ifelse(y==1,log(ey),log(1-ey)))
      return(out)
    }
  }


  nlfit <- nlm(logl,bstart,hessian=TRUE,iterlim=1000)
  vmat <- solve(nlfit$hessian)
  stderr <- sqrt(abs(diag(vmat)))
  tval <- nlfit$estimate/stderr
  pval <- 2*(1-pnorm(abs(tval)))
  
  nall <- length(nlfit$estimate)
  nk1 = nk+1
  amat <- pnorm(nlfit$estimate[nk1:nall])
  tmat <- tval[nk1:nall]
  smat <- abs(amat/tmat)
  
  if (!identical(a0,FALSE)&!identical(a1,FALSE)) {nvect <- c("a0","a1")}
  if (identical(a0,FALSE)) {nvect <- "a1"}
  if (identical(a1,FALSE)) {nvect <- "a0"}
  names(nlfit$estimate[nk1:nall]) <- nvect
  names(amat) <- nvect
  names(tmat) <- nvect
  names(smat) <- nvect

  cname <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  if (print.summary==TRUE) {
    outmat <- cbind(nlfit$estimate, stderr, tval, pval )
    rownames(outmat) <- c(names(probit$coef), paste("qnorm(", nvect, ")", sep=""))
    colnames(outmat) <- cname
    cat("nlm estimates","\n")
    cat(" ", "\n")
    print(round(outmat,5))
    cat(" ", "\n")

    cat("Estimated misclassification probabilities:", "\n")
    outmat <- cbind(amat, smat, tmat, 2*(1-pnorm(abs(tmat))))
    rownames(outmat) <- nvect
    colnames(outmat) <- cname
    print(round(outmat,5))
    cat(" ", "\n")
   
  }

  noa0 = identical(a0,FALSE)
  noa1 = identical(a1,FALSE)
  a0 <- ifelse(noa0==TRUE, 0, amat[1])
  a1 <- ifelse(noa1==TRUE, 0, ifelse(noa0==TRUE, amat[1], amat[2]))
  out1 <- list(a0,a1,nlfit$estimate,stderr,vmat,nlfit$iterations,-nlfit$minimum,nlfit$gradient)
  names(out1) <- c("a0","a1","estimate","stderr","vmat","iterations","minimum","gradient")
  return(out1)
}


