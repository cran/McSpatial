subgs <-function(shpfile,dens,emp,mind=10,totemp=10000, wmat=0) {
  library(spdep)
  if (identical(wmat,0)) {
    library(spdep)
    neighbors <- poly2nb(shpfile,queen=TRUE)
    wmat <- nb2mat(neighbors,style="B",zero.policy=TRUE)
  }
  dens <- ifelse(is.na(dens),0,dens)
  obs <- seq(1:length(dens))
  densobs <- obs[dens>mind]
  wmat <- wmat[densobs,densobs]
  n = nrow(wmat)
  amat <- matrix(0,nrow=n,ncol=n)
  amat[row(amat)==col(amat)] <- 1
  bmat <- wmat
  wmat1 <- wmat
  newnum = sum(bmat)
  cnt = 1
  while (newnum>0) {
    amat <- amat+bmat
    wmat2 <- wmat1%*%wmat
    bmat <- ifelse(wmat2>0&amat==0,1,0)
    wmat1 <- wmat2
    newnum = sum(bmat)
    cnt = cnt+1
  }
  emat <- emp[dens>mind]
  tmat <- amat%*%emat
  obsmat <- densobs[tmat>totemp]

  subemp <- array(0,dim=length(dens))
  subemp[obsmat] <- tmat[tmat>totemp]
  subobs <- ifelse(subemp>0,1,0)

  tab <- tabulate(factor(subemp))
  numsub = sum(tab>0)-1

  cat("Number of Subcenters = ",numsub,"\n")
  cat("Total Employment and Number of Tracts in each Subcenter","\n")
  print(table(subemp))

  out <- list(subemp,subobs)
  names(out) <- c("subemp","subobs")
  return(out)
}

