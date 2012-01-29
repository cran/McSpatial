mapoverlay <- function(shpfile,pointfile,shpvar,pointvar=NULL,func="sum") {
  library(maptools)
  longvar <- pointfile$longitude
  latvar  <- pointfile$latitude
  pointdata <- SpatialPoints(cbind(longvar,latvar))
  n = nrow(pointfile)
  obs <- overlay(pointdata,shpfile)
  pointout <- array(NA,dim=n)
  x <- data.frame(shpfile)[,shpvar]
  if (class(x)=="factor") {x <- levels(x)}
  pointout[!is.na(obs)] <- x[obs[!is.na(obs)]]

  shpout <- NULL
  shpdata <- data.frame(x)
  if (!identical(pointvar,NULL)) {
    y <- as.numeric(pointfile[,pointvar])
    shpout <- aggregate(y, by=list(pointout), func)
    colnames(shpout) <- c("x","y")
    shpdata <- merge(shpdata,shpout,by="x",all.x=TRUE)
    shpout <- shpdata[,"y"]
    
  }

  out <- list(pointout,shpout)
  names(out) <- c("pointout","shpout")
  return(out)
}

