mappoints <- function(shpfile,longvar,latvar,pointvar,sampvar=NULL,nclass=8,col="RdBu",pointsize=.6,colreverse=TRUE,
  xlimfact=c(1,1),ylimfact=c(1,1),legendloc=NULL,legendsize=.75,legenddigits=4,title=NULL) {
  library(RColorBrewer)
  library(classInt)
  if (!identical(sampvar,NULL)) {shpfile <- shpfile[sampvar==TRUE,]}
  box <- bbox(shpfile)
  box[1,] <- box[1,]*xlimfact
  box[2,] <- box[2,]*ylimfact

  plot(shpfile,main=title,xlim=c(box[1,]),ylim=c(box[2,]))
  pcolor <- brewer.pal(nclass,col)
  if (colreverse==TRUE){pcolor <- pcolor[nclass:1] }
  class <- classIntervals(pointvar,nclass,style="quantile")
  class$brks <- signif(class$brks,digits=legenddigits)
  colcode <- findColours(class, pcolor)
  points(longvar,latvar,col=findColours(class,pcolor),cex=pointsize)
  if (!identical(legendloc,NULL)) {
    legend(legendloc,legend=names(attr(colcode,"table")),fill=attr(colcode,"palette"),cex=legendsize)
  }
  title(main=title)
}

