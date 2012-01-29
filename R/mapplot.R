mapplot <- function(shpfile,varname,brks=8,col="RdBu",sampvar=NULL,print=TRUE,title=NULL) {
  library(RColorBrewer)
  library(maptools)
  if (!identical(sampvar,NULL)) {shpfile <- shpfile[sampvar==TRUE,]}
  x <- data.frame(shpfile[,which(names(shpfile)==varname)])
  nbrks = ifelse(length(brks)==1, brks, length(brks)-1)
  if (length(brks)==1) {brks <- seq(min(x[,1],na.rm=TRUE),max(x[,1],na.rm=TRUE),length=(nbrks+1)) }
  if (print==TRUE) {cat("breaks = ",brks,"\n") }
  spplot(shpfile,varname,at=brks,col.regions=rev(brewer.pal(nbrks,col)),main=title)
}

