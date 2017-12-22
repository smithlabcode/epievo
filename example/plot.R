#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Rscript plot.R <paths.txt> <paths.pdf>.n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  paths <- args[1]	
  pathplot <- args[2]
}


dat <- read.table(paths, as.is=T)

colors <- c("cadetblue", "black")
pdf (pathplot, width=120, height=3, pointsize=8 )
for (node in unique(dat$V1)) {
  x <- dat[dat$V1==node,-1]
  tot <- as.numeric(unlist(strsplit(x[1, 3], ",")))[2]
  plot(x=c(0,nrow(x)), y=c(0,tot), pch="", xlab= "Position", ylab="Time", main=paste("Node", node))
  rect(xleft=0, ybottom=0, xright=nrow(x), ytop=tot)
  for (i in 1:nrow(x)) {
    xval <- x[i, 1]
    s <- x[i, 2]
    yvals <- as.numeric(unlist(strsplit(x[i, 3], ",")))

    for (k in 1:(length(yvals)-1)) {
      if (s==0)
        rect(xleft=xval, ybottom=yvals[k], xright=xval+1, ytop=yvals[k+1], 
             col = colors[s+1], border = NA)
      s = 1-s
    }
  }
}

dev.off()

#y <- read.table(freqs)
#pdf(freqplot, width=4, height=4, pointsize=8)
#par(mfrow=c(2,2), mar=c(5,4,1,1))
#plot(y$V2, y$V3, xlab="Time", ylab="Freq 00")
#plot(y$V2, y$V6, xlab="Time", ylab="Freq 11")
#plot(y$V2, y$V4, xlab="Time", ylab="Freq 01")
#plot(y$V2, y$V5, xlab="Time", ylab="Freq 10")
#dev.off()
