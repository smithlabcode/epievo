#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Rscript plot.R <paths.txt> <paths.pdf> <freq.txt> <freq.pdf>.n", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  paths <- args[1]	
  pathplot <- args[2]
  freqs <- args[3]
  freqplot <- args[4]
}


x <- read.table(paths, as.is=T)
tot <- as.numeric(unlist(strsplit(x$V3[1], ",")))[2]

pdf (pathplot, width=120, height=3, pointsize=8 )
colors <- c("cadetblue", "black")
plot(x=c(0,nrow(x)), y=c(0,tot), pch="", xlab= "Position", ylab="Time")
rect(xleft=0, ybottom=0, xright=nrow(x), ytop=tot)
for (i in 1:nrow(x)) {
xval <- x$V1[i]
s <- x$V2[i]
yvals <- as.numeric(unlist(strsplit(x$V3[i], ",")))

for (k in 1:(length(yvals)-1)) {
  if (s==0)
    rect(xleft=xval, ybottom=yvals[k], xright=xval+1, ytop=yvals[k+1], 
         col = colors[s+1], border = NA)
  s = 1-s
}
}
dev.off()


y <- read.table(freqs)
pdf(freqplot, width=4, height=4, pointsize=8)
par(mfrow=c(2,2), mar=c(5,4,1,1))
plot(y$V2, y$V3, xlab="Time", ylab="Freq 00")
plot(y$V2, y$V6, xlab="Time", ylab="Freq 11")
plot(y$V2, y$V4, xlab="Time", ylab="Freq 01")
plot(y$V2, y$V5, xlab="Time", ylab="Freq 10")
dev.off()
