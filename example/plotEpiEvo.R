#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("plotEpiEvo.R <input.states> <n.cols> <output.pdf>", call.=FALSE)
}

library(raster)

color.breakpoints <- c(0.0,0.5,1.0)
colors <- c("red","blue")

X <- readLines(args[1])
the.width <- as.numeric(args[2])

m <- matrix(nrow=length(X), ncol=the.width)

for (i in 1:length(X)) {
  m[i, ] <- as.numeric(unlist(strsplit(substr(X[i], 1, the.width), '')))
}

pdf(args[3])
plot(raster(m),breaks=color.breakpoints,col=colors,legend=FALSE)
dev.off()
