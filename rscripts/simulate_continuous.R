args = commandArgs(trailingOnly=TRUE)
dat <- read.table(args[1])
dat <- dat[,order(-1:-ncol(dat))]
#dat <- log(dat[,order(-1:-ncol(dat))]+1)

n.sample <- 100
time <- 1

pdf(file=args[2],
    width=7, height=3, pointsize=8)

layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), 
       heights=c(1, 5, 1))

root <- dat[, ncol(dat)]
leaf <- dat[, 1]

# root
par(mar=c(0.5, 4, 3.5, 4))
print(dim(dat))
image(as.matrix(root), col=gray((n.sample:0)/n.sample),
      useRaster=TRUE, axes=FALSE, xlab="", ylab="")
box()

# paths
par(mar=c(0, 4, 0, 4))
image(as.matrix(dat), col=gray((n.sample:0)/n.sample),
      useRaster=TRUE, axes=FALSE, xlab="", ylab="Time")
axis(2, at = seq(0, 1, by = 1/10),
     labels = seq(0, time, by = time/10))
box()

# leaf
par(mar=c(3.5, 4, 0.5, 4))
image(as.matrix(leaf), col=gray((n.sample:0)/n.sample),
      useRaster=TRUE, axes=FALSE, xlab="Sites", ylab="", mgp=c(2.5, 0.1, 0))
axis(1, at = seq(0, 1, by = 1/10),
     labels = seq(0, nrow(dat), by = nrow(dat)/10))
box()

dev.off()
