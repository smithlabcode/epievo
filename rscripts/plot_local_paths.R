#!/usr/local/bin/Rscript
/* Copyright (C) 2019 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith, Jianghan Qu and Xiaojing Ji
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Rscript plot.R <prefix> <nodeidx 0-based> <position_lo> <position_hi> <pdf>" , call.=FALSE)
} else if (length(args)==5) {
  # default output file
  paths <- paste(args[1], ".paths", sep="")
  which_node <- as.numeric(args[2])
  xrange <- as.numeric(c(args[3], args[4]))
  pathplot <- args[5]
}

print(paste("Reading from file", paths))
dat <- scan(paths, what="", sep="\n")
name_lines <- grep("NODE", dat)
n_nodes <- length(name_lines)

if (which_node >= n_nodes)
  stop("nodeidx outof range", call.=FALSE)

begin <- name_lines[which_node+1] + 1
if (which_node >=1 && which_node < n_nodes - 1 ) {
  end <- name_lines[which_node+2] - 1
} else {
  end <- length(dat)
}
branch_path <- dat[begin:end]

records <- strsplit(branch_path, "\t")

nsites  <- length(records)
pdf(pathplot, width=7, height=3, pointsize=8 )
colors <- c("white", "black")
#colors <- c("black", "white")
# assuming first site has no state change
tot <- as.numeric(records[[1]][3])
plot(x=xrange, y=c(0,tot), ylim=c(tot, 0), pch="", xlab= "Position", ylab="Time",
     main=dat[name_lines[which_node+1]])
rect(xleft=xrange[1], ybottom=0, xright=xrange[2], ytop=tot)
for (i in xrange[1]:xrange[2]) {
  xval <- i
  s <- as.numeric(records[[i]][2])
  #print(s)
  yvals <- c(0, as.numeric(records[[i]][-c(1:3)]), tot)
  for (k in 1:(length(yvals) - 1)) {
    if (s==1)
      rect(xleft=xval, ybottom=yvals[k], xright=xval+1, ytop=yvals[k+1],
           col = colors[s+1], border = NA)
    s = 1-s
  }
}
dev.off()




### old code wrapped in unused function
.f = function(){
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Rscript plot.R <prefix>.n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  paths <- paste(args[1], "_path.txt", sep="")
  pathplot <- paste(args[1], "_path.pdf", sep="")
  hmrplot <- paste(args[1], "_hmrsizes.pdf", sep="")
}

dat <- read.table(paths, as.is=T)
names(dat) <- c("node", "pos", "s0", "jumps")

pdf(pathplot, width=120, height=3, pointsize=8 )
colors <- c("cadetblue", "black")
for (node in unique(dat$node)) {
  x <- dat[dat[,"node"]==node,-1]
  # assuming first site has no state change
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

get_states <- function(dat, node, time) {
  subdat <- dat[which(dat$node == node), ]
  states <- rep(0, nrow(subdat))
  tot <- as.numeric(unlist(strsplit(subdat$jumps[1], ",")))[2]
  for (pos in 1:nrow(subdat)) {
    endpoints <- as.numeric(unlist(strsplit(subdat$jumps[pos], ",")))
    jumps <- endpoints[-c(1, length(endpoints))]
    states[pos] <- (subdat$s0[pos] + sum(jumps < time)) %%2
  }
  return (states)
}

hmr_sizes <- function(states) {
  x <- numeric(0)
  in_hmr <- F
  y <- 0
  for (i in 1:length(states)) {
    if (states[i] == 0) {
      y <- y + 1
    } else {
      if (y > 0) {x <- c(x, y)}
      y <- 0
    }
  }
  if (y > 0) { x <- c(x, y)}
  return (x)
}

get_hmr_sizes <- function(dat, node, time) {
  states <- get_states(dat, node, time)
  sizes <- hmr_sizes(states)
  return (sizes)
}

get_hmr_frac <- function(dat, node, time) {
  states <- get_states(dat, node, time)
  frac <- 1- sum(states)/length(states)
  return (frac)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(ggplot2)

pdf(hmrplot, width=7, height=4, pointsize=8)
par(mfrow =c(1,2))
for (node in unique(dat$node)) {
  subdat <- dat[which(dat$node == node), ]
  tot <- as.numeric(unlist(strsplit(subdat$jumps[1], ",")))[2]
  time_step <- tot/10
  time_points <- seq(time_step, tot, by=time_step)
  size <- data.frame(time=numeric(0), hmrsize=numeric(0))
  frac <- data.frame(time=numeric(0), hmrfrac=numeric(0))
  for (t in time_points) {
    hs <- get_hmr_sizes(dat, node, t)
    size <- rbind(size, data.frame(time=rep(t, length(hs)), hmrsize=hs))
    frac <- rbind(frac,data.frame(time=t, hmrfrac=get_hmr_frac(dat, node, t)))
  }
  size$time <- as.factor(size$time)

  p1 <- ggplot(size, aes(x=time, y=hmrsize)) +  geom_violin() +
       geom_boxplot(width=0.1, outlier.size=0.2) + theme_bw() +
       ggtitle(paste("branch", node)) + ylab("HMR size") +
       theme(panel.border = element_rect(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line = element_line(colour = "black"),
             plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(angle = 90, hjust = 0))
  p2 <- ggplot() + geom_line(aes(y = hmrfrac, x = time),data = frac, stat="identity") +
      ggtitle(paste("branch", node)) + ylab("HMR fraction") +theme_bw() +
      theme(panel.border = element_rect(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line = element_line(colour = "black"),
             plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(angle = 90, hjust = 0))
 multiplot(p1, p2, cols=2)

}
dev.off()
}
