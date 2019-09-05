#!/usr/local/bin/Rscript
# Copyright (C) 2019 University of Southern California
# Xiaojing Ji, Jianghan Qu and Andrew D Smith
#
# Author: Xiaojing Ji
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301 USA


library("optparse")

option.list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="input/output file prefix"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input directory"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output pdf file")
); 

opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser, positional_arguments=T);

prefix <- opt$options$prefix
din <- opt$options$input
output <- opt$options$output
#prefix <- "test"
#din <- "~/codes/epievo/test/output/mcmc"
#output <- "~/codes/epievo/test/output/mcmc/test.pdf"

files.forward <- list.files(path=din, pattern=paste0(prefix, ".forward$"))
nodes <- sapply(strsplit(files.forward, ".", fixed=T), function(x) x[2])
names(files.forward) <- nodes

files.mcmc <- list.files(path=din, pattern=paste0(prefix, ".mcmc$"))
names(files.mcmc) <- sapply(strsplit(files.mcmc, ".", fixed=T), function(x) x[2])

pdf(file=output, width=6.7, height=6.7)
par(mfrow=c(4, 4))
par(mar=c(2, 2, 1, 1))
par(mgp = c(0.5, 0, 0))
par(tcl = -0.17)
par(cex.lab=0.8, cex.axis=0.8, cex.main=1)

node <- 1
dat.forward <- read.table(file.path(din, files.forward[node]), header=T)
dat.mcmc <- read.table(file.path(din, files.mcmc[node]), header=T)

itrs <- unique(dat.mcmc$ITR)
batch <- nrow(dat.forward)
stats.name <- colnames(dat.forward)[3:ncol(dat.forward)]

for (stat in stats.name) {
  # plot batch mean over iterations
  mean.forward <- mean(dat.forward[, stat])
  mean.mcmc <- sapply(itrs,
                      function (x) mean(dat.mcmc[dat.mcmc$ITR==x, stat]))
  
  # plot(mean.mcmc, type="l", xlab="itr", ylab=stat, main=stat,
  #      ylim=range(mean.forward, mean.mcmc))
  # abline(h=mean.forward, lty=2, col="red")
  
  # plot final distributions
  samples.forward <- dat.forward[, stat]
  samples.mcmc <- dat.mcmc[(nrow(dat.mcmc)-batch+1):nrow(dat.mcmc), stat]
  if (substr(stat, 1, 1) == "J") {
    max_jumps <- max(samples.forward, samples.mcmc)
    hist(samples.forward, col=rgb(1,0,0,0.5), xlim=c(0, max_jumps),
         breaks=seq(0-0.5, max_jumps+0.5, by=1), ylim=c(0, batch),
         main=stat, xlab=stat)
    hist(samples.mcmc, col=rgb(0,0,1,0.5), add=T,
         breaks=seq(0-0.5, max_jumps+0.5, by=1))
    legend("topright", c("true", "sampled"),
           bty="n", fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
  } else {
    hist(samples.forward, col=rgb(1,0,0,0.5), main=paste0(stat, ".forward"), xlab=stat)
    hist(samples.mcmc, col=rgb(0,0,1,0.5), main=paste0(stat, ".mcmc"), xlab=stat)
  }
}

dev.off()
