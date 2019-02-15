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
  make_option(c("-t", "--true"), type="character", default=NULL, 
              help="true parameter file"),
  make_option(c("-e", "--est"), type="character", default=NULL, 
              help="estimated parameter file prefix"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output pdf file")
); 

opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser, positional_arguments=T);

filepath.true <- opt$options$true
filepath.est <- opt$options$est
output <- opt$options$output

################################################################################

# load true parameters
load.param <- function(file) {
  table <- read.table(file, header=F)
  return(list(st0=table[1, 2], st1=table[1, 3],
              bl0=table[2, 2], bl1=table[2, 3],
              rt0=table[3, 2], rt1=table[3, 3]))
}
params.true <- load.param(filepath.true)
parameters <- names(params.true)

pdf(file=output, width=6.7, height=6.7/2)
par(mar=c(3, 3, 2, 2))
par(mgp = c(1, 0, 0))
par(tcl = -0.25)
par(cex.lab=1, cex.axis=1, cex.main=1)
par(mfrow=c(1, 6))

for (par in parameters) {
  data <- read.table(paste0(filepath.est, "_", par, ".stats"), header=F)
  avg <- mean(data[, 1])
  sdev <- sd(data[, 1])
  yl <- c(min(0, avg-2*sdev), max(0, avg+2*sdev))
  barCenters <- barplot(avg, 
                        ylim=yl,
                        xlab=par, ylab="", main="")
  segments(barCenters, avg-sdev, barCenters, avg+sdev, lwd=2)

  arrows(barCenters, avg-sdev, barCenters, avg+sdev, length=0.05, angle=90, code=3)
  abline(h=params.true[par], col="red", lty=2)
}

dev.off()