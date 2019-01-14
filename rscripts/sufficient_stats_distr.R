#!/usr/local/bin/Rscript
# Copyright (C) 2016-2018 University of Southern California
#                         Andrew D Smith
# Author: Andrew D. Smith, Xiaojing Ji
# 
#  This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundat.corrion; either version 2 of the License, or
# (at your option) any later version.
# 
#   This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
# along with this software; if not, write to the Free Software
# Foundat.corrion, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301 USA


library("optparse")

option.list = list(
  make_option(c("-t", "--true"), type="character", default=NULL, 
              help="true stats file prefix"),
  make_option(c("-e", "--est"), type="character", default=NULL, 
              help="estimated stats file prefix"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output pdf file")
); 

opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser, positional_arguments=T);

filepath.true <- opt$options$true
filepath.est <- opt$options$est
output <- opt$options$output
proposal.codes <- as.numeric(opt$args)+1

proposal.texts <- c("direct", "unif", "poisson", "forward")
proposals <- proposal.texts[proposal.codes]

# load the summary stats data
input <- sprintf("%s.stats", filepath.true)
input2 <- sprintf("%s.stats", filepath.est)

ftrace.mean <- sprintf("%s.trace.mean", filepath.est)
ftrace.var <- sprintf("%s.trace.mean", filepath.est)

tab <- read.table(input, header=T, sep="\t")
tab2 <- read.table(input2, header=T, sep="\t")
trace.mean <- read.table(ftrace.mean, header=T, sep="\t")
trace.var <- read.table(ftrace.var, header=T, sep="\t")

labels <- colnames(tab)
histlim <- nrow(tab)

################################################################################

pdf(file=output, width=6.7, height=6.7)
par(mfrow=c(4, 4))
par(mar=c(2, 2, 1, 1))
par(mgp = c(0.5, 0, 0))
par(tcl = -0.17)
par(cex.lab=0.8, cex.axis=0.8, cex.main=1)

# total distribution
total <- sapply(1:nrow(tab), function (x) sum(tab[x, 1:8]))
total2 <- sapply(1:nrow(tab2), function (x) sum(tab2[x, 1:8]))

max_total <- max(max(total), max(total2))

hist(total, col=rgb(1,0,0,0.5), xlim=c(0, max_total),
     breaks=seq(0-0.5, max_total+0.5, by=1),
     main="total jumps", font.main = 1,
     xlab="# jumps", ylim=c(0, histlim))
hist(total2, col=rgb(0,0,1,0.5), freq=T, add=T, xlim=c(0,max_total),
     breaks=seq(0-0.5, max_total+0.5, by=1))
legend("topright", c("true", "sampled"),
       bty="n", fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))


# J/D distribution
for (i in 1:(length(labels) / 2)) {
  if (sum(tab[, i]) > 0) {
    max_jumps <- max(max(tab[, i], tab2[, i]))
    
    hist(tab[, i], col=rgb(1,0,0,0.5), xlim=c(0, max_jumps),
         breaks=seq(0-0.5, max_jumps+0.5, by=1),
         ylim=c(0, histlim),
         #main=sprintf("%s, p.chi=%11.3e", labels[i], chi$p.value), font.main = 1,
         main="", font.main = 1,
         xlab=labels[i])
    hist(tab2[, i], col=rgb(0,0,1,0.5), add=T,
         breaks=seq(0-0.5, max_jumps+0.5, by=1)
         )
    legend("topright", c("true", "sampled"),
           bty="n", fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
  }
  if (sum(tab[, 8+i]) > 0) {
    d1 <- density(tab[, 8+i])
    d2 <- density(tab2[, 8+i])
    maxD <- max(max(d1$y), max(d2$y))
    
    plot(d1, col=rgb(1,0,0,0.5),
        main="", font.main = 1,
         ylim=c(0, maxD), lwd=2,
         xlab=labels[8+i])
    lines(d2, col=rgb(0,0,1,0.5), lwd=2)
    legend("topright", c("true", "sampled"),
           bty="n", fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
  }
}

dev.off()