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
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="input/output file prefix"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input directory"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output pdf file"),
  make_option(c("-P", "--param"), type="character", default=NULL, 
              help="true parameter file")
); 

opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser, positional_arguments=T);

prefix <- opt$options$prefix
din <- opt$options$input
output <- opt$options$output
true.param.file <- opt$options$param
proposal.codes <- as.numeric(opt$args)+1

proposal.texts <- c("direct", "unif", "poisson", "forward")
proposals <- proposal.texts[proposal.codes]

stats <- list()
trace.mean <- list()
trace.var <- list()
params <- list()

# LOAD DATA
for (proposal in proposals) {
  stats[[proposal]] <- read.table(file.path(din, proposal,
                                          sprintf("%s.stats", prefix)),
                                header=T, sep="\t")
  trace.mean[[proposal]] <- read.table(file.path(din, proposal,
                                               sprintf("%s.trace.mean", prefix)),
                                     header=T, sep="\t")
  trace.var[[proposal]] <- read.table(file.path(din, proposal,
                                              sprintf("%s.trace.var", prefix)),
                                    header=T, sep="\t")
  
  if (!is.null(opt$options$param)) {
    params[[proposal]] <- read.table(file.path(din, proposal,
                                              sprintf("%s.stats.param", prefix)),
                                    header=T, sep="\t")
  }
}

labels <- colnames(stats[[proposals[1]]])
################################################################################

pdf(file=output, width=6.7, height=6.7)
par(mfrow=c(4, 4))
par(mar=c(2, 2, 1, 1))
par(mgp = c(0.5, 0, 0))
par(tcl = -0.17)
par(cex.lab=0.8, cex.axis=0.8, cex.main=1)

library(RColorBrewer)
cols.stats <- brewer.pal(8, "Set1")[1:length(proposals)]
cols.mean <- brewer.pal(8, "Dark2")[1:length(proposals)]
cols.var <- brewer.pal(8, "Paired")[4:(4+length(proposals))]
cols.param <- rev(brewer.pal(8, "Set1"))[1:length(proposals)]

ltys <- rep(1, length(proposals))

################################################################################
# distribution
################################################################################

for (i in 1:(length(labels) / 2)) {
  J <- list()
  D <- list()
  maxJ <- 0
  maxD <- 0
  for (k in 1:length(proposals)) {
    J[[k]] <- density(stats[[proposals[k]]][, i])
    D[[k]] <- density(stats[[proposals[k]]][, 8+i])
    maxJ <- max(maxJ, max(J[[k]]$y))
    maxD <- max(maxD, max(D[[k]]$y))
  }
  
  if (maxJ != 0) {
    plot(J[[1]], type="l",
         main="", font.main = 1, ylim=c(0, maxJ),
         lwd=1, col=cols.stats[1], lty=ltys[1],
         xlab=labels[i])
    for (k in 1:length(proposals)) {
      lines(J[[k]], lwd=1, col=cols.stats[k], lty=ltys[k])
    }
    legend("topright", proposals, col=cols.stats, lty=ltys, bty="n")
  }
  
  if (maxD != 0) {
    plot(D[[1]], type="l",
         main="", font.main = 1, ylim=c(0, maxD),
         lwd=1, col=cols.stats[1], lty=ltys[1],
         xlab=labels[8+i])
    for (k in 1:length(proposals)) {
      lines(D[[k]], lwd=1, col=cols.stats[k], lty=ltys[k])
    }
    legend("topright", proposals, col=cols.stats, lty=ltys, bty="n")
  }
}

################################################################################
# trace.mean
################################################################################

for (i in 1:(length(labels) / 2)) {
  means.J <- matrix(trace.mean[[proposals[1]]][, i],
                    ncol=length(proposals),
                    nrow=length(trace.mean[[proposals[1]]][, i]))
  means.D <- matrix(trace.mean[[proposals[1]]][, 8+i],
                    ncol=length(proposals),
                    nrow=length(trace.mean[[proposals[1]]][, 8+i]))
  if (k > 1)
    for (k in 2:length(proposals)) {
      means.J[, k] <- trace.mean[[proposals[k]]][, i]
      means.D[, k] <- trace.mean[[proposals[k]]][, 8+i]
    }

  if (diff(range(means.J)) != 0) {
    matplot(means.J, type="l",
            col=cols.mean, lty=ltys,
            xlim=c(0, nrow(means.J)), ylim=range(means.J),
            xlab="MCMC time", ylab=paste0("Mean ", labels[i]), main="")
    legend("topright", proposals,
           col=cols.mean, lty=ltys, bty="n")
  }

  if (diff(range(means.D)) != 0) {
    matplot(means.D, type="l",
            col=cols.mean, lty=ltys,
            xlim=c(0, nrow(means.D)), ylim=range(means.D),
            xlab="MCMC time", ylab=paste0("Mean ", labels[8+i]), main="")
    legend("topright", proposals,
           col=cols.mean, lty=ltys, bty="n")
  }
}

################################################################################
# trace.var
################################################################################
# for (i in 1:(length(labels) / 2)) {
#   vars.J <- matrix(trace.mean[[proposals[1]]][, i],
#                    ncol=length(proposals),
#                    nrow=length(trace.mean[[proposals[1]]][, i]))
#   vars.D <- matrix(trace.mean[[proposals[1]]][, 8+i],
#                    ncol=length(proposals),
#                    nrow=length(trace.mean[[proposals[1]]][, 8+i]))
#   
#   if (k > 1)
#     for (k in 2:length(proposals)) {
#       vars.J[, k] <- trace.var[[proposals[k]]][, i]
#       vars.D[, k] <- trace.var[[proposals[k]]][, 8+i]
#     }
#   
#   if (diff(range(vars.J)) != 0) {
#     matplot(vars.J, type="l",
#             col=cols.var, lty=ltys,
#             xlim=c(0, nrow(vars.J)), ylim=range(vars.J),
#             xlab="MCMC time", ylab=paste0("Var ", labels[i]), main="")
#     legend("topright", proposals,
#            col=cols.var, lty=ltys, bty="n")
#   }
#   
#   if (diff(range(vars.D)) != 0) {
#     matplot(vars.D, type="l",
#             col=cols.var, lty=ltys,
#             xlim=c(0, nrow(vars.D)), ylim=range(vars.D),
#             xlab="MCMC time", ylab=paste0("Var ", labels[8+i]), main="")
#     legend("topright", proposals,
#            col=cols.var, lty=ltys, bty="n")
#   }
# }

################################################################################
# parameter
################################################################################
par(mfrow=c(3, 2))

load.param <- function(file) {
  table <- read.table(file, header=F)
  return(list(st0=table[1, 2], st1=table[1, 3],
              bl0=table[2, 2], bl1=table[2, 3],
              rt0=table[3, 2], rt1=table[3, 3]))
}

if (!is.null(opt$options$param)) {
  params.true <- load.param(true.param.file)
  parameters <- names(params.true)

  titles <- c("stationary T00", "stationary T11",
              "log-baseline l00", "log-baseline l11",
              "root T00", "root T11")
  names(titles) <- parameters

  for (par in parameters) {
    param.est <- matrix(params[[proposals[1]]][, par],
                        ncol=length(proposals),
                        nrow=length(params[[proposals[1]]][, par]))
    if (k > 1)
      for (k in 2:length(proposals)) {
        param.est[, k] <- params[[proposals[k]]][, par]
      }
    matplot(param.est, type="l",
            col=cols.param, lty=ltys,
            xlim=c(0, nrow(param.est)),
            ylim=range(c(range(param.est), params.true[[par]])),
            xlab="MCMC time", ylab=par, main=titles[par])
    abline(h=params.true[[par]], lwd=2, col="red", lty=2)
    legend("topright", proposals,
           col=cols.param, lty=ltys, bty="n")
  }
}

dev.off()