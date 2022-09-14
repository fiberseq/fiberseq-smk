#!/bin/bash
# author : sjn
# date : Aug 25, 2022

# Histogram of per read length/#nucleosomes (bin sz=1, range 0-2000bp)
#  Collapse all reads > 2000 to 2000
#  Vertical line at median value between 0-300 bp/nuc 
# Can you actually make this with xlim 0-1500, and say on the chart the % of fibers in the bin >1500
# Get % reads with readlength/#nucs < 300 and put in stats file

set -euo pipefail

if [[ $# != 4 ]]; then
  printf "Expect $0 <sample-name> <input-file> <output-pdf> <output-stat.txt>\n"
  exit -1
fi

samplenm=$1
inp=$2 # "*_unaligned.fiberseq.bam"
outpdf=$3
outstat=$4

ftype=readlength.div.number.nucleosomes
tmpd=/tmp/`whoami`/$$

if [[ ! -s $inp ]]; then
  printf "Problem finding 1 file: %s\n" $inp
  exit -1
fi

if [[ -s $tmpd ]]; then
  rm -rf $tmpd
fi
mkdir -p $tmpd
mkdir -p $(dirname "${outpdf}")

ft extract --all - $inp \
  | cutnm fiber_length,nuc_starts \
  | awk '{ lng=gsub(/,/, "", $2); if ( lng > 1 ) { print int($1/lng) } }' \
  > $tmpd/$samplenm.$ftype

R --no-save --quiet << __R__
  ss <- scan("$tmpd/$samplenm.$ftype")
  mxh <- 1500
  s <- ss
  mxv <- 300
  nlarge <- length(s[s>mxh])
  s[s>mxh] <- mxh
  m <- median(s[s<=mxv])
  plarge <- nlarge/length(s)

  stats_file <- "$outstat"
  prop = round(100*length(s[s<=mxv])/length(s), 1)
  cat("# Note: ***readlength divided by #nucs stats***\n", file=stats_file, append=FALSE)
  cat("Proportion(readlength/#nucs<=", mxv, ")=", prop, "%\n", file=stats_file, sep="", append=TRUE)
  cat("Median(readlength/#nucs<=", mxv, ")=", m, "\n", file=stats_file, sep="", append=TRUE)

  mycol <- "darkgreen"
  pdf("$outpdf")
  h <- hist(s, xlim=c(0, mxh), breaks=mxh, axes=F, main=paste("$samplenm"), xlab="Read length (bp)/# nucleosome footprints per read", ylab="Count")
  abline(v=m, col=mycol, lty=1)
  msg <- paste(m)
  msg2 <- paste(round(plarge*100,1), "% >", mxh, "bp", sep="")

  xv <- 75
  offxv <- 8
  rtoff <- 300
  div <- 3
  lines(c(m+rtoff, m), c(max(h[["counts"]])/(div*1.05), 0), col=mycol)
  text(m+rtoff+3*offxv, max(h[["counts"]])/div, msg, col=mycol)
  text(mxh-100, max(h[["counts"]])/div, msg2, col=mycol)

  axis(1, seq(0, mxh, 100))
  axis(2)

  dev.off()
__R__

rm -rf $tmpd

exit 0
