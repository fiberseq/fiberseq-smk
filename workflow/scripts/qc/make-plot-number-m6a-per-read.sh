#!/bin/bash
# author : sjn
# date : Aug 22, 2022

# histogram of per read #m6A on read / # of A/T bases on read
# histogram of the ec field's distribution

# this depends upon rust's cargo and fibertools developed by Mitchell
#   module load gcc/10.2.0
#   module load cmake/3.21.1
#   cargo install --git https://github.com/mrvollger/fibertools-rs
#   ft extract --m6a - <BAM> | cut -f 4,10 | head
# The tools takes into account M6A vs CpG methylations

# Can you change the name of "qc_m6A.ec.pdf" to "qc_ccs_passes.pdf" - now factored out to Snakefile

set -exuo pipefail

if [[ $# != 5 ]]; then
    printf "Expect $0 <sample-name> <input-file> <output-pdf> <output-ec-pdf> <output-stat.txt>\n"
    exit 1
fi

samplenm=$1
inp=$2 # "*_unaligned.fiberseq.bam"
outpdf=$3
ec_outpdf=$4
outstat=$5

ftype=m6a
#ec_outpdf=$(basename "${outpdf}" | sed 's/\(.*\)\..*/\1/')".ec.pdf"
#ec_outpdf=$(dirname "${outpdf}")"/${ec_outpdf}"
#ec_outpdf=$(dirname "${outpdf}")"/qc_ccs_passes.pdf"
tmpd=/tmp/$(whoami)/$$
if [[ ! -s $inp ]]; then
    printf "Problem finding 1 file: %s\n" $inp
    exit 1
fi

if [ -s $tmpd ]; then
    rm -rf $tmpd
fi
mkdir -p $tmpd
mkdir -p $(dirname "${outpdf}")

hck -F total_m6a_bp -F total_AT_bp -F ec -z $inp >$tmpd/$samplenm.$ftype

R --no-save --quiet <<__R__
  maxx <- 0.30
  max_coverage <- 30
  s <- read.table("$tmpd/$samplenm.$ftype", header=TRUE, sep="\t", row.names=NULL)
  p <- s[,1]/s[,2]
  coverage <- s[,3]
  pcov = 100*length(coverage[coverage>max_coverage])/length(coverage)
  coverage[coverage > max_coverage] <- max_coverage

  mycol <- "darkgreen"
  pdf("$outpdf")
  df <- as.data.frame(cbind(coverage, p))
  colnames(df) <- c("cvg", "prop")
  ndf <- subset(df, cvg>0)
  p <- ndf[,2]
  mp <- median(p) # median calculated after removing those with low coverage
  u <- subset(ndf, prop>maxx)
  v <- subset(ndf, prop<=maxx)
  if ( dim(u)[1] > 0 ) { u[,2] <- maxx }
  ndf <- as.data.frame(rbind(v, u))
  xlab <- "Proportion of adenines with m6A per read (for reads with 0+ CCS coverage)"
  h <- hist(ndf[,2], axes=FALSE, main="$samplenm", xlab=xlab, ylab="Count", breaks=1000, xlim=c(0,maxx))
  mycol <- "darkgreen"
  abline(v=mp, col=mycol, lty=1)
  axis(1)
  axis(2)

  mxhist <- max(h[["counts"]])
  msg1 <- paste(round(mp, 2))
  rtoff <- 0.025
  text(mp+rtoff, mxhist, msg1, col=mycol)

  b <- 100*dim(u)[1]/dim(ndf)[1]
  msg2 <- paste(round(b,1), "%>", maxx, sep="")
  text(maxx-rtoff, mxhist/2, msg2, col=mycol)
  dev.off()

  rtoff <- 1.5
  pdf("$ec_outpdf")
  h <-  hist(coverage, xlim=c(0,max_coverage), main="$samplenm", xlab="Subread coverage of each CCS read", ylab="Count", breaks=1000)
  mc <- median(coverage)
  abline(v=mc, col=mycol, lty=1)
  text(mc+rtoff, max(h[["counts"]]), round(mc, 3), col=mycol)
  text(max(coverage)-rtoff*3, max(h[["counts"]])/2, paste(round(pcov,2), "%>", max_coverage, sep=""), col=mycol)
  dev.off()

  # use p and mp which have been filtered for coverage>0
  pv1 <- round(100 * length(p[p<0.01])/length(p), 1)
  pv2 <- round(100 * length(p[p<0.2])/length(p), 1)
  stats_file <- "$outstat"
  cat("# Note: ***#m6A/#ATs stats***\n", file=stats_file, sep="", append=FALSE)
  cat("# Note: #m6A/#ATs filtered to reads with CCS Coverage > 0\n", file=stats_file, sep="", append=TRUE)
  cat("Proportion(#m6A/#ATs)<0.01=", pv1, "%\n", file=stats_file, sep="", append=TRUE)
  cat("Proportion(#m6A/#ATs)<0.2=", pv2, "%\n", file=stats_file, sep="", append=TRUE)
  cat("Median(#m6A/#ATs)=", mp, "\n", file=stats_file, sep="", append=TRUE)
  cat("Percent(#m6A/#ATs)>", maxx, "=", paste(round(b,1)), "%\n", file=stats_file, sep="", append=TRUE)
  cat("Median(EQ)", "=", round(mc,3), "\n", file=stats_file, sep="", append=TRUE)
__R__

rm -rf $tmpd

exit 0
