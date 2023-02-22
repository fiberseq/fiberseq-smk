#!/bin/bash
# author : sjn
# date : Dec.2022

set -euox pipefail

# For a background bin size along the genome, calculate the proportion of di-nucleosome footprints
#   Numerator: nucleosome footprints overlapping that bin by at least 1 bp that are 220-450 bp in length
#   Denominator: nucleosome footprints overlapping that bin by at least 1 bp that are 80-450 bp in length
# Shift the bin by 50bp and repeat.

if [[ $# != 3 ]]; then
  printf "Require a reference's fai file, a zipped nucleosome BED file from fiberseq-smk, and an output file\n"
  exit -1
fi

ref_fai=$1
exp_nuc=$2 # fiber-all-table.tbl.gz
output=$3

smallest=80 # bp
dinuc_smallest=300 # bp
largest=1000 # bp
winsz=100 # looking +/- $winsz bp

tmpd=/tmp/`whoami`/$$
rm -rf ${tmpd}
mkdir -p ${tmpd}

cut -f1-2 ${ref_fai} \
  | awk 'BEGIN {OFS="\t"} ; { \
          if ( $2 > 0 ) { \
            print $1, 0, $2-1; \
          } \
        }' \
  | sort-bed - \
  > ${tmpd}/ref.sizes

# ignore last comma in ref_nuc_lengths/ref_nuc_starts lists
hck -z -F \#ct -F st -F ref_nuc_starts -F ref_nuc_lengths ${exp_nuc} \
  | awk 'NR > 1' \
  | awk '$3 != "."' \
  | awk -v m=${smallest} -v n=${largest} \
      'BEGIN {OFS="\t"} ; { \
        lng=split($4,a,","); \
        split($3,b,","); \
        for(i=1;i<lng;++i) { \
          if ( a[i] >= m && a[i] <= n ) { \
            print $1, b[i], b[i]+a[i]; \
          } \
        } \
      }' \
  | sort-bed --max-mem 25G - \
  > ${tmpd}/f

bedops --chop 1 --stagger 50 ${tmpd}/ref.sizes \
  | bedops -u --range -1:-1 - \
  | bedmap --ec --sweep-all --range ${winsz} --echo --echo-map-size - ${tmpd}/f \
  | awk -v m=${smallest} -v n=${largest} -v dm=${dinuc_smallest} \
      'BEGIN {FS="|"; OFS="\t"} ; { \
        num=0; denom=0; \
        lng=split($NF, a, ";"); \
        for(i=1;i<=lng;++i) { \
          if (a[i]>=m && a[i]<=n) { \
            denom+=1; \
            if (a[i]>=dm) { num+=1; } \
          } \
        } \
        if ( denom > 0 ) { print $1, num"/"denom, num/denom } \
        else { print $1, "0/0", 0 } \
      }' \
  > ${output}

rm -rf ${tmpd}

exit 0
