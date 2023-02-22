#!/bin/bash
# author : sjn
# date : Dec.2022

set -euox pipefail

# For a 1kb bin along the genome, calculate the proportion of MSPs>150bp
#   Numerator: Number of MSPs overlapping that bin by at least 1 bp that are >=150bp in length
#   Denominator: Number of MSPs overlapping that bin by at least 1 bp
# Shift the bin by 50bp and repeat.

if [[ $# != 3 ]]; then
  printf "Require a reference's .fai file, a zipped 'all table output' file from fiberseq-smk, and an output file\n"
  exit -1
fi

ref_fai=$1 # first 2 columns show chr/size
exp_msp=$2
output=$3
winsz=100 # +/- bp
thold=150 # bp

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

# skip last ',' in ref_msp_*
hck -z -F \#ct -F st -F ref_msp_starts -F ref_msp_lengths ${exp_msp} \
  | awk 'NR > 1' \
  | awk '$3 != "."' \
  | awk \
      'BEGIN {OFS="\t"} ; { \
        lng=split($4,a,","); \
        split($3,b,","); \
        for(i=1;i<lng;++i) { \
          if ( a[i] > 0 ) { \
            print $1, b[i], b[i]+a[i]; \
          } \
        } \
      }' \
  | sort-bed --max-mem 25G - \
  > ${tmpd}/f

bedops --chop 1 --stagger 50 ${tmpd}/ref.sizes \
  | bedops -u --range -1:-1 - \
  | bedmap --ec --sweep-all --range ${winsz} --echo --echo-map-size - ${tmpd}/f \
  | awk -v t=${thold} \
      'BEGIN {FS="|"; OFS="\t"} ; { \
        num=0; \
        denom=split($NF, a, ";"); \
        for(i=1;i<=denom;++i) { \
          if (a[i]>=t) { \
            num+=1; \
          } \
        } \
        if ( denom > 0 ) { print $1, num"/"denom, num/denom } \
        else { print $1, "0/0", 0 } \
      }' \
  > ${output}

rm -rf ${tmpd}

exit 0
