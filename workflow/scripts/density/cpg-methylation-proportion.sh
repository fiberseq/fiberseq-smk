#!/bin/bash
# author : sjn
# date : Dec.2022

set -euox pipefail

if [[ $# != 3 ]]; then
  printf "Require a reference's cpg locations, fiberseq-smk's all table output, and an output-dir\n"
  printf "Some references are available at /mmfs1/gscratch/stergachislab/sjn/data/cpg/\n"
  exit -1
fi

ref_cpg=$1
exp_cpg=$2
output=$3

tmpd=/tmp/`whoami`/$$
if [[ -d $tmpd ]]; then
  rm -rf $tmpd
fi
mkdir -p $tmpd

# use field 4 too, because cpg field gets too big for typical BEDOPS ID field
hck -z -F \#ct -F st -F en -F fiber -F ref_5mC $exp_cpg \
  | awk 'NR > 1 && $2 != $3' \
  | sort-bed - \
  | tee $tmpd/f \
  | bedmap --ec --sweep-all --count $ref_cpg - \
  > $tmpd/g

# ignore last comma
awk 'BEGIN {OFS="\t"} ; { \
      lng=split($NF,a,","); \
      for(i=2;i<lng;++i) { \
        if ( a[i] >= 0 ) { \
          print $1, a[i], a[i]+1; \
        } \
      } \
    }' $tmpd/f \
  | sort-bed --max-mem 25G - \
  | bedmap --faster --sweep-all --exact --delim "\t" --echo --count $ref_cpg - \
  | paste - $tmpd/g \
  | awk 'BEGIN {OFS="\t"} ; { \
      if ( $5 > 0 ) { print $1, $2, $3, $4"/"$5, $4/$5; } \
      else { print $1, $2, $3, 0"/"0, 0; } \
    }' \
  > $output

rm -rf $tmpd

exit 0
