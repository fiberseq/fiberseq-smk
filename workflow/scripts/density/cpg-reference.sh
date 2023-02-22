#!/bin/bash
# author : sjn
# date : Dec.2022

set -euox pipefail

if [[ $# != 2 ]]; then
  printf "Require a reference fasta file and an output file\n"
  exit -1
fi

inp=$1
outp=$2

awk 'BEGIN {OFS="\t"; pos=0; prev=""; chr=""} ; { \
      if ( $1 ~ /^>/ ) { chr=substr($1, 2); prev=""; pos=0; next; } \
      lng=split(toupper($1), a, ""); fst=a[1]; \
      if ( fst == "G" && prev == "C" ) { \
        print chr, pos-1, pos; \
      } \
      prev=fst; \
      pos+=1; \
      for(i=2;i<=lng;++i) { \
        if ( a[i] == "G" && prev == "C" ) { \
          print chr, pos-1, pos; \
        } \
        pos+=1; \
        prev=a[i]; \
      } \
    }' ${inp} \
  | sort-bed - \
  > ${outp}

exit 0
