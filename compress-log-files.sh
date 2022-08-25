#!/usr/bin/env bash
set -euo pipefail

DT=$(date +"%Y-%m-%d_%T")
DEST_DIR=archived_logs


echo "Adding current datatime as a tag: $DT"
mkdir -p $DEST_DIR

for dir in logs cluster_logs benchmarks .snakemake; do
  archive=${DEST_DIR}/${DT}_${dir}.tar.gz
  if [ -d "$dir" ]; then
    echo "Compressing $dir into $archive"
    tar --checkpoint=10000 -czf $archive $dir && rm -r $dir 
  else
    echo "$dir does not exist, skipping."
  fi
done 

