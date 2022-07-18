# `fiberseq-smk`: A Snakemake for calling **fiberseq**

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/StergachisLab/fiberseq-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/StergachisLab/fiberseq-smk/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for making **fiberseq** calls.



# Installation

```bash
conda create -n fiberseq-smk
mamba env update -n fiberseq-smk --file workflow/envs/env.yml 
# Add the SMRTlink tools to your path
export PATH=$PATH:/path/to/smrtlink/tools
```

# Usage

Send your jobs to the cluster with:
```bash
snakemake \
    --profile profile/ 
    --config env="fiberseq-smk" \
    --configfile .test/config.yml \
    --set-scatter chunks=400 \
    -p
```
# TODO
- [x] Make an extract tool
- [ ] Add a sample tag to the bam header.
- [ ] Add a pipeline version to the bam header (git commit).
- [x] Add primrose
- [x] Add unknown case to end of fiber calls

# Workflow

![alt text](./images/dag.png)