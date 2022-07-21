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
    --profile profile/checkpoint \
    --config env="fiberseq-smk" \
    --config \
        ref=.test/ref.fa \
        subreads=.test/subreads.bam \
    -p
```
## Multiplexed data
If you have multiplexed data, you have to pass a `ccs` bam into the pipeline that has already been processed with `lima`. e.g.
```bash
snakemake \
    --profile profile/checkpoint \
    --config env="fiberseq-smk" \
    --config \
        ref=.test/ref.fa \
        subreads=.test/subreads.bam \
        ccs=.test/ccs.bam \
    -p
```
In general, you can pass a pre-generated ccs bam file to save on compute.


## Example config file
Instead of passing configuration options on the command line, you can pass a config file. Here is an example, you often call this file `config.yml`:
```yaml
# sample name and then path to the subreads
test: .test/subreads.bam
# path to the reference genome
ref: /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.fa
```

# TODO
- [x] Make an extract tool
- [ ] Add a sample tag to the bam header.
- [ ] Add a pipeline version to the bam header (git commit).
- [x] Add primrose
- [x] Add unknown case to end of fiber calls
- [x] Add a ccs input option

# Workflow

![alt text](./images/dag.png)