# `fiberseq-smk`: A Snakemake for calling **fiberseq**

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/StergachisLab/fiberseq-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/StergachisLab/fiberseq-smk/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for making **fiberseq** calls.



# Installation

```bash
git clone https://github.com/StergachisLab/fiberseq-smk
cd fiberseq-smk
conda create -n fiberseq-smk
mamba env update -n fiberseq-smk --file workflow/envs/env.yml 
conda activate fiberseq-smk
# Add the SMRTlink tools to your path
export PATH=$PATH:/path/to/smrtlink/tools
```

# Test case
**Before running your own data** please run this small test case included in the repository.

```bash
## Local test case
```bash
snakemake \
  --profile profile/local \
  --config \
    env="fiberseq-smk" \
    test=.test/subreads.bam \
    ccs=.test/ccs.bam \
    ref=.test/ref.fa \
  -p 
```

## Cluster test case
You can also try submitting this test case to the cluster. Note the profile is currently configured for the Stergachis Lab cluster, you may need to modify it.
```bash
snakemake \
  --profile profile/compute \
  --config \
    env="fiberseq-smk" \
    test=.test/subreads.bam \
    ccs=.test/ccs.bam \
    ref=.test/ref.fa \
  -p 
```

# Usage

Send your jobs to the cluster with:
```bash
snakemake \
    --profile profile/checkpoint \
    --config \
        env="fiberseq-smk" \
        ref=.test/ref.fa \
        subreads=.test/subreads.bam \
    -p
```
## Multiplexed data
If you have multiplexed data, you have to pass a `ccs` bam into the pipeline that has already been processed with `lima`. e.g.
```bash
snakemake \
    --profile profile/checkpoint \
    --config \
        env="fiberseq-smk" \
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
# path to the reference genome, optional input
ref: .test/hg38.analysisSet.fa
# path to the ccs bam file, optional input
ccs: .test/ccs.bam
```

# Output files
Example output files if your sample name is `test`:
<table>
<tr>
<th> Unaligned outputs </th>
<th> Aligned outputs </th>
</tr>
<tr>
<td>

```bash 
# bam files 
results/test/unaligned.fiberseq.bam
results/test/unaligned.fiberseq.bam.pbi
# bed files
results/test/unaligned.cpg.bed.gz
results/test/unaligned.m6a.bed.gz
results/test/unaligned.msp.bed.gz
results/test/unaligned.nuc.bed.gz
```
</td>
<td>

```bash
# aligned bam files 
results/test/aligned.fiberseq.bam
results/test/aligned.fiberseq.bam.bai
# aligned bed files
results/test/aligned.cpg.bed.gz
results/test/aligned.m6a.bed.gz
results/test/aligned.msp.bed.gz
results/test/aligned.nuc.bed.gz
# aligned big bed file
results/test/aligned.nuc.bed.bb
results/test/aligned.m6a.bed.bb
results/test/aligned.cpg.bed.bb
results/test/aligned.msp.bed.bb
```
</td>
</tr>
</table>



# TODO
- [x] Make an extract tool
- [x] Add primrose
- [x] Add unknown case to end of fiber calls
- [x] Add a ccs input option
- [ ] Add a sample tag to the bam header.
- [ ] Add a pipeline version to the bam header (git commit).
- [ ] Add env version to the output somewhere. 
- [ ] Can we check for de-multiplexing automatically?


# Workflow

![alt text](./images/dag.png)