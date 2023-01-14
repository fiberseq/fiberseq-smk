# `fiberseq-smk`: A Snakemake workflow for making **fiberseq** data

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/StergachisLab/fiberseq-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/StergachisLab/fiberseq-smk/actions?query=branch%3Amain+workflow%3ATests)

# Alpha warning!
This code is in an alpha state and will be changing without notice or versioning.

# Installation

## Install dependencies
### `fibertools-rs`
See the [fibertools-rs install instructions](https://github.com/mrvollger/fibertools-rs#install) and make sure `ft` is in your path.

If you are on `hyak` you can add my copy to your path by adding this to your `.bashrc`:
```bash
PATH=$PATH:/mmfs1/gscratch/stergachislab/mvollger/projects/large_home/.cargo/bin/
```
You will also need to add my PyTorch lib to your environment:
```bash
  export LIBTORCH=/mmfs1/gscratch/stergachislab/mvollger/projects/large_home/libs/libtorch-static-with-deps-1.13.0_cu116
  export LIBTORCH_CXX11_ABI=0
  export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
  export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
```

### `SMRTLINK` tools
```
# Add the SMRTlink tools to your path
export PATH=$PATH:/path/to/smrtlink/tools
```
If you are on `hyak` you can add our copy to your path by adding this to your `.bashrc`:
```bash
PATH=$PATH:/gscratch/stergachislab/install_dir/smrtlink/smrtcmds/bin/
```

## Install the workflow
```bash
git clone https://github.com/StergachisLab/fiberseq-smk
cd fiberseq-smk
conda create -n fiberseq-smk
mamba env update -n fiberseq-smk --file workflow/envs/env.yml 
conda activate fiberseq-smk
```


# Usage
You can run data using the following command, read the comments to learn more about the config options:
```bash
snakemake \
  --profile profile/local `# sets up where and how jobs are submitted` \
  --config \
    env="fiberseq-smk" `# sets the conda env for the jobs, always the same` \
    test=.test/ccs.bam `# path to the ccs reads with HiFi kinetics, and the key sets the sample name` \
    ref=.test/ref.fa `# reference to align results to`  
```
If you find this too verbose you can instead include the config options in a configuration file:
```bash
snakemake --profile profile/local --configfile config.yaml 
```
See `config/config.yml` for an example configuration file.

# Starting with SUBREADS instead of CCS
Add to the config line the following option:
```bash
    input_type=subreads
```
And replace the ccs bam file(s) with subread bam file(s).

# Running on a cluster with distributed resources
Replace the `profile/local` argument with `profile/checkpoint` or `profile/compute` or make your own snakemake profile.

# Multiplexed data
If you have multiplexed data, you have to pass a `ccs` bam into the pipeline that has already been processed with `lima`. In general, you can pass a pre-generated ccs bam file to save on compute. For directions on this see [section above](#running-with-precomputed-ccs-data) on `ccs` inputs.

# Use deprecated m6A calling with `ipdSummary`
If you want to use `ipdSummary` for m6A predictions, you can add the following to your config options:
```bash
    ipdsummary=True
```

# Test case
**Before running your own data** please run this small test case included in the repository shown in the `Usage` section. In addition, please clear this test case and then try again with the distributed resources (e.g. `profile/compute/`). 


# Output files
Example output files if your sample name is `test`:
```bash
# aligned bam files 
results/test/test.fiberseq.bam
# fiber table 
results/test/test.fiberseq.all.tbl.gz 
# Quality control directory with figures and html reports
results/test/qc/
```

## Added bam tags and definitions.
- MSP: methylation sensitive patch, defined as being any stretch of sequence between nucleosomes that has methylation.
- as: A bam tag with an array of MSP start sites
- al: A bam tag with an array of MSP lengths
- ns: A bam tag with an array of nucleosome start sites
- nl: A bam tag with an array of nucleosome lengths
- MM/ML: Bam tags for sorting m6a and 5mC methylation information. See the SAM spec for details.

# Making bed files to supplement your fiberdata
This optional step allows you to make some common bed files for fiberseq data.  These bed files are not required for the pipeline to run, but they are useful for visualizing your data. 
## Dependencies
For this pipeline the only additional dependency is UCSC tools, so make sure you have UCSC tools installed and in your path.  You can get them from here: [https://hgdownload.soe.ucsc.edu/admin/exe/](https://hgdownload.soe.ucsc.edu/admin/exe/)
If you are on `hyak` you can add my copy to your path by adding this to your `.bashrc`:
```
PATH=$PATH:/mmfs1/home/mvollger/software/ucsc_bin
```
## Generating bed files
To run the bed file generation pipeline grab your fiberseq bam file and the reference it was aligned to. Then you can run the pipeline with a command based on this example:
```bash
snakemake \
  --profile profile/local `# sets up where and how jobs are submitted` \
  --config \
    make_beds=True `# Tells the pipeline to make bed files instead of a fiberseq bam` \
    env="fiberseq-smk" `# sets the conda env for the jobs, always the same` \
    GM12878=GM12878.fiberseq.bam `# path to the fiberseq bam, and the key sets the sample name` \
    ref=hg38.fa `# reference that the fiberseq.bam`  
```
Separating these bed steps into their own pipeline allows you to merge multiple SMRT cells of fiberseq data from the same sample before making bed files. Which is the more common use case.

## Outputs of the bed pipeline
Example output files if your sample name is `test`:
<table>
<tr>
<th> Unaligned outputs </th>
<th> Aligned outputs </th>
</tr>
<tr>
<td>

```bash 
# bed files
results/test/bed/test.unaligned.cpg.bed.gz
results/test/bed/test.unaligned.m6a.bed.gz
results/test/bed/test.unaligned.msp.bed.gz
results/test/bed/test.unaligned.nuc.bed.gz
```
</td>
<td>

```bash
# aligned bed files
results/test/bed/test.aligned.cpg.bed.gz
results/test/bed/test.aligned.m6a.bed.gz
results/test/bed/test.aligned.msp.bed.gz
results/test/bed/test.aligned.nuc.bed.gz
# aligned big bed file
results/test/bigbed/test.aligned.nuc.bed.bb
results/test/bigbed/test.aligned.m6a.bed.bb
results/test/bigbed/test.aligned.cpg.bed.bb
results/test/bigbed/test.aligned.msp.bed.bb
```
</td>
</tr>
</table>

All the bed output files are [bed12 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). If the bed output beings with `aligned` all the bed file coordinates are in reference space, conversely if it begins with `unaligned` then all the coordinates are with respect to the unaligned fiberseq read. For both `aligned` and `unaligned` there are four types of bed files:

- `cpg`: Uses the block starts in bed12 format to describe the positions of `5mC` calls from `primrose`
- `m6a`: Uses the block starts in bed12 format to describe the positions of `m6a` calls from the pipeline
- `msp`: Uses the block starts and block lengths in bed12 format to describe the start and end of Methylation Sensitive Patches (MSPs)
- `nuc`: Uses the block starts and block lengths in bed12 format to describe the start and end of nucleosomes

For all records in all the bed files the first and last block position do not reflect real data, they are only there because bed12 format requires the first and last block to match the first and last position of the entire read.

# TODO
- [ ] Add a sample tag to the bam header.
- [ ] Add a pipeline version to the bam header (git commit).
- [ ] Add env version to the output somewhere. 
- [ ] Split out the bed part of the pipeline.

# Workflow

![alt text](./images/dag.png)
