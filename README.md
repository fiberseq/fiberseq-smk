# `fiberseq-smk`: A Snakemake for calling **fiberseq**

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/StergachisLab/fiberseq-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/StergachisLab/fiberseq-smk/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for making **fiberseq** calls.



# Installation

## Install dependencies
### UCSC tools
Make sure you have UCSC tools installed and in your path.  You can get them from here: [https://hgdownload.soe.ucsc.edu/admin/exe/](https://hgdownload.soe.ucsc.edu/admin/exe/)

If you are on `hyak` you can add my copy to your path by adding this to your `.bashrc`:
```
PATH=$PATH:/mmfs1/home/mvollger/software/ucsc_bin
```
### `fibertools-rs`
See the [fibertools-rs install instructions](https://github.com/mrvollger/fibertools-rs#install) and make sure `ft` is in your path.

If you are on `hyak` you can add my copy to your path by adding this to your `.bashrc`:
```bash
PATH=$PATH:/mmfs1/gscratch/stergachislab/mvollger/projects/large_home/.cargo/bin/
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

### `cutnm`
Download and add cutnm to your path
```bash
wget https://raw.githubusercontent.com/sjneph/cutnm/master/src/cutnm
chmod 755 cutnm
export PATH=$PATH:/path/to/cutnm
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
    test=.test/subreads.bam `# path to the subreads, and the key sets the sample name` \
    ref=.test/ref.fa `# reference to align results to`  
```
If you find this too verbose you can instead include the config options in a configuration file:
```bash
snakemake --profile profile/local --configfile config.yaml 
```
And then specify the options in the `config.yaml`, for example:
```yaml
# sets the conda env for the jobs, always the same
env: fiberseq-smk
# choose any sample name followed by a path to the subreads
test: .test/subreads.bam
# optional path to a reference fasta to align the results to
ref: .test/hg38.analysisSet.fa
```

# Running with precomputed CCS data
Add to the config line the following option:
```bash
    input_type=ccs
```
And replace the subread bam file(s) with ccs file(s).

# Running on a cluster with distributed resources
Replace the `profile/local` argument with `profile/checkpoint` or `profile/compute` or make your own snakemake profile.

# Multiplexed data
If you have multiplexed data, you have to pass a `ccs` bam into the pipeline that has already been processed with `lima`. In general, you can pass a pre-generated ccs bam file to save on compute. For directions on this see [section above](#running-with-precomputed-ccs-data) on `ccs` inputs.

# Use `fibertools-rs` for m6A predictions
If you want to use `fibertools-rs` for m6A predictions, you can add the following to your config options:
```bash
    predict_with_hifi=True
```
You will also likely have a precomputed `ccs` bam file, in which case you will also add:
```bash
    input_type=ccs
```

# Test case
**Before running your own data** please run this small test case included in the repository. 
```bash
snakemake \
  --profile profile/local `# sets up where and how jobs are submitted` \
  --config \
    env="fiberseq-smk" `# sets the conda env for the jobs` \
    test=.test/subreads.bam `# path to the subreads, and the key sets the sample name` \
    ref=.test/ref.fa `# reference to align results to` \
  -p 
```
In addition, please clear this test case and then try again with the distributed resources (e.g. `profile/compute/`). 


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
results/test/test.unaligned.fiberseq.bam
results/test/test.unaligned.fiberseq.bam.pbi
# bed files
results/test/bed/test.unaligned.cpg.bed.gz
results/test/bed/test.unaligned.m6a.bed.gz
results/test/bed/test.unaligned.msp.bed.gz
results/test/bed/test.unaligned.nuc.bed.gz
```
</td>
<td>

```bash
# aligned bam files 
results/test/test.aligned.fiberseq.bam
results/test/test.aligned.fiberseq.bam.bai
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

## Added bam tags and definitions.
- MSP: methylation sensitive patch, defined as being any stretch of sequence between nucleosomes that has methylation.
- as: A bam tag with an array of MSP start sites
- al: A bam tag with an array of MSP lengths
- ns: A bam tag with an array of nucleosome start sites
- nl: A bam tag with an array of nucleosome lengths
- MM/ML: Bam tags for sorting m6a and 5mC methylation information. See the SAM spec for details.

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