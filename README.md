# `fiberseq-smk`: A Snakemake workflow for making **fiberseq** data

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/StergachisLab/fiberseq-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/StergachisLab/fiberseq-smk/actions?query=branch%3Amain+workflow%3ATests)

# Alpha warning!
This code is in an alpha state and will be changing without notice or versioning.

# Installation
See [docs/INSTALL.md](docs/INSTALL.md) for installation instructions.

# Input data
The input data is a `ccs` bam file(s) **with HiFi kinetics**. You can generate this file from PacBio subreads using `pbccs`. 

### Multiplexed ccs data
If the data is multiplexed you must first process it with `lima` and pass in demultiplexed bam file(s).

### Subread input (deprecated)
You can pass subread bam(s) to `fiberseq-smk` by adding `input_type=subreads` to your config options. This will run `pbccs` on the subreads and then run the rest of the pipeline. This feature is deprecated and will be removed in the future.
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

# Running on a cluster with distributed resources
You can configure `fiberseq-smk` to run on a distributed resource manager (e.g. SLURM, PBS, SGE, etc.) by creating a [snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). An example configuration for a SLURM cluster is included in the [profile/compute](profile/compute) directory. To use this profile you can add the following to your command:
```bash
--profile profile/compute
```
Several additional profiles are made available under [profile/](profile/) which can be used as a starting point for your own cluster configuration.
# Test case
**Before running your own data** please run this small test case included in the repository shown in the `Usage` section. In addition, please clear this test case and then try again with the distributed resources (e.g. `profile/compute/`). 

# Examples 
You can find examples of running the pipeline in [docs/EXAMPLES.md](docs/EXAMPLES.md).

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
See [docs/bed.md](docs/bed.md) for instructions on how to make bed files from your fiberseq bam(s).

# TODO
- [ ] Add a pipeline version to the bam header (git commit).
- [ ] Add env version to the output somewhere. 
