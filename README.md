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
```bash
snakemake --cores 10 -p --set-scatter chunks=2  --show-failed-logs --configfile .test/config.yml -n
```

# Cluster usage
### Example:
```bash
mkdir -p logs_cluster
partition="compute"
partition="ckpt"

snakemake \
  --cores 8 \
  --jobs 500 \
  -p -r \
  --use-conda \
  --config env="fiberseq-smk" \
  --configfile .test/config.yml \
  --set-scatter chunks=8 \
  --cluster " sbatch --account=stergachislab --partition=$partition --nodes=1 --ntasks-per-node={threads} --time=20:00 --mem={resources.mem_mb} --export=all -o ./logs_cluster/slurm-%j.out -e ./logs_cluster/slurm-%j.err " \
  --show-failed-logs \
  --rerun-incomplete  \
  --rerun-triggers mtime 
```


# TODO
* Make an extract tool
# Workflow

![alt text](./images/dag.png)