# Snakemake workflow: `fiberseq-smk`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/StergachisLab/fiberseq-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/StergachisLab/fiberseq-smk/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for making *fiberseq* calls.
# Usage 

```
snakemake --cores 10 -p --set-scatter chunks=2  --show-failed-logs --configfile .test/config.yml -n
```
# TODO
* Make an extract tool
# Workflow

![alt text](./images/dag.png)