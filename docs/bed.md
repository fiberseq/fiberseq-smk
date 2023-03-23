# Generating bed files from fiberseq bam(s)
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
# aligned pseudo-bed files
results/test/density/test.CpG.bed.gz
results/test/density/test.CpG.reference.bed.gz
results/test/density/test.dinuc.bed.gz
results/test/density/test.msp.bed.gz
# aligned big bed file
results/test/bigbed/test.aligned.nuc.bed.bb
results/test/bigbed/test.aligned.m6a.bed.bb
results/test/bigbed/test.aligned.cpg.bed.bb
results/test/bigbed/test.aligned.msp.bed.bb
```
</td>
</tr>
</table>

All the bed output files are [bed12 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). If the bed output file name includes `aligned`, its bed coordinates are in reference space, conversely if it begins with `unaligned` then all the coordinates are with respect to the unaligned fiberseq read. For both `aligned` and `unaligned` there are four types of bed files:

- `cpg`: Uses the block starts in bed12 format to describe the positions of `5mC` calls from `primrose`
- `m6a`: Uses the block starts in bed12 format to describe the positions of `m6a` calls from the pipeline
- `msp`: Uses the block starts and block lengths in bed12 format to describe the start and end of Methylation Sensitive Patches (MSPs)
- `nuc`: Uses the block starts and block lengths in bed12 format to describe the start and end of nucleosomes

For all records in all the bed files the first and last block position do not reflect real data, they are only there because bed12 format requires the first and last block to match the first and last position of the entire read.

For the 4 pseudo-bed output files, the first 3 fields are consistent with the bed12 format.

- `CpG.reference.bed.gz`: demarcates CG dinucleotides in the reference genome
- `CpG.bed.gz`: gives the fraction of reads at each CG location that indicate methylation
- `dinuc.bed.gz`: gives the proportion of di-nucleosome calls (300-1000 bp) relative to all nucleosome calls (80-1000 bp) within a +/- 100 bp window, centered every 50 bp across the genome
- `msp.bed.gz`: gives the proportion of MSPs that are 150+ bp relative to the total number of MSPs within a +/- 100 bp window, centered every 50 bp across the genome
