# Install

All dependencies for the primary workflow are installed via conda. To set up your environment first make sure you channel priorities are set correctly according to the [bioconda documentation](https://bioconda.github.io/#usage):
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
Then install the conda environment used by the workflow with these commands:
```bash
git clone https://github.com/StergachisLab/fiberseq-smk
cd fiberseq-smk
conda create -n fiberseq-smk
mamba env update -n fiberseq-smk --file workflow/envs/env.yml 
conda activate fiberseq-smk
```

# Install dependencies for ipdSummary pipeline
If you want to use the `ipdSummary` pipeline you will need to install some additional dependencies available through `SMRTLINK`. 
```
# Add the SMRTlink tools to your path
export PATH=$PATH:/path/to/smrtlink/tools
```
If you are on `hyak` you can add our copy to your path by adding this to your `.bashrc`:
```bash
PATH=$PATH:/gscratch/stergachislab/install_dir/smrtlink/smrtcmds/bin/
```

