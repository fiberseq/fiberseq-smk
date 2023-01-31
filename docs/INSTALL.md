# Install

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
First make sure you channel priorities are set correctly according to the [bioconda documentation](https://bioconda.github.io/#usage):
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Then install the conda environment used by the workflow:
```bash
git clone https://github.com/StergachisLab/fiberseq-smk
cd fiberseq-smk
conda create -n fiberseq-smk
mamba env update -n fiberseq-smk --file workflow/envs/env.yml 
conda activate fiberseq-smk
```
