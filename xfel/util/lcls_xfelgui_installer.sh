#!/bin/bash

# Script to download latest DIALS sources in the provided folder, install
# conda, configure the psana environment using the xfel_dependencies
# packages and MySQL, and compile the sources.

if [ $# -lt 1 ]
 then
   echo "Not enough arguments supplied"
   echo "Please specify the following arguments to correctly run the installation:
   {WORKING}"
   exit
fi

WORKING=$1
NPROC=$2
if [ -z $NPROC ]; then
  NPROC=1
fi

cd $WORKING

# Download DIALS sources
mkdir dialsBuild; cd dialsBuild
wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
python bootstrap.py --builder=dials hot update

# Download and install dependencies using conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $WORKING/dialsBuild/miniconda3
# psana1 needs some discontinued packages
echo "\
restore_free_channel: true
" > $WORKING/dialsBuild/miniconda3/.condarc
source $WORKING/dialsBuild/miniconda3/etc/profile.d/conda.sh
conda env create -f $WORKING/dialsBuild/modules/cctbx_project/xfel/conda_envs/psana_environment.yml
conda install -y -c conda-forge -n psana_env --no-deps gemmi
conda activate psana_env

# Build DIALS sources
python bootstrap.py --builder=dials build --nproc $NPROC --use-conda $WORKING/dialsBuild/miniconda3/envs/psana_env

# Create setup file
echo "\
#!/bin/bash
source $WORKING/dialsBuild/miniconda3/etc/profile.d/conda.sh
conda activate psana_env
source $WORKING/dialsBuild/build/setpaths.sh" > $WORKING/setup.sh

echo "Download complete. Add the sources to your path with 'source $WORKING/setup.sh'"
