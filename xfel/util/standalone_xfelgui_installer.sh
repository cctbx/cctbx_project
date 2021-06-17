#!/bin/bash

# Script to download latest DIALS bundle, install it in the provided folder, and
# add MySQL.

if [ $# -lt 1 ]
 then
   echo "Not enough arguments supplied"
   echo "Please specify the following arguments to correctly run the installation:
   {WORKING}"
   exit
fi

# Make sure an absolute path was provided. If not, print a warning
target_file=$(basename $1)
pushd $(dirname $1) >/dev/null
target_dir_phys=$(pwd -P)
popd >/dev/null

if [ "$target_dir_phys/$target_file" != "$1" ]
then
  echo "Warning: Please give an absolute path to the installation directory."
  read -p "Proceed anyway? y/[n] " run_anyway
  if [ "$run_anyway" != "y" ]; then exit; fi
fi

WORKING=$1
cd $WORKING

# Download and install the latest build of DIALS
wget $(curl -s https://api.github.com/repos/dials/dials/releases/latest | grep browser_download_url |grep "linux.*\.xz" |tail -1 |cut -d'"' -f 4)
tar -xvf dials*.tar.xz
cd dials-installer
# Replace the conda-forge iota with the full version from github
git clone git@github.com:ssrl-px/iota modules/iota
rm -r conda_base/lib/python3*/site-packages/iota*
./install --prefix=$WORKING
cd -
rm -rf dials-installer dials*.tar.xz

DIALS_DIR=$(ls -d dials*)

# Download and install conda and mysql
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $WORKING/miniconda3
source miniconda3/etc/profile.d/conda.sh
conda activate $WORKING/$DIALS_DIR/conda_base
conda install -y mysql mysql-connector-python mysqlclient -c conda-forge

# Clean up
rm Miniconda3-latest-Linux-x86_64.sh

# Create setup file
echo "\
#!/bin/bash
source $WORKING/miniconda3/etc/profile.d/conda.sh
conda activate $WORKING/$DIALS_DIR/conda_base
source $WORKING/$DIALS_DIR/build/setpaths.sh" > $WORKING/setup.sh

echo "Download complete. Add the sources to your path with 'source $WORKING/setup.sh'"
