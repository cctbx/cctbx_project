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
wget -r http://cci.lbl.gov/dials/installers/current/ -A *linux* -l1
tar -xvf cci.lbl.gov/dials/installers/current/*linux*
DEV=`ls cci.lbl.gov/dials/installers/current | xargs python -c "import sys; print (sys.argv[1].split('-')[3])"`
cd dials-installer-*
./install --prefix=$WORKING
cd -

# Download and install conda and mysql
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $WORKING/miniconda3
source miniconda3/etc/profile.d/conda.sh
conda activate $WORKING/dials-dev-$DEV/conda_base
conda install -y mysql mysql-python -c conda-forge --no-deps

# Clean up
rm -rf cci.lbl.gov dials-installer-*

# Create setup file
echo "\
#!/bin/bash
source $WORKING/miniconda3/etc/profile.d/conda.sh
conda activate $WORKING/dials-dev-$DEV/conda_base
source $WORKING/dials-dev-$DEV/build/setpaths.sh" > $WORKING/setup.sh

echo "Download complete. Add the sources to your path with 'source $WORKING/setup.sh'"
