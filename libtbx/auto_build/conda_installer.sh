#!/bin/bash

function install_pkgs_conda() {
  pkgs=($@)
  for pkg in ${pkgs[@]}
  do
    echo 'PACKAGE BEING INSTALLED = ',$pkg
    conda install -y --dry-run $pkg
  done
}

function install_pkg_base() {
  pkg=$1
  echo 'BASE PACKAGE BEING INSTALLED =' $pkg
#  python modules/cctbx_project/libtbx/auto_build/install_base_packages.py --with-python=base/bin/python $pkg
}

############################################################
# Doing the base step with a conda environment which provides all dependecies
# The conda installation directory will be called base for now
echo $PWD '########### ---- installing base packages with conda ---- ######### '
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh --no-check-certificate
chmod +x Miniconda2-latest-Linux-x86_64.sh
./Miniconda2-latest-Linux-x86_64.sh -b -p $PWD/base 
base_pkgs=($@)
export PATH="$PWD/base/bin:$PATH"
source activate
install_pkgs_conda ${base_pkgs[@]} 

