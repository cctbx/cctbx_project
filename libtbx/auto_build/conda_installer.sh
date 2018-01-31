#!/bin/bash

function install_pkgs_conda() {
  pkgs=($@)
  for pkg in ${pkgs[@]}
  do
    echo 'PACKAGE BEING INSTALLED = ',$pkg
    conda install -y $pkg
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
while test $# -gt 0; do
  case "$1" in 
    --h|--help)
      echo "Script for installing miniconda and software using conda in cctbx"
      echo "options are"
      echo "--install-miniconda"
      echo "--install-packages {package list}"
      break
      ;;

    --install-miniconda)
      echo "install Miniconda"
      wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh --no-check-certificate
      chmod +x Miniconda2-latest-Linux-x86_64.sh
      ./Miniconda2-latest-Linux-x86_64.sh -b -p $PWD/base 
      break
      ;;

    --install-packages)
      echo "install packages using conda"
      shift
      base_pkgs=($@)
      echo ${base_pkgs[@]}
      export PATH="$PWD/../newconda/bin:$PATH"
      source activate base
      conda clean --index-cache
      install_pkgs_conda ${base_pkgs[@]} 
      break
      ;;
  esac
done
#wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh --no-check-certificate
#chmod +x Miniconda2-latest-Linux-x86_64.sh
#./Miniconda2-latest-Linux-x86_64.sh -b -p $PWD/base 
#base_pkgs=($@)
#export PATH="$PWD/base/bin:$PATH"
#source activate
#install_pkgs_conda ${base_pkgs[@]} 

