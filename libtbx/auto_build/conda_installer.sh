#!/bin/bash

function install_pkgs_conda() {
  pkgs=($@)
  for pkg in ${pkgs[@]}
  do
    echo 'PACKAGE BEING INSTALLED = ',$pkg
    return 22
    conda install -y $pkg
  done
}

function install_psanaconda_lcls_channel() {
  rhel_version=($@)
  echo 'PSANA CONDA CHANNEL BEING INSTALLED = lcls-rhel'${rhel_version}
  conda install -y --channel lcls-rhel${rhel_version} psana-conda | grep 'PackagesNotFoundError'
  ret_code=$?
  echo $ret_code
  if [ $ret_code -gt 0 ]; then
    exit $ret_code
  fi
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
      echo "--install-psanaconda-lcls {Specify Red Hat Linux Version: 5,6 or 7}"
      break
      ;;

    --install-miniconda)
      echo "install Miniconda"
      wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh --no-check-certificate
      chmod +x Miniconda2-latest-Linux-x86_64.sh
      ./Miniconda2-latest-Linux-x86_64.sh -b -p $PWD/newconda 
      export PATH="$PWD/newconda/bin:$PATH"
      # Make sure no old installation of conda and softlinks exist
      rm -r base
      conda env remove -y --name myEnv
      conda create -y -n myEnv
      ln -s $PWD/newconda/envs/myEnv base
      break
      ;;

    --install-packages)
      echo "install packages using conda"
      shift
      base_pkgs=($@)
      echo ${base_pkgs[@]}
      export PATH="$PWD/newconda/bin:$PATH"
      echo $PATH
      source activate myEnv
      conda clean --index-cache
      install_pkgs_conda ${base_pkgs[@]} 
      break
      ;;
    --install-psanaconda-lcls)
      echo "install packages through psana-conda lcls channel. Note it doesn't contain all base packages and should only be used for XFEL on a Linux system"
      shift
#      rhel_version=`cat /etc/redhat-release | sed 's/.*release \([0-9]*\).*/\1/'`
      rhel_version=`cat /etc/redhat-release | sed 's/.*release \([0-9]*\).*/\1/'`
      export PATH="$PWD/newconda/bin:$PATH"
      echo $PATH
      source activate myEnv
      conda clean --index-cache
      install_psanaconda_lcls_channel ${rhel_version} 
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

