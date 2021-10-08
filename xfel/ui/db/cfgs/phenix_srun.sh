#!/bin/bash

#cctbx
source /global/common/software/cctbx/phenix-1.19.2-4158/build/setpaths.sh

#for experiment database
export SIT_DATA=/global/common/software/lcls/psdm/data

#for psana
export SIT_PSDM_DATA=/global/cscratch1/sd/psdatmgr/data/psdm

#needed to open h5 files from compute nodes
export HDF5_USE_FILE_LOCKING=FALSE

# run
<command>
