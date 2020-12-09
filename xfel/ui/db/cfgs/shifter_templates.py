from __future__ import division
sbatch_template = \
"""#!/bin/bash -l
#SBATCH --qos=<queue>
#SBATCH --job-name=<jobname>
#SBATCH --reservation=<reservation>
#SBATCH --time=<walltime>
#SBATCH --nodes=<nnodes>
#SBATCH --tasks-per-node=<nproc_per_node>
#SBATCH --constraint=<constraint>
#SBATCH --image=<shifter_image>
#SBATCH --mail-type=NONE
#SBATCH -A <project>
#SBATCH -o <out_log>
#SBATCH -e <err_log>
#DW jobdw capacity=10GB access_mode=striped type=scratch
#DW stage_out source=$DW_JOB_STRIPED/stdout destination=<output_dir>/stdout type=directory
# base directory

# submit jobs
mkdir ${DW_JOB_STRIPED}/stdout    #DW <-Tagged so we can delete this line if not using DW

echo -n "Starting job @ t="; date +%s
srun shifter <srun_script>
echo -n "Finished job @ t="; date +%s
"""

srun_template = \
"""#!/bin/bash

#cctbx
source /img/activate.sh

#for experiment database
export SIT_DATA=/global/cscratch1/sd/blaschke/lv95/psdm/data

#for psana
export SIT_PSDM_DATA=/global/cscratch1/sd/psdatmgr/data/psdm

#needed to open h5 files from compute nodes
export HDF5_USE_FILE_LOCKING=FALSE

# run
<command>
"""
