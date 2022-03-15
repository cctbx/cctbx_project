#!/bin/bash -l
#SBATCH --qos=<queue>
#SBATCH --job-name=<jobname>
#SBATCH --reservation=<reservation>
#SBATCH --time=<walltime>
#SBATCH --ntasks=2
#SBATCH --constraint=<constraint>
#SBATCH --mail-type=NONE
#SBATCH -A <project>
#SBATCH -o <out_log>
#SBATCH -e <err_log>
#DW jobdw capacity=10GB access_mode=striped type=scratch
#DW stage_out source=$DW_JOB_STRIPED/stdout destination=<output_dir>/stdout type=directory
# base directory

# submit jobs
mkdir ${DW_JOB_STRIPED}/stdout    #DW <-Tagged so we can delete this line if not using DW
cd <output_dir>

echo -n "Starting job @ t="; date +%s
srun -n 1 <srun_script>
echo -n "Finished job @ t="; date +%s

