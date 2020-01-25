=== Instructions for merging multiple specific isoforms within a serial crystallography dataset ===

The purpose of this workflow is to separately group crystal isoforms with very similar unit cells.
This is distinct from tools like dials.cosym, which deals with the assignment of Laue symmetry,
and with the resolution of indexing ambiguity for polar space groups.

General steps in the pipeline are:

1) Process the serial crystallography diffraction patterns with dials.stills_process, 
producing files that contain refined experiment models and integrated Bragg spots.

2) Run a short version of the program cctbx.xfel.merge with two steps only: 
reading the data in, and producing a text list of the unit cells and space groups of
all integrated lattices.

cores=32

export effective_params="dispatch.step_list=input \
dispatch.step_list=tdata \
input.path=${DIALS_OUTPUT}/r*/*/out \
input.experiments_suffix=_integrated_experiments.json \
input.reflections_suffix=_integrated.pickle \
input.parallel_file_load.method=uniform \
tdata.output_path=${TRIAL}_cells \
output.prefix=${TRIAL} \
output.output_dir=${OUT_DIR}/${TRIAL}/out 
output.tmp_dir=${OUT_DIR}/${TRIAL}/tmp \
output.do_timing=True \
output.log_level=0"
 
mpiexec -n $cores cctbx.xfel.merge  ${effective_params}
sample output:
78.643232 78.643232 38.831412 90.000000 90.000000 90.000000 P4/mmm
78.889991 78.889991 38.179862 90.000000 90.000000 90.000000 P4/mmm

