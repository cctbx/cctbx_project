from __future__ import absolute_import, division, print_function

from xfel.ui import master_phil_scope

def get_help(path):
  return master_phil_scope.get(path).objects[0].help

tooltips = {
  # Settings dialog
  'db_cred_ctr': get_help('experiment_tag'),
  'db_cred_btn_big': 'Set up database credentials the GUI will use to connect with',
  'facility_ctr': get_help('facility.name'),
  'btn_facility_options': 'Facility specific options',
  'experiment_ctr': get_help('facility.lcls.experiment'),
  'output_ctr': get_help('output_folder'),
  'advanced': 'Multiprocessing, queueing, and other options',
  # Advanced settings dialog
  'mp_option_ctr': get_help('mp.method'),
  'queue_ctr': get_help('mp.queue'),
  'nproc_ctr': get_help('mp.nproc'),
  'nnodes_ctr': get_help('mp.nnodes'),
  'nproc_per_node': get_help('mp.nproc_per_node'),
  'wall_time_ctr': get_help('mp.wall_time'),
  'mpi_command_ctr': get_help('mp.mpi_command'),
  'env_script_ctr': get_help('mp.env_script'),
  'htcondor_executable_path_ctr': get_help('mp.htcondor.executable_path'),
  'htcondor_filesystemdomain_ctr': get_help('mp.htcondor.filesystemdomain'),
  'nnodes_index_ctr': get_help('mp.nnodes_index'),
  'nnodes_scale_ctr': get_help('mp.nnodes_scale'),
  'nnodes_merge_ctr': get_help('mp.nnodes_merge'),
  'extra_options': get_help('mp.extra_options'),
  'shifter_image_ctr': get_help('mp.shifter.shifter_image'),
  'shifter_srun_template_ctr': get_help('mp.shifter.srun_script_template'),
  'shifter_sbatch_template_ctr': get_help('mp.shifter.sbatch_script_template'),
  'shifter_jobname_ctr': get_help('mp.shifter.jobname'),
  'shifter_project_ctr': get_help('mp.shifter.project'),
  'shifter_reservation_ctr': get_help('mp.shifter.reservation'),
  'shifter_constraint_ctr': get_help('mp.shifter.constraint'),
  'staging_ctr': get_help('mp.shifter.staging'),
  'back_end_ctr': get_help('dispatcher'),
}

def setup_tooltip(obj):
  obj.SetToolTip(tooltips.get(obj.Name))

