from __future__ import absolute_import, division, print_function

from xfel.ui import master_phil_scope
try:
  from dxtbx.format.FormatXTC import locator_scope
except (ImportError, TypeError):
  locator_scope = None
from xfel.ui.command_line.plot_run_stats import phil_scope as rs_scope

def get_help(path, scope = master_phil_scope):
  return scope.get(path).objects[0].help

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
  'phenix_script_ctr': get_help('mp.phenix_script'),
  'htcondor_executable_path_ctr': get_help('mp.htcondor.executable_path'),
  'htcondor_filesystemdomain_ctr': get_help('mp.htcondor.filesystemdomain'),
  'nnodes_index_ctr': get_help('mp.nnodes_index'),
  'nnodes_tder_ctr': get_help('mp.nnodes_tder'),
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
  # DBCredentialsDialog
  'db_host_ctr': get_help('db.host'),
  'db_port_ctr': get_help('db.port'),
  'db_name_ctr': get_help('db.name'),
  'db_user_ctr': get_help('db.user'),
  'db_password_ctr': get_help('db.password'),
  'web_location_ctr': get_help('facility.lcls.web.location'),
  # StandaloneOptions
  'data_dir_ctr': get_help('facility.standalone.data_dir'),
  'monitor_for': get_help('facility.standalone.monitor_for'),
  'folders_options': get_help('facility.standalone.folders.method'),
  'n_files_needed_ctr': get_help('facility.standalone.folders.n_files_needed'),
  'last_modified_ctr': get_help('facility.standalone.files.last_modified'),
  'minimum_file_size_ctr': get_help('facility.standalone.files.minimum_file_size'),
  'template_ctr': get_help('facility.standalone.template'),
  # Main GUI
  'btn_persistent_tags': 'Auto-tag new runs as they arrive',
  'btn_manage_tags': 'Add/rename/delete tags',
  'btn_view_phil': 'View trial parameters',
  'rs_d_min': get_help('d_min', rs_scope),
  'rs_multiples': 'Number of multiple lattices before a hit is counted as a multiple hit',
  'rs_ratio': 'Ratio of 2θ high to 2θ low needed for an image to be a solvent hit',
  'rs_n_strong': get_help('n_strong_cutoff', rs_scope),
  'rs_isigi': get_help('i_sigi_cutoff', rs_scope),
  'rs_n_dump': 'Number of images to convert to cbf and then display',
  'uc_selection_type': 'Union: include runs matching any of these tags\n' + \
                       'Intersection: include runs matching all of these tags',
  # Trial dialog
  'trial_throttle_ctr': 'Percent of images (events) to process',
  'trial_num_bins_ctr': 'Used for logging only',
  'trial_d_min_ctr': 'Used for logging only', # XXX doesn't appear
  # Run group dialog
  'rg_end_type': 'Auto add runs: new data will be added to this block as it arrives\nSpecify end run: set the last run for this block explicitly.',
  'rg_address_ctr': 'Detector address in XTC stream (use detnames to list available detectors)',
  'rg_beam_xyz': 'Beam center in pixels, and detector distance in mm (overridden by the phil parameter input.reference_geometry)',
  'rg_bin_nrg_gain_binning': 'Rayonix binning (2, 3, 4, etc.)',
  'rg_bin_nrg_gain_energy': 'Energy override for all images (eV)',
  'rg_wavelength_offset': 'Offset applied to wavelength of each image (Å)',
  'rg_spectrum_calibration': get_help('spectrum_eV_per_pixel', locator_scope) if locator_scope else '',
  'rg_energy_ctr': 'Energy override for all images (eV)',
  'rg_two_thetas': 'Two 2θ values (deg). The ratio of high/low is used to check for presence of solvent on each image. ' + \
                   'Defaults are the water ring and a low resolution ring',
}

def setup_tooltip(obj):
  obj.SetToolTip(tooltips.get(obj.Name))

