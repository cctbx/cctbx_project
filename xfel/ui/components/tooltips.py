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
  'db_cred_btn_big': 'Set up database credentials the GUI will use to connect with. '
                     'Disabled while the GUI is running — restart to change connection settings.',
  'load_project': 'Load a saved project bundle (experiment tag, database, output folder, and '
                  'multiprocessing settings) from ~/.cctbx.xfel/settings_<name>.phil. '
                  'While the GUI is running, changing the database connection or experiment tag '
                  'is not allowed; a full restart is required for those changes.',
  'save_project': 'Save the current settings as a named project bundle to '
                  '~/.cctbx.xfel/settings_<name>.phil.',
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
  # StartDBDialog / DBCredentialsDialog
  'basedir': 'Directory where the MySQL server data files will be stored. '
             'Defaults to <output_folder>/MySql.',
  'db_root_password': 'Root password for the MySQL server. Only required when '
                      'starting a local server via "Start DB Server".',
  'start_db': 'Launch a local MySQL server process using the directory and '
              'password specified above.',
  'db_OK': 'Confirm database credentials and close this dialog.',
  # Ensemble refinement / MergingStats dialogs
  'nnodes_tder': 'Number of nodes for Time-Dependent Ensemble Refinement.',
  # EnergyDialog
  'skip_images': 'Number of images to skip at the start of the dataset.',
  'num_images': 'Maximum number of frames to average.',
  # TrialDialog — trial meta
  'trial_info': 'Trial number for this processing run. Use "Import PHIL" to load '
                'parameters from a file, or "Edit PHIL" to modify them directly.',
  'trial_comment': 'Optional note attached to this trial (stored in the database).',
  'copy_runblocks': 'Initialise this trial\'s run blocks by copying them from '
                    'an existing trial.',
  # TrialDialog — overall / spotfinding
  'min_spots': 'Minimum number of strong spots required to attempt indexing an image.',
  'min_spot_size': 'Minimum number of pixels a connected region must have to be '
                   'counted as a spot (filters single-pixel noise).',
  'max_spot_size': 'Maximum number of pixels a spot may have (filters ice rings '
                   'and bad pixels).',
  'sigma_background': 'Local background standard-deviation multiplier for the '
                      'dispersion spotfinder. Lower values are more sensitive.',
  'sigma_strong': 'Minimum signal-to-noise ratio above background for a pixel to '
                  'be counted as part of a spot.',
  'global_threshold': 'Absolute intensity threshold; pixels below this value are '
                      'never counted as signal regardless of local background.',
  'gain': 'Detector gain in ADU/photon. Used to convert pixel values to photon '
          'counts for the dispersion algorithm.',
  'kernel_size': 'Size of the local background estimation kernel in pixels '
                 '(the neighbourhood used to estimate mean and standard deviation).',
  'threshold_algorithm': 'Spotfinding threshold algorithm. "dispersion" is the '
                         'standard method; "dispersion_extended" is more aggressive '
                         'at low resolution; "radial_profile" uses radial background subtraction.',
  # TrialDialog — indexing
  'unit_cell': 'Target unit cell parameters (a b c α β γ) to guide indexing.',
  'space_group': 'Target space group for indexing (e.g. "P 21 21 21").',
  'd_min_indexing': 'High-resolution cutoff used during indexing (Å). Reflections '
                    'beyond this limit are ignored.',
  'max_lattices': 'Maximum number of lattices to find per image. Values greater '
                  'than 1 enable multi-lattice indexing.',
  # DatasetDialog — shared settings applied to all stages
  'shared_model': 'Path to the reference model (MTZ or PDB) used for scaling and '
                  'merging. Click Browse to select a file.',
  'shared_unit_cell': 'Unit cell parameters (a b c α β γ) applied to all stages '
                      'in this dataset.',
  'shared_space_group': 'Space group applied to all stages in this dataset '
                        '(e.g. "P 21 21 21").',
  'shared_d_min': 'High-resolution cutoff for scaling and merging (Å).',
  'shared_resolution_scalar': 'Scalar multiplied by d_min to set the internal '
                              'resolution limit. Values less than 1 extend slightly '
                              'beyond d_min.',
  'shared_n_bins': 'Number of resolution shells used in statistics output tables '
                   'and plots.',
  'shared_merge_anomalous': 'When checked, Friedel mates are merged and anomalous '
                            'signal is not preserved. Uncheck for anomalous phasing.',
  # DatasetDialog scaling stage friendly controls
  'scale_min_corr': 'Minimum Pearson correlation coefficient between a lattice\'s '
                    'intensities and the reference model. Lattices below this '
                    'threshold are rejected before merging.',
  'scale_rel_tol': 'Maximum fractional difference in unit cell edge length allowed '
                   'versus the reference cell (value-mode filter). Superseded by '
                   'the cluster filter when that mode is enabled.',
  'scale_sigma': 'Significance filter: for each image, select the highest resolution bin '
                 'with I/σ above this cutoff.',
  # DatasetTab filter
  'filter': 'Filter the dataset list by name. Type any substring to narrow the list.',
  # UnitCellTab
  'uc_plot_eps': 'DBSCAN epsilon: maximum distance between two unit cells to be '
                 'considered neighbours when forming clusters. Smaller values '
                 'produce tighter, more selective clusters.',
  # Dataset scaling stage — unit-cell cluster filter
  'chk_use_cluster': 'When checked, filter lattices using a covariance model previously '
                     'computed on the Unit Cells tab. When unchecked, the relative-length '
                     'tolerance filter (above) is used instead.',
  'cluster_file': 'Covariance pickle file written by the Unit Cells tab '
                  '(output_folder/cluster/cluster_<name>.pickle). '
                  'Use Browse to select a file outside the default location.',
  'browse_cluster': 'Browse for a covariance pickle file outside the default cluster directory.',
  'cluster_component': 'Index of the Gaussian mixture component to filter on (0 = largest cluster). '
                       'The Unit Cells tab labels components in the same order.',
  'cluster_mahalanobis': 'Maximum Mahalanobis distance from the cluster centre for a lattice to be '
                         'accepted (essentially a sigma cutoff for the multivariate Gaussian). '
                         'Default is 4.0.',
  # Dataset stage checkboxes and radio groups
  'enable_chk': 'Include this stage in the dataset processing pipeline. Uncheck to skip '
                'the stage without losing its configured settings.',
  'chk_pre_split': 'Chunk the processing results locally for better chunk sizes or on the '
                   'cluster. On the cluster is faster but less accurate for counting '
                   'chunk sizes.',
  'chk_expand_nave': 'During reintegration, expand the mosaic (Nave) parameters to '
                     'catch more integrated reflections',
  'selection_type_radio': 'How to combine the selected tags: "intersection" keeps runs '
                          'matching ALL selected tags; "union" keeps runs matching ANY '
                          'selected tag.',
  'model_mode_radio': 'Whether a known reference model is supplied for scaling and merging. '
                      'Choose "No reference model" to scale and merge without one.',
}

def setup_tooltip(obj):
  # Look up the tooltip by the widget's name. Two naming conventions coexist:
  # inner widgets of composite controls are named "<panel>_ctr"/"<panel>_btn_big"
  # and register their own keys, while composite panels register a key under the
  # bare panel name. When a bare-name key matches, propagate it down to every
  # child widget so the tip shows no matter which sub-widget the pointer is over
  # (on wxGTK a tooltip set only on the parent panel is not shown over a child
  # native control). Children that already carry their own tip are left alone.
  tip = tooltips.get(obj.Name)
  if not tip:
    return
  _apply_tooltip(obj, tip, is_root=True)

def _apply_tooltip(widget, tip, is_root=False):
  if is_root or not widget.GetToolTipText():
    widget.SetToolTip(tip)
  for child in widget.GetChildren():
    _apply_tooltip(child, tip)

