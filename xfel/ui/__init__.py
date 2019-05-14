from __future__ import absolute_import, division, print_function
import os
from iotbx.phil import parse
from libtbx.utils import Sorry

master_phil_str = """
dispatcher = cctbx.xfel.xtc_process
  .type = str
  .help = Which program to run. cxi.xtc_process is for module only based processing, \
          such as mod_hitfind. cctbx.xfel.xtc_process uses the DIALS back end.
dry_run = False
  .type = bool
  .help = If True, the program will create the trial directory but not submit the job, \
          and will show the command that would have been executed.
facility {
  name = *lcls standalone
    .type = choice
    .help = Facility for the XFEL gui. LCLS configures the GUI to use LCLS services \
            for data monitoring, job submission, and so forth. Standlone runs the \
            GUI for all other data sources.
  lcls {
    experiment = ""
      .type = str
      .help = Experiment name, eg cxid9114
    web {
      user = ""
        .type = str
        .help = Username for LCLS run database web service
      password = ""
        .type = str
        .help = Web password. Will be cached in plain text!
      enforce80 = False
        .type = bool
        .help = report only on stream 81, FEE spectrometer
      enforce81 = False
        .type = bool
        .help = report only on stream 81, FEE spectrometer
    }
    use_ffb = False
      .type = bool
      .help = Run on the ffb if possible. Only for active users!
    dump_shots = False
      .type = bool
      .help = Write images to disk whether they index or not. \
              Helpful for tuning spotfinding and indexing parameters, and necessary \
              for the "Should have indexed" feature of the Run Stats tab.
  }
  standalone {
    data_dir = None
      .type = path
      .help = Folder to monitor for new data
    monitor_for = *files folders
      .type = choice
      .help = Whether to monitor for new files or new folders with files in them
    template = None
      .type = str
      .help = File matching pattern for new data. Example: *.h5
    composite_files = True
      .type = bool
      .help = If True, files are submitted as individual runs. Otherwise, groups \
              of files are submitted as single runs
    }
}
output_folder = ""
  .type = path
  .help = Processing results will go in this folder
average_raw_data = False
  .type = bool
  .help = If True, don't use any psana corrections (dark, common mode, etc.)

include scope xfel.command_line.cxi_mpi_submit.mp_phil_scope
"""
db_phil_str = """
experiment_tag = None
  .type = str
  .help = User defined tag to describe the set of trials being performed. All database tables will \
          be pre-pended with this string
db {
  host = psdb-user.slac.stanford.edu
    .type = str
    .help = Host name for mysql databse server
  name = ""
    .type = str
    .help = Database name
  user = ""
    .type=str
    .help = Database user name
  password = ""
    .type = str
    .help = Database password. Will be cached as plain text!
  verbose = False
    .type = bool
    .expert_level = 2
    .help = Print to the terminal all database queries
  port = 3306
    .type = int
    .help = Port number to be used for connection
}
"""
master_phil_scope = parse(master_phil_str + db_phil_str, process_includes=True)

settings_dir = os.path.join(os.path.expanduser('~'), '.cctbx.xfel')
settings_file = os.path.join(settings_dir, 'settings.phil')

known_dials_dispatchers = {
  'cctbx.xfel.xtc_process': 'xfel.command_line.xtc_process',
  'cctbx.xfel.process': 'xfel.command_line.xfel_process',
  'dials.stills_process': 'dials.command_line.stills_process',
  'cctbx.xfel.small_cell_process': 'xfel.small_cell.command_line.xfel_small_cell_process',
}

def load_phil_scope_from_dispatcher(dispatcher):
  import importlib
  try:
    phil_scope = importlib.import_module(known_dials_dispatchers[dispatcher]).phil_scope
  except KeyError:
    import imp
    mod = imp.load_source('module', dispatcher)
    phil_scope = mod.phil_scope
  return phil_scope

def load_cached_settings():
  if os.path.exists(settings_file):
    user_phil = parse(file_name = settings_file)
    params, unused = master_phil_scope.fetch(source = user_phil, track_unused_definitions=True)
    params = params.extract()
    if len(unused) > 0:
      from .phil_patch import sync_phil
      sync_phil(params, unused)
    return params
  else:
    return master_phil_scope.extract()

def save_cached_settings(params):
  if not os.path.exists(settings_dir):
    os.makedirs(settings_dir)

  working_phil = master_phil_scope.format(python_object = params)
  diff_phil = master_phil_scope.fetch_diff(source = working_phil)

  try:
    f = open(settings_file.encode('utf8'), 'wb')
    f.write(diff_phil.as_str())
    f.close()
  except IOError:
    raise Sorry('Unable to write %s.' % settings_file)
