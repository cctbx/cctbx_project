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
}

def setup_tooltip(obj):
  obj.SetToolTip(tooltips.get(obj.Name))

