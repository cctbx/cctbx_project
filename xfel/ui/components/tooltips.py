from __future__ import absolute_import, division, print_function

from xfel.ui import master_phil_scope

def get_help(path):
  return master_phil_scope.get(path).objects[0].help

tooltips = {
  'db_cred_ctr': get_help('experiment_tag'),
  'db_cred_btn_big': 'Set up database credentials the GUI will use to connect with',
}

def setup_tooltip(obj):
  obj.SetToolTip(tooltips.get(obj.Name))

