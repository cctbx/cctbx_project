from __future__ import absolute_import, division, print_function

tooltips = {
  'db_cred_ctr': 'Identifier for these processing results. Similar to a project name.',
  'db_cred_btn_big': 'Set up database credentials the GUI will use to connect with',
}

def setup_tooltip(obj):
  obj.SetToolTip(tooltips.get(obj.Name))

