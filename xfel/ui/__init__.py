from __future__ import division
import os
from iotbx.phil import parse
from libtbx.utils import Sorry

master_phil_str = """
experiment = ""
  .type = str
experiment_tag = ""
  .type = str
db {
  host = psdb-user.slac.stanford.edu
    .type=str
  name = ""
    .type = str
  user = ""
    .type=str
  password = ""
    .type = str
}
output_folder = ""
  .type = path
web {
  user = None
    .type = str
  password = None
    .type = str
}
"""
master_phil_scope = parse(master_phil_str, process_includes=True)

settings_file = os.path.join(os.path.expanduser('~'), '.cctbx.xfel/settings.phil')

def load_cached_settings():
  if os.path.exists(settings_file):
    user_phil = parse(file_name = settings_file)
    return master_phil_scope.fetch(source = user_phil).extract()
  else:
    return master_phil_scope.extract()

def save_cached_settings(params):
  settings_dir = os.path.dirname(settings_file)
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
