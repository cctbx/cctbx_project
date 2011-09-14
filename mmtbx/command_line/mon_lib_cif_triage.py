from mmtbx.monomer_library import cif_triage
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, Usage
import libtbx.load_env
import sys

def run(args, command_name=libtbx.env.dispatcher_name):
  if (len(args) == 0):
    raise Usage("%s cif [...]" % command_name)
  for file_name in args:
    obj_count = cif_triage.check_comp(file_name=file_name)
    if (obj_count == 0):
      raise Sorry("No data found in file: %s" % show_string(file_name))
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
