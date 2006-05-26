from iotbx import pdb
from iotbx.option_parser import iotbx_option_parser
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
import sys, os

def run(args):
  if (len(args) == 0): args = ["--help"]
  command_line = (iotbx_option_parser(
    usage="iotbx.pdb.hierarchy file...")
    .option(None, "--details",
      action="store",
      type="string",
      default=None,
      dest="details",
      help="level of detail",
      metavar="|".join(pdb.hierarchy_level_ids))
    .option(None, "--prefix",
      action="store",
      type="string",
      default="",
      dest="prefix",
      help="prefix for all output lines",
      metavar="STRING")
  ).process(args=args)
  level_id = command_line.options.details
  prefix = command_line.options.prefix
  file_names = command_line.args
  for file_name in file_names:
    if (not os.path.isfile(file_name)): continue
    execute(file_name=file_name, level_id=level_id, prefix=prefix)
    print prefix

def execute(file_name, level_id=None, prefix=""):
  try:
    pdb.show_summary(
      file_name=file_name,
      level_id=level_id,
      level_id_exception=Sorry,
      prefix=prefix)
  except KeyboardInterrupt: raise
  except Exception, e:
    print "Exception: file %s: %s: %s" % (
      show_string(file_name), e.__class__.__name__, str(e))

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
