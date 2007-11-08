# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.hierarchy

from iotbx import pdb
from iotbx.option_parser import option_parser
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
import sys, os

def run(args, command_name="phenix.pdb.hierarchy"):
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage="%s file..." % command_name)
    .option(None, "--details",
      action="store",
      type="string",
      default=None,
      help="level of detail",
      metavar="|".join(pdb.hierarchy_level_ids))
    .option(None, "--duplicate_max_show",
      action="store",
      type="int",
      default=10,
      help="maximum number of groups of duplicate atom labels to be listed",
      metavar="INT")
    .option(None, "--prefix",
      action="store",
      type="string",
      default="",
      help="prefix for all output lines",
      metavar="STRING")
  ).process(args=args)
  co = command_line.options
  for file_name in command_line.args:
    if (not os.path.isfile(file_name)): continue
    execute(
      file_name=file_name,
      level_id=co.details,
      duplicate_max_show=co.duplicate_max_show,
      prefix=co.prefix)
    print co.prefix.rstrip()

def execute(file_name, level_id=None, duplicate_max_show=10, prefix=""):
  try:
    pdb.show_summary(
      file_name=file_name,
      level_id=level_id,
      level_id_exception=Sorry,
      duplicate_max_show=duplicate_max_show,
      prefix=prefix)
  except KeyboardInterrupt: raise
  except Exception, e:
    print "Exception: file %s: %s: %s" % (
      show_string(file_name), e.__class__.__name__, str(e))

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
