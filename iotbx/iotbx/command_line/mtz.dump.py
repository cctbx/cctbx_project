from iotbx import mtz
import iotbx.mtz.wrapper
from iotbx.option_parser import iotbx_option_parser
import sys, os

def process(file_name):
  print "Processing:", file_name
  mtz.wrapper.object(file_name=file_name).show_summary()
  print
  sys.stdout.flush()

def walk_callback(arg, top, names):
  for name in names:
    if (not name.lower().endswith(".mtz")): continue
    file_name = os.path.normpath(os.path.join(top, name))
    process(file_name=file_name)

def run():
  command_line = (iotbx_option_parser(
    usage="iotbx.mtz.dump [options] file_name [...]")
    .option(None, "--walk",
      action="store",
      type="string",
      dest="walk",
      metavar="ROOT_DIR",
      help="Find and process all MTZ files under ROOT_DIR")
  ).process(args=sys.argv[1:])
  for file_name in command_line.args:
    process(file_name=file_name)
  if (command_line.options.walk is not None):
    os.path.walk(top=command_line.options.walk, func=walk_callback, arg=None)

if (__name__ == "__main__"):
  run()
