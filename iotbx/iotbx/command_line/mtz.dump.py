from iotbx import mtz
from iotbx.option_parser import iotbx_option_parser
import sys, os

def process(file_name, show_column_data, show_batches):
  print "Processing:", file_name
  mtz_object = mtz.object(file_name=file_name)
  mtz_object.show_summary()
  print
  if (show_column_data):
    mtz_object.show_column_data()
    print
  if (show_batches):
    for batch in mtz_object.batches():
      batch.show()
      print "-" * 79
    print
  sys.stdout.flush()

def walk_callback(arg, top, names):
  for name in names:
    if (not name.lower().endswith(".mtz")): continue
    file_name = os.path.normpath(os.path.join(top, name))
    process(
      file_name=file_name,
      show_column_data=arg.show_column_data,
      show_batches=arg.show_batches)

def run():
  command_line = (iotbx_option_parser(
    usage="iotbx.mtz.dump [options] file_name [...]")
    .option("-v", "--verbose",
      action="store_true",
      default=False,
      dest="verbose",
      help="Enable CMTZ library messages.")
    .option(None, "--show_column_data",
      action="store_true",
      dest="show_column_data")
    .option(None, "--show_batches",
      action="store_true",
      dest="show_batches")
    .option(None, "--walk",
      action="store",
      type="string",
      dest="walk",
      metavar="ROOT_DIR",
      help="Find and process all MTZ files under ROOT_DIR")
  ).process(args=sys.argv[1:])
  if (command_line.options.verbose):
    mtz.ccp4_liberr_verbosity(1)
  for file_name in command_line.args:
    process(
      file_name=file_name,
      show_column_data=command_line.options.show_column_data,
      show_batches=command_line.options.show_batches)
  if (command_line.options.walk is not None):
    os.path.walk(
      top=command_line.options.walk,
      func=walk_callback,
      arg=command_line.options)

if (__name__ == "__main__"):
  run()
