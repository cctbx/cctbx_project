"""Summarize contents of a reflection file or show each reflection"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.mtz.dump
# LIBTBX_SET_DISPATCHER_NAME iotbx.mtz.dump

from iotbx import mtz
from iotbx.option_parser import option_parser
import sys, os

def process(file_name, show_column_data, column_data_format, show_batches):
  if (show_column_data):
    column_data_format = mtz.tidy_show_column_data_format_keyword(
      input=column_data_format)
  else:
    column_data_format = ""
  if (column_data_format != "spreadsheet"):
    print("Processing:", file_name)
  mtz_object = mtz.object(file_name=file_name)
  if (column_data_format != "spreadsheet"):
    mtz_object.show_summary()
    print()
  if (show_column_data):
    mtz_object.show_column_data(format=column_data_format)
    if (column_data_format != "spreadsheet"):
      print()
  if (show_batches):
    for batch in mtz_object.batches():
      batch.show()
      print("-" * 79)
    print()
  sys.stdout.flush()

def walk_callback(arg, top, names):
  for name in names:
    if (not name.lower().endswith(".mtz")): continue
    file_name = os.path.normpath(os.path.join(top, name))
    process(
      file_name=file_name,
      show_column_data=arg.show_column_data,
      show_batches=arg.show_batches)

def run(args, command_name="phenix.mtz.dump"):
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage=command_name+" [options] file_name [...]")
    .option("-v", "--verbose",
      action="store_true",
      default=False,
      help="Enable CMTZ library messages.")
    .option("-c", "--show_column_data",
      action="store_true")
    .option("-f", "--column_data_format",
      action="store",
      type="string",
      metavar="KEYWORD",
      help="Valid keywords are: %s."
             % ", ".join(mtz.show_column_data_format_keywords)
          +" Human readable is the default. The format keywords can be"
          +" abbreviated (e.g. -f s).")
    .option("-b", "--show_batches",
      action="store_true")
    .option(None, "--walk",
      action="store",
      type="string",
      metavar="ROOT_DIR",
      help="Find and process all MTZ files under ROOT_DIR")
  ).process(args=args)
  if (len(command_line.args) == 0):
    print(command_line.parser.format_help())
    return
  if (command_line.options.verbose):
    mtz.ccp4_liberr_verbosity(1)
  for file_name in command_line.args:
    process(
      file_name=file_name,
      show_column_data=command_line.options.show_column_data,
      column_data_format=command_line.options.column_data_format,
      show_batches=command_line.options.show_batches)
  if (command_line.options.walk is not None):
    os.path.walk(
      top=command_line.options.walk,
      func=walk_callback,
      arg=command_line.options)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
