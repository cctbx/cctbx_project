"""Validate an mmCIF file"""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME iotbx.cif.validate

import glob, os, sys
from six.moves import urllib
from six.moves import cStringIO as StringIO

from iotbx import cif
from iotbx.cif import validation
from iotbx.option_parser import option_parser
from libtbx.utils import time_log
import libtbx.load_env

def run(args, out=sys.stdout):
  if len(args) == 0: args = ["--help"]
  command_line = (option_parser(
                  usage="iotbx.cif.validate filepath|directory [options]")
                  .option(None, "--file_ext",
                          action="store",
                          default="cif")
                  .option(None, "--dic",
                          action="append",
                          dest="dictionaries")
                  .option(None, "--show_warnings",
                          action="store_true")
                  .option(None, "--show_timings",
                          action="store_true")
                  .option(None, "--strict",
                          action="store",
                          type="bool",
                          default="true")).process(args=args)
  if len(command_line.args) != 1:
    command_line.parser.show_help()
    return
  total_timer = time_log("total").start()
  filepath = command_line.args[0]
  if not os.path.isabs(filepath):
    abs_path = libtbx.env.find_in_repositories(relative_path=filepath)
    if abs_path is None:
      abs_path = libtbx.env.find_in_repositories(
        relative_path=filepath, test=os.path.isfile)
    if abs_path is not None: filepath = abs_path
  cif_dics = command_line.options.dictionaries
  if cif_dics is None:
    cif_dics = ["cif_core.dic"]
  cif_dic = validation.smart_load_dictionary(name=cif_dics[0])
  if len(cif_dics) > 1:
    [cif_dic.update(
      validation.smart_load_dictionary(name=d)) for d in cif_dics[1:]]
  show_warnings = command_line.options.show_warnings == True
  show_timings = command_line.options.show_timings == True
  strict = command_line.options.strict
  if os.path.isdir(filepath):
    file_ext = command_line.options.file_ext
    crawl(filepath, file_ext=file_ext,
          cif_dic=cif_dic, show_warnings=show_warnings,
          show_timings=show_timings, strict=strict)
  elif os.path.isfile(filepath):
    cm = cif.reader(file_path=filepath, strict=strict).model()
    cm.validate(cif_dic, show_warnings=show_warnings)
  else:
    try:
      file_object = urllib.request.urlopen(filepath)
    except urllib.error.URLError as e:
      pass
    else:
      cm = cif.reader(file_object=file_object, strict=strict).model()
      cm.validate(cif_dic, show_warnings=show_warnings)
  if show_timings:
    total_timer.stop()
    print(total_timer.report())

def crawl(directory, file_ext, cif_dic, show_warnings, show_timings, strict):
  timer = time_log("parsing")
  validate_timer = time_log("validate")
  for root, dirs, files in os.walk(directory):
    cif_g = glob.glob(os.path.join(root, "*.%s" %file_ext))
    files_to_read = cif_g
    for path in files_to_read:
      timer.start()
      try:
        cm = cif.reader(file_path=path, strict=strict).model()
      except AssertionError:
        continue
      timer.stop()
      s = StringIO()
      validate_timer.start()
      cm.validate(cif_dic, show_warnings=show_warnings, out=s)
      validate_timer.stop()
      if s.getvalue():
        print(path)
        print(s.getvalue())
  if show_timings:
    print(timer.legend)
    print(timer.report())
    print(validate_timer.report())

if __name__ == '__main__':
  run(sys.argv[1:])
