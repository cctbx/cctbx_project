import glob, os, sys
import urllib2
from cStringIO import StringIO

from libtbx.utils import time_log
import libtbx.load_env
from iotbx.option_parser import option_parser
from iotbx import cif

class crawl(object):
  def __init__(self, directory, file_ext,
               build_miller_arrays=False,
               build_xray_structure=False):
    timer = time_log("parsing")
    error_count = 0
    self.parsing_error_count = 0
    for root, dirs, files in os.walk(directory):
      cif_g = glob.glob(os.path.join(root, "*.%s" %file_ext))
      files_to_read = cif_g
      for path in files_to_read:
        timer.start()
        try:
          reader = self.run_once(path, build_miller_arrays=build_miller_arrays,
                            build_xray_structure=build_xray_structure)

        except Exception, e:
          print "error reading %s" %path
          print e
          error_count += 1
        timer.stop()
    print
    print "%i files read (%i with building errors and %i with parsing errors)" %(
      timer.n, error_count, self.parsing_error_count)
    print timer.legend
    print timer.report()
    sys.stdout.flush()

  def run_once(self, file_path, build_miller_arrays=False,
               build_xray_structure=False):
    reader = cif.reader(file_path=file_path, max_errors=10)
    if reader.error_count(): self.parsing_error_count += 1
    if build_xray_structure:
      xs = reader.build_crystal_structures()
    elif build_miller_arrays:
      ma = reader.build_miller_arrays()

def run_once(file_path, dic=None, build_miller_arrays=False, build_xray_structure=False):
  timer_p = time_log("parse")
  timer_w = time_log("write")
  timer_v = time_log("validate")
  from iotbx.cif import validation
  if dic is None: dic = "cif_core.dic"
  name = dic_path = None
  if os.path.isfile(dic): file_path = dic
  else: name = dic
  dic = validation.smart_load_dictionary(name=name, file_path=dic_path)
  #dic = validation.smart_load_dictionary("cif_core.dic")
  #dic = validation.smart_load_dictionary("cif_mm.dic")
  #dic.update(validation.smart_load_dictionary("mmcif_pdbx.dic"), mode="overlay")
  #dic = validation.smart_load_dictionary("mmcif_pdbx_v40.dic")
  #dic = validation.smart_load_dictionary("mmcif_pdbx.dic")
  #dic = validation.smart_load_dictionary("ddl_core.dic")
  #dic = validation.smart_load_dictionary("mmcif_ddl_2.1.6.dic")
  for i in range(10):
    timer_p.start()
    reader = cif.reader(file_path=file_path, max_errors=10)
    cif_model = reader.model()
    if build_xray_structure:
      xs = reader.build_crystal_structures()
    elif build_miller_arrays:
      ma = reader.build_miller_arrays()
    timer_p.stop()
    timer_w.start()
    print >> StringIO(), cif_model
    timer_w.stop()
    timer_v.start()
    s = StringIO()
    cif_model.validate(dic, out=s)
    timer_v.stop()

  print s.getvalue()
  print
  print timer_p.legend
  print timer_p.report()
  print timer_w.report()
  print timer_v.report()

def run(args, out=sys.stdout):
  command_line = (option_parser()
                  .option(None, "--file_ext",
                          action="store",
                          default="cif")
                  .option(None, "--dic",
                          action="store",
                          default="cif_core.dic")
                  .option(None, "--build_xray_structure",
                          action="store_true")
                  .option(None, "--build_miller_arrays",
                          action="store_true")).process(args=args)
  filepath = command_line.args[0]
  if not os.path.isabs(filepath):
    abs_path = libtbx.env.find_in_repositories(relative_path=filepath)
    if abs_path is not None: filepath = abs_path
  file_ext = command_line.options.file_ext
  build_miller_arrays = command_line.options.build_miller_arrays == True
  build_xray_structure = command_line.options.build_xray_structure == True

  if os.path.isdir(filepath):
    crawl(filepath, file_ext=file_ext,
          build_miller_arrays=build_miller_arrays,
          build_xray_structure=build_xray_structure)
  elif os.path.isfile(filepath):
    run_once(filepath, dic=command_line.options.dic,
             build_miller_arrays=build_miller_arrays,
             build_xray_structure=build_xray_structure)
  else:
    try:
      file_object = urllib2.urlopen(filepath)
    except urllib2.URLError, e:
      pass
    else:
      cm = reader(file_object=file_object).model()

if __name__ == '__main__':
  run(sys.argv[1:])
  print "OK"
