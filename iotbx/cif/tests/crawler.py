import glob, os, sys
import urllib2

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
      xs = reader.build_crystal_structure()
    elif build_miller_arrays:
      ma = reader.build_miller_arrays()

def run_once(file_path, build_miller_arrays=False, build_xray_structure=False):
  reader = cif.reader(file_path=file_path, max_errors=10)
  if build_xray_structure:
    xs = reader.build_crystal_structure()
  elif build_miller_arrays:
    ma = reader.build_miller_arrays()

def run(args, out=sys.stdout):
  assert len(args) > 0
  command_line = (option_parser()
                  .option(None, "--file_ext",
                          action="store",
                          default="cif")
                  .option(None, "--build_xray_structure",
                          action="store_true")
                  .option(None, "--build_miller_arrays",
                          action="store_true")).process(args=args[1:])
  filepath = args[0]
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
    run_once(filepath, build_miller_arrays=build_miller_arrays,
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
