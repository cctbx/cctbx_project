
from __future__ import absolute_import, division, print_function
import libtbx.load_env
from libtbx.utils import multi_out
from optparse import OptionParser
from six.moves import cStringIO as StringIO
import os.path as op
import os
import sys

# XXX these will be skipped (for various reasons)
ignore_modules = set([
  "libtbx.crossmingw",
  "libtbx.start_print_trace",   # produces lots of output
  "libtbx.command_line.epydoc_run",
  "libtbx.command_line.ipython_shell_start",
  "mmtbx.pdb_distances",        # too much executed code
  "gltbx.wx_viewer_leapmotion", # crashes if Leap installed
  # non-CCTBX modules
  "PyQuante.ThomasFermi",       # uses reserved keyword 'with'
  "elbow.example_script",
  "phenix_regression",
  "phenix_dev",
  "phenix.autosol.bayes_5",     # too much executed code
])

def has_init_py(path_name):
  for file_name in os.listdir(path_name):
    if (file_name == "__init__.py"):
      return True
  return False

def split_all(path_name):
  paths = []
  while True :
    leading_path, dir_name = op.split(path_name)
    if (dir_name != ''):
      paths.append(dir_name)
      path_name = leading_path
    else :
      break
  paths.reverse()
  return paths

def run(args):
  verbose = False
  parser = OptionParser()
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
    help="Turn on verbose output")
  parser.add_option("--skip-tests", dest="skip_tests", action="store_true",
    help="Don't import modules beginning with 'tst'")
  options, args = parser.parse_args(args)
  module_list = []
  if (len(args) == 0):
    module_list.extend([ m.name for m in libtbx.env.module_list ])
  else :
    for arg in args :
      assert (arg in libtbx.env.module_dict), arg
      module_list.append(arg)
  has_stdout = []
  stdout_old = sys.stdout
  for module_name in module_list :
    if (module_name in ignore_modules):
      continue
    try :
      module = __import__(module_name)
    except ImportError as e:
      print(e, file=sys.stderr)
      continue
    assert len(module.__path__) == 1
    mod_path = module.__path__[0]
    path_fields = split_all(mod_path)
    n_leading_dirs = len(path_fields) - 1
    for dirname, dirnames, filenames in os.walk(mod_path):
      for file_name in filenames :
        if file_name.endswith(".py") and (file_name != "libtbx_refresh.py"):
          py_mod_name, ext = op.splitext(file_name)
          if (ext != '.py') or ("." in py_mod_name):
            if (options.verbose):
              print("skipping %s" % file_name, file=sys.stderr)
            continue
          py_path = split_all(dirname)[n_leading_dirs:]
          import_name = ".".join(py_path)
          if (not has_init_py(dirname)) or (import_name in ignore_modules):
            continue
          top_level_module = py_path[0]
          if (file_name != "__init__.py"):
            import_name += "." + file_name[:-3]
          if (import_name in ignore_modules):
            continue
          elif ((file_name.startswith("tst_") or file_name.startswith("test_"))
                and options.skip_tests):
            continue
          if (options.verbose):
            print(import_name)
          try :
            sys.stdout = multi_out()
            sys.stdout.register("stdout", stdout_old)
            out = StringIO()
            sys.stdout.register("stringio", out)
            submodule = __import__(import_name)
            if (out.getvalue() != ''):
              has_stdout.append(import_name)
              if (options.verbose):
                print(out.getvalue(), file=sys.stderr)
          except ImportError as e :
            print(e, file=sys.stderr)
          finally :
            sys.stdout = stdout_old
  print("")
  print("*" * 80)
  print("ALL MODULES IMPORTED SUCCESSFULLY")
  print("*" * 80)
  print("")
  if (len(has_stdout) > 0):
    print("Warning: %d modules print to stdout on import" % \
      len(has_stdout), file=sys.stderr)
    for import_name in has_stdout :
      print(import_name, file=sys.stderr)

if (__name__ == "__main__"):
  run(sys.argv[1:])
