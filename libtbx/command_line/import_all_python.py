
from __future__ import division
import libtbx.load_env
from libtbx.utils import multi_out
from cStringIO import StringIO
import os.path as op
import os
import sys

# XXX these will be skipped (for various reasons)
ignore_modules = set([
  "libtbx.crossmingw",
  "libtbx.start_print_trace",   # produces lots of output
  "libtbx.command_line.epydoc_run",
  "libtbx.command_line.ipython_shell_start",
  "libtbx.pythonpath.deep_thought",
  "libtbx.pythonpath.stdlib",
  "libtbx.pythonpath.optik",
  "mmtbx.pdb_distances",        # too much executed code
  "gltbx.wx_viewer_leapmotion", # crashes if Leap installed
  # non-CCTBX modules
  "PyQuante.ThomasFermi",       # uses reserved keyword 'with'
  "elbow.example_script",
])

def has_init_py (path_name) :
  for file_name in os.listdir(path_name) :
    if (file_name == "__init__.py") :
      return True
  return False

def split_all (path_name) :
  paths = []
  while True :
    leading_path, dir_name = op.split(path_name)
    if (dir_name != '') :
      paths.append(dir_name)
      path_name = leading_path
    else :
      break
  paths.reverse()
  return paths

def run (args) :
  verbose = False
  if ("--verbose" in args) :
    verbose = True
    args.remove("--verbose")
  module_list = []
  if (len(args) == 0) :
    module_list.extend([ m.name for m in libtbx.env.module_list ])
  else :
    for arg in args :
      assert (arg in libtbx.env.module_dict), arg
      module_list.append(arg)
  has_stdout = []
  stdout_old = sys.stdout
  for module_name in module_list :
    if (module_name in ignore_modules) :
      continue
    try :
      module = __import__(module_name)
    except ImportError, e :
      print >> sys.stderr, e
      continue
    assert len(module.__path__) == 1
    mod_path = module.__path__[0]
    path_fields = split_all(mod_path)
    n_leading_dirs = len(path_fields) - 1
    for dirname, dirnames, filenames in os.walk(mod_path) :
      for file_name in filenames :
        if file_name.endswith(".py") and (file_name != "libtbx_refresh.py") :
          py_mod_name, ext = op.splitext(file_name)
          if (ext != '.py') or ("." in py_mod_name) :
            if (verbose) :
              print >> sys.stderr, "skipping %s" % file_name
            continue
          py_path = split_all(dirname)[n_leading_dirs:]
          import_name = ".".join(py_path)
          if (file_name != "__init__.py") :
            import_name += "." + file_name[:-3]
          if (import_name in ignore_modules) :
            continue
          if (verbose) :
            print import_name
          try :
            sys.stdout = multi_out()
            sys.stdout.register("stdout", stdout_old)
            out = StringIO()
            sys.stdout.register("stringio", out)
            submodule = __import__(import_name)
            if (out.getvalue() != '') :
              has_stdout.append(import_name)
              if (verbose) :
                print >> sys.stderr, out.getvalue()
          except ImportError, e :
            print >> sys.stderr, e
          finally :
            sys.stdout = stdout_old
  if (len(has_stdout) > 0) :
    print >> sys.stderr, "Warning: %d modules print to stdout on import" % \
      len(has_stdout)
    for import_name in has_stdout :
      print >> sys.stderr, import_name

if (__name__ == "__main__") :
  run(sys.argv[1:])
