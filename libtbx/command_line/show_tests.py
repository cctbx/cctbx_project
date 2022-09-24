from __future__ import absolute_import, division, print_function
import os
import importlib

def match_dir_with_called_in_main(dir_list, f):
  lines_after_main = []
  found_main = False
  with open(f, "r") as fo:
    for l in fo.readlines():
      l = l.strip()
      if(l.count("__main__")): found_main = True
      if(found_main): lines_after_main.append(l)
  matches = 0
  for lid in dir_list:
    for lam in lines_after_main:
      if(lam.startswith(lid)): matches += 1
  if(matches != 1):
    print(lines_after_main)
    print(dir_list)
    assert matches == 1, matches

def get_dirs(realpath):
  modules_found = False
  result = []
  for d in os.path.dirname(realpath).split("/"):
    if(d == "modules"):
      modules_found = True
      continue
    if(modules_found):
      result.append(d)
  result = ".".join(result)
  return result

def run():
  """
  FACTS about this utility:
  0) To execute: run libtbx.show_tests with no args in the current directory.
  1) Show doc strings of test files in the current directory where it is
     executed.
  2) Does not explore sub-directories recursively.
  3) Assume test file name begins with 'tst_' and ends with '.py'.
  4) Tries to assert one test file contains one function that runs the test.
  5) modules/phenix_regression/real_space_refine was used as a model for
     developemnt of this utility.
  """
  for f in os.listdir("."):
    if(f.startswith("tst") and f.endswith(".py")):
      realpath = os.path.realpath(f)
      print( f )
      module = get_dirs(realpath)
      i = importlib.import_module("%s.%s"%(module,f.strip(".py")))
      doc = i.run.__doc__
      print(doc)
      d = dir(i)
      match_dir_with_called_in_main(d, f)

if(__name__ == "__main__"):
  run()
