from __future__ import division
import os
import libtbx.load_env

include = ["cctbx", "iotbx", "mmtbx", "scitbx"]

def find_all_tst_dot_py_files(root_dir):
  result = []
  for root, dirs, files in os.walk(root_dir):
    for f in files:
      if([f.startswith("tst"),f.startswith("test")].count(True)==1 and
          f.endswith(".py")):
        result.append("/".join([root,f]))
  return result

def extract_tests_from_run_tests_dot_py(file_name, full_path):
  fo = open(file_name,"r")
  unique_files = []
  for line in fo.readlines():
    line = line.strip()
    if(line.count("$D") and line.count(".py") and
       (line.count("") or line.count("test")) and not line.count("#")):
      line = line.split()[0]
      for c in ["[","$D",",",'"',"]"]: line=line.replace(c,"")
      if(not line in unique_files): unique_files.append(line)
  fo.close()
  result = []
  for uf in unique_files:
    result.append(full_path+uf)
  for r in result:
    assert os.path.isfile(r), r
  return result

def find_mismatch(from_run_tests, from_dirs, full_path):
  print "  Number of tests in:"
  print "    run_tests.py         :", len(from_run_tests)
  print "    actual in folder-tree:", len(from_dirs)
  if(len(from_run_tests)!=len(from_dirs)):
    print "  Tests listed below are never executed:"
  for f1 in from_dirs:
    if(not f1 in from_run_tests):
      if(1): print "      ", f1

def run():
  root_dir = libtbx.env.find_in_repositories("cctbx_project")
  for subdir in os.listdir(root_dir):
    full_path = "/".join([root_dir,subdir])
    if(subdir in include):
      assert os.path.isdir(full_path)
      print full_path
      run_tests_file = "/".join([full_path, "run_tests.py"])
      assert os.path.isfile(run_tests_file)
      #
      from_run_tests = extract_tests_from_run_tests_dot_py(
        file_name=run_tests_file, full_path=full_path)
      from_dirs = find_all_tst_dot_py_files(root_dir=full_path)
      find_mismatch(from_run_tests, from_dirs, full_path)

if (__name__ == "__main__"):
  run()
