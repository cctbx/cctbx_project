from __future__ import division, print_function
import os
import sys
import shutil
import libtbx.load_env
from libtbx import easy_run

def exercise(script, tmp_path, use_pdb_file=False):
  regression_dir = os.path.dirname(os.path.abspath(__file__))
  root_dir = os.path.dirname(regression_dir)
  examples_dir = os.path.join(root_dir, 'examples')
  tmp_dir = os.path.join(regression_dir, tmp_path)

  # make temporary directy "tmp_dir"
  if (not os.path.isdir(tmp_dir)):
    os.makedirs(tmp_dir)
  os.chdir(tmp_dir)

  results = list()
  skipped = False

  if script in ["script_ideal_ss.py"] and not libtbx.env.has_module('phenix'):
    print("phenix not available, skipping test")
    skipped = True
  if script in ["script_compare_ss.py"] and not libtbx.env.has_module('ksdssp'):
    if (not libtbx.env.has_module("ksdssp")):
      print("ksdssp not available, skipping test)")
      skipped = True
  if script in ["doc_programming_tips_3.py"]:
    if sys.platform == 'win32':
      print("Not structured for Windows multiprocessing, skipping test")
      skipped = True
  if not skipped:
    # Some scripts use a PDB file, use one from phenix_regression if available
    if use_pdb_file:
      if script in ['doc_models_hierarchy.py']: #
        path = "phenix_regression/mmtbx/ions/3e0f.pdb"
      else:
        path = "phenix_regression/pdb/1ywf.pdb"
      pdb_file = libtbx.env.find_in_repositories(
        relative_path=path,
        test=os.path.isfile)
      if (pdb_file is None):
        print("phenix_regression not available, skipping test")
        return
      cmd = 'libtbx.python ' + os.path.join(examples_dir, script) + ' ' + pdb_file
    else:
      cmd = 'libtbx.python ' + os.path.join(examples_dir, script)
    # run script from html file
    r = easy_run.fully_buffered(cmd)
    results = [script, r.return_code, r.stdout_lines, r.stderr_lines]

    # go back up to "regression" and delete tmp directory
    os.chdir(root_dir)
    shutil.rmtree(tmp_dir)
    #for f in os.listdir(tmp_dir):
    #  os.remove(os.path.join(tmp_dir, f))

  # parse results to see if it failed
  return_code = 0
  if results:
    if (results[1] == 0): re = 'ran successfully'
    else: re = 'failed'
    print('%s %s  ' % (results[0], re))
    if results[1] != 0:
      return_code = 1
      for line in results[3]:
        print('\t', line, file=sys.stderr)
  if skipped:
    print('%s skipped' % script)

  return return_code
