from __future__ import division, print_function
import os
from libtbx import easy_run

def run():
  '''
  Loop through scripts in cctbx_website/examples and check if they run
  correctly
  '''
  regression_dir = os.path.dirname(os.path.abspath(__file__))
  root_dir = os.path.dirname(regression_dir)
  examples_dir = os.path.join(root_dir, 'examples')
  tmp_dir = os.path.join(regression_dir, 'tmp_files')

  # make temporary directy "tmp_dir"
  if (not os.path.isdir(tmp_dir)):
    os.makedirs(tmp_dir)
  os.chdir(tmp_dir)

  results = list()
  # loop through all .py files in "examples"
  for script in os.listdir(examples_dir):
    cmd = 'libtbx.python ' + os.path.join(examples_dir, script)
    r = easy_run.fully_buffered(cmd)
    results.append([script, r.return_code, r.stdout_lines, r.stderr_lines])
    # remove files once done
    # (the scripts sometimes on the same input files, so need to delete each time)
    for f in os.listdir(tmp_dir):
      os.remove(os.path.join(tmp_dir, f))
  #
  # go back up to "regression" and delete tmp directory
  os.chdir(root_dir)
  os.rmdir(tmp_dir)

  # print info if fail or success; print stderr if failed
  for l in results:
    if (l[1] == 0): re = 'ran successfully'
    else: re = 'failed'
    print('%s %s  ' % (l[0], re))
    if l[1] != 0:
      for line in l[3]:
        print('\t', line)

if __name__ == '__main__':
  run()
