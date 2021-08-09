from __future__ import division
import os
import tempfile
import libtbx
import libtbx.load_env
from libtbx.easy_run import fully_buffered


def run():
  lrl_path = libtbx.env.find_in_repositories('cctbx_project/cctbx/uctbx/lrl')
  ml_path = os.path.join(lrl_path, 'match_lattices.py')

  inputs = [
      "r 57.98 57.98 57.98 92.02 92.02 92.02",
      "r 57.98 57.98 57.98 92.02 92.02 92.02",
      "h 80.36 80.36 99.44 90 90 120",
      "c 80.95 80.57 57.1 90 90.35 90",
      "h 80.36 80.36 99.44 90 90 120",
      "r 57.1 57.1 57.1 89.75 89.75 89.75",
      ]
  inp_file = tempfile.NamedTemporaryFile(mode='w')
  for line in inputs: inp_file.write(line + '\n')
  inp_file.flush()

  outputs = [
      "0    57.98000  57.98000  57.98000  92.02000  92.02000  92.02000",
      "1    57.98000  57.98000  57.98000  92.02000  92.02000  92.02000",
      "2    57.01998  57.01998  57.01998  89.60502  90.39498  90.39498",
      "3    57.10000  57.10610  57.10610  90.26959  89.75193  90.24807",
      "4    57.01998  57.01998  57.01998  89.60502  90.39498  90.39498",
      "5    57.10000  57.10000  57.10000  90.25000  89.75000  90.25000",
  ]

  command = "libtbx.python %s %s" %(ml_path, inp_file.name)
  result = fully_buffered(command).raise_if_errors().stdout_lines[-6:]

  for l, ref in zip(result, outputs):
    assert l.strip()==ref

  inp_file.close()

if __name__ == "__main__":
  run()
