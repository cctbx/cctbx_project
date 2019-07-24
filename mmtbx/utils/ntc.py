from __future__ import absolute_import, division, print_function

import json
import gzip
import sys

class ntc_validation_results(object):
  def __init__(self, results_fname=None):
    self.js_data = None
    if results_fname is None:
      # Not clear if we ever want to have an object with empty data
      return
    f = gzip.open(results_fname)
    self.js_data = json.load(f)
    f.close()

  def get_number_of_residues(self):
    if self.js_data is not None:
      return self.js_data.get('num_std_residues', None)
    return None

  def get_number_of_steps(self):
    if self.js_data is not None:
      return self.js_data.get('num_steps', None)
    return None

  def print_number_of_steps(self, out=None):
    if out is None:
      out=sys.stdout
    if self.js_data is not None:
      val = self.js_data.get('num_steps', None)
      if val:
        print("Number of steps: %d" % val, file=out)
      else:
        print("Number of steps is not available.", file=out)

  def print_rmsds_for_steps(self, out=None):
    if out is None:
      out=sys.stdout
    for step_name, rmsds in self.js_data['rmsd_dict'].iteritems():
      print(step_name, file=out) # can be cleverer and parse the name out
      curated_data = []
      for step_type, rmsd in rmsds.iteritems():
        curated_data.append((step_type, float(rmsd))) # why they are in string?
      curated_data.sort(key=lambda tup:tup[1])
      for step_type, rmsd in curated_data:
        print("  %s: %.3f" % (step_type, rmsd), file=out)

# Usage example

# from mmtbx.utils.ntc import ntc_validation_results
# ntc_vr = ntc_validation_results("5j02.json.gz")
# ntc_vr.get_number_of_residues()
# ntc_vr.print_rmsds_for_steps()
