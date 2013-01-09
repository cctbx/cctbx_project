
from __future__ import division

def exercise () :
  from mmtbx.regression import tst_fmodel_twin_law
  from libtbx import easy_run
  pdb, mtz = tst_fmodel_twin_law.exercise("tmp_model_vs_data_twinned")
  result = easy_run.fully_buffered(
    "phenix.model_vs_data %s %s" % (pdb, mtz)).raise_if_errors()
  assert ("""    twinned                              : l,-k,h """ in
    result.stdout_lines)

if (__name__ == "__main__") :
  exercise()
  print "OK"
