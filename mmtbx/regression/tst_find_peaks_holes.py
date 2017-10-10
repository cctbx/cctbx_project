
from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from libtbx.utils import null_out
from libtbx import easy_pickle
from io import StringIO
import os

def exercise () :
  from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
  from mmtbx.command_line import find_peaks_holes
  mtz_file, pdb_file = generate_calcium_inputs(
    file_base = "tst_find_peaks_holes", anonymize = True)
  out = StringIO()
  peaks_holes = find_peaks_holes.run(
    args=[pdb_file, mtz_file],
    out=out)
  peaks_holes.save_pdb_file(file_name="tst_fph_peaks.pdb", log=null_out())
  p = easy_pickle.dumps(peaks_holes)
  s = peaks_holes.get_summary()
  sp = easy_pickle.dumps(s)
  out2 = StringIO()
  s.show(out=out2)
  lines = out2.getvalue().splitlines()
  assert ("""  anomalous H2O (anomalous > 3):      1""" in lines)
  assert ("""  anomalous non-water atoms:          0""" in lines)
  assert ("""  mFo-DFc >  9:                       0""" in lines)
  peaks_holes = find_peaks_holes.run(
    args=[pdb_file, mtz_file, "filter_peaks_by_2fofc=1.0"],
    out=null_out())
  out3 = StringIO()
  peaks_holes.get_summary().show(out=out3)
  lines = out3.getvalue().splitlines()
  assert ("""  anomalous > 3:                      0""" in lines)
  out3 = StringIO()
  peaks_holes = find_peaks_holes.run(
    args=[pdb_file, mtz_file, "include_peaks_near_model=True",],
    out=out3)
  lines = out3.getvalue().splitlines()
  assert ("""  mFo-DFc >  9:                       1""" in lines)
  os.remove(mtz_file)
  os.remove(pdb_file)

if (__name__ == "__main__") :
  exercise()
  print("OK")
