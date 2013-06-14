
from __future__ import division
from cStringIO import StringIO

def exercise () :
  from mmtbx.regression import make_fake_anomalous_data
  from mmtbx.command_line import validate_waters
  mtz_file, pdb_file = make_fake_anomalous_data.generate_calcium_inputs(
    file_base="tst_validate_waters")
  out = StringIO()
  waters = validate_waters.run(
    args=[pdb_file, mtz_file, "skip_twin_detection=True"],
    out=out)
  assert len(waters) == 3
  assert (waters[0].b_iso == 2.94)
  print "OK"

if (__name__ == "__main__") :
  exercise()
