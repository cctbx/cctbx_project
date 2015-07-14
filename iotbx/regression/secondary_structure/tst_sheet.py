from __future__ import division
import sys
from iotbx.pdb import secondary_structure as ss
import StringIO
from libtbx.utils import Sorry

ann_1 = """\
SHEET    1   B 4 ARG A  90  PHE A  95  0
SHEET    2   B 4 GLY A 179  ARG A 184 -1  N  GLY A 179   O  PHE A 95
SHEET    3   B 4 ARG A 128  THR A 132 -1  N  ALA A 130   O  ARG A 184
SHEET    4   B 4 VAL A 148  PRO A 152  0  O  ILE A 151   N  VAL A 129
"""

ann_2 = """\
SHEET    1   B 4 ARG A  90  PHE A  95  0
SHEET    2   B 4 GLY A 179  ARG A 184 -1  N  GLY A 179   O  PHE A 95
SHEET    3   B 4 ARG A 128  THR A 132 -1  N  ALA A 130   O  ARG A 184
SHEET    4   B 4 VAL A 148  PRO A 152 -1  O  ILE A 151   N  VAL A 129
"""

def exercise_01():
  restr_groups = StringIO.StringIO()
  try:
    annot = ss.annotation.from_records(records=ann_1.split('\n'))
  except Sorry, e:
    m = str(e)
    assert m.find("should be 1 or -1 for non-first strand") > 0

  annot = ss.annotation.from_records(records=ann_2.split('\n'))
  assert len(annot.helices) == 0
  assert len(annot.sheets) == 1
  sheet = annot.sheets[0]
  assert sheet.n_strands == 4
  assert len(sheet.strands) == 4
  assert len(sheet.registrations) == 4
  assert sheet.registrations[0] is None

def exercise(args):
  exercise_01()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
