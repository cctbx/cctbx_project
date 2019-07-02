from __future__ import absolute_import, division, print_function
import sys
from iotbx.pdb import secondary_structure as ss
from libtbx.utils import Sorry, format_cpu_times
from libtbx import test_utils

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

ann_3 = """\
SHEET    1   3 3 MET A 157  VAL A 161  0
SHEET    2   3 3 ILE A 114  ALA A 122 -1
SHEET    3   3 3 THR A 193  LEU A 201  1
"""

def exercise_01():
  try:
    annot = ss.annotation.from_records(records=ann_1.split('\n'))
  except Sorry as e:
    m = str(e)
    assert m.find("should be 1 or -1 for non-first strand") > 0

def exercise_02():
  annot = ss.annotation.from_records(records=ann_2.split('\n'))
  assert len(annot.helices) == 0
  assert len(annot.sheets) == 1
  sheet = annot.sheets[0]
  assert sheet.n_strands == 4
  assert len(sheet.strands) == 4
  assert len(sheet.registrations) == 4
  assert sheet.registrations[0] is None

def exercise_03():
  answer_phil_str = """\

protein.sheet {
  sheet_id = "  3"
  first_strand = chain 'A' and resid  157  through  161
  strand {
    selection = chain 'A' and resid  114  through  122
    sense = antiparallel
    bond_start_current = None
    bond_start_previous = None
  }
  strand {
    selection = chain 'A' and resid  193  through  201
    sense = parallel
    bond_start_current = None
    bond_start_previous = None
  }
}"""
  annot = ss.annotation.from_records(records=ann_3.split('\n'))
  assert len(annot.helices) == 0
  assert len(annot.sheets) == 1
  ann_sheet = annot.sheets[0]
  assert ann_sheet.n_strands == 3
  assert ann_sheet.registrations == [None]*3
  only_sheet = annot.sheets[0]
  phil_string = annot.as_restraint_groups()
  assert not test_utils.show_diff(phil_string, answer_phil_str,
      strip_trailing_whitespace=True)

def exercise(args):
  exercise_01()
  exercise_02()
  exercise_03()
  print("OK")
  print(format_cpu_times())

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
