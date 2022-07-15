from __future__ import absolute_import, division, print_function
from mmtbx.command_line.find_residue_in_pdb import run
from six.moves import cStringIO as StringIO

#
# May be unstable due to RCSB availability.
#

def exercise_1():
  out = StringIO()
  run(['nag'], out=out)
  v = out.getvalue()
  assert v.find('PDB IDs retrieved')>0

if (__name__ == "__main__"):
  exercise_1()
  print('OK')
