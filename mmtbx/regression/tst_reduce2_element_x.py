"""reduce2 removes atoms with element X (unknown element) at loading: the log
and the description file name them; they are not in the output.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from six.moves import cStringIO as StringIO
from libtbx.utils import null_out
from iotbx.cli_parser import CCTBXParser
from mmtbx.programs import reduce2

pdb_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       6.458   5.000   5.000  1.00 10.00           C
ATOM      3  C   ALA A   1       7.009   6.420   5.000  1.00 10.00           C
ATOM      4  O   ALA A   1       6.252   7.390   5.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       6.986   4.235   6.210  1.00 10.00           C
ATOM      6  N   GLY A   2       8.330   6.540   5.000  1.00 10.00           N
ATOM      7  CA  GLY A   2       8.960   7.850   5.000  1.00 10.00           C
ATOM      8  C   GLY A   2      10.470   7.740   5.000  1.00 10.00           C
ATOM      9  O   GLY A   2      11.050   6.650   5.000  1.00 10.00           O
HETATM   10  UNK UNX A 101      15.000  15.000  15.000  1.00 20.00           X
HETATM   11  UNKAUNX B 102      20.000  20.000  20.000  0.50 20.00           X
END
"""

LINE = ('Removed 2 atom(s) with element X: the element is unknown, so H cannot '
  'be placed around them; they are not in the output.')
ATOMS = ['  UNX A 101 UNK', '  UNX B 102 UNK altloc A']

def test_element_x():
  prefix = 'tst_reduce2_element_x'
  fn_in, out, txt = prefix + '_in.pdb', prefix + '_out.pdb', prefix + '.txt'
  with open(fn_in, 'w') as f:
    f.write(pdb_str)
  for n in (out, txt):
    if os.path.exists(n): os.remove(n)
  log = StringIO()
  parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
  parser.parse_args([fn_in, 'approach=add', 'output.filename=' + out,
    'output.description_file_name=' + txt])
  task = reduce2.Program(parser.data_manager, parser.working_phil.extract(),
    master_phil=parser.master_phil, logger=log)
  task.validate()
  task.run()
  for text in (log.getvalue(), open(txt).read()):
    lines = text.splitlines()
    assert LINE in lines, text
    i = lines.index(LINE)
    assert lines[i+1:i+3] == ATOMS, lines[i:i+4]
  h = iotbx.pdb.input(out).construct_hierarchy()
  assert not [a for a in h.atoms() if a.element.strip().upper() == 'X']
  assert not [ag for ag in h.atom_groups() if ag.resname == 'UNX']
  assert sum(1 for a in h.atoms() if a.element.strip() == 'H') > 0

if __name__ == '__main__':
  test_element_x()
  print('OK')
