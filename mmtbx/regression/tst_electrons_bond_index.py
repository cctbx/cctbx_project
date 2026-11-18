import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.ligands import electrons

# One Cys whose SG bonds its own symmetry mate across the 2-fold: the bond
# proxies have six simple bonds and one symmetry bond (5, 5).
pdb_str = '''
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1 2 1
ATOM      1  N   CYS A   1       4.373   7.170   1.950  1.00 10.00           N
ATOM      2  CA  CYS A   1       3.540   6.237   1.247  1.00 10.00           C
ATOM      3  C   CYS A   1       4.199   4.875   1.150  1.00 10.00           C
ATOM      4  O   CYS A   1       5.420   4.750   1.150  1.00 10.00           O
ATOM      5  CB  CYS A   1       2.144   6.170   1.877  1.00 10.00           C
ATOM      6  SG  CYS A   1       1.020   5.000   0.000  1.00 10.00           S
ATOM      7  OXT CYS A   1       3.450   3.870   1.100  1.00 10.00           O
'''

def bonds_by_scan(ed, i_seq):
  """Return the bonds of one atom, found by scanning every bond."""
  return [i_seqs for i_seqs in ed.generate_bond_i_seqs() if i_seq in i_seqs]

def exercise():
  """get_bonds_containing_i_seq() equals a scan of every bond, in order,
  for every atom and for a bond key; the returned list is a copy."""
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(lines=pdb_str, source_info=None),
    log=null_out())
  model.process(make_restraints=True)
  ed = electrons.electron_distribution(model.get_hierarchy(),
    model.get_restraints_manager().geometry, log=null_out())
  assert ed.simple.size() == 6, ed.simple.size()
  assert ed.asu.size() == 1, ed.asu.size()
  for i_seq in range(model.get_number_of_atoms()):
    got, expected = ed.get_bonds_containing_i_seq(i_seq), bonds_by_scan(ed, i_seq)
    assert got == expected, (i_seq, got, expected)
  got = ed.get_bonds_containing_i_seq(5)
  assert got == [(4, 5), (5, 5)], got
  # the distribution's keys include bonds, which contain no atom
  got = ed.get_bonds_containing_i_seq((4, 5))
  assert got == [], got
  ed.get_bonds_containing_i_seq(1).append((0, 6))
  got, expected = ed.get_bonds_containing_i_seq(1), bonds_by_scan(ed, 1)
  assert got == expected, (got, expected)

if __name__ == '__main__':
  exercise()
  print('OK')
