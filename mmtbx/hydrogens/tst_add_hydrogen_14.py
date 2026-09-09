from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen

def run():
  test_001()
  test_002()

# ------------------------------------------------------------------------------

pdb_str_001 = """
CRYST1   64.290   64.290   32.490  90.00  90.00 120.00 H 3           9
ATOM     57  N   TYR A   4      15.190  27.150  14.710  1.00  7.30           N
ATOM     58  CA  TYR A   4      14.520  26.160  15.530  1.00  9.70           C
ATOM     59  C   TYR A   4      14.080  26.840  16.820  1.00 10.30           C
ATOM     60  O   TYR A   4      14.850  27.500  17.430  1.00 12.10           O
ATOM     61  CB  TYR A   4      15.340  24.840  15.890  1.00  9.80           C
ATOM     62  CG  TYR A   4      15.430  24.000  14.680  1.00  9.60           C
ATOM     63  CD1 TYR A   4      14.750  22.680  14.730  1.00 10.90           C
ATOM     64  CD2 TYR A   4      16.080  24.240  13.530  1.00 11.30           C
ATOM     65  CE1 TYR A   4      14.800  21.820  13.580  1.00 13.80           C
ATOM     66  CE2 TYR A   4      16.060  23.380  12.290  1.00 13.20           C
ATOM     67  CZ  TYR A   4      15.410  22.150  12.520  1.00 12.60           C
ATOM     68  OH  TYR A   4      15.490  21.440  11.230  1.00 15.40           O
ATOM     77  N   THR A   5      12.890  26.450  17.390  1.00  9.50           N
ATOM     78  CA  THR A   5      12.550  27.000  18.750  1.00  9.60           C
ATOM     79  C   THR A   5      12.370  25.920  19.750  1.00  7.50           C
ATOM     80  O   THR A   5      11.900  24.740  19.380  1.00 10.00           O
ATOM     81  CB  THR A   5      11.140  27.660  18.450  1.00 18.30           C
ATOM     82  OG1 THR A   5      10.930  28.360  19.760  1.00 20.90           O
ATOM     83  CG2 THR A   5      10.520  27.760  17.470  1.00 24.90           C
ATOM     86  N   CYS A   6      12.800  26.250  20.890  1.00  8.10           N
ATOM     87  CA  CYS A   6      12.630  25.410  22.060  1.00  8.30           C
ATOM     88  C   CYS A   6      11.040  25.550  22.390  1.00 10.30           C
ATOM     89  O   CYS A   6      10.590  26.610  22.770  1.00 10.10           O
ATOM     90  CB  CYS A   6      13.450  25.910  23.210  1.00  7.70           C
ATOM     91  SG  CYS A   6      13.130  24.890  24.650  1.00 12.50           S
END
"""

def get_model_with_h(pdb_str):
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None),
    log         = null_out())
  obj = reduce_hydrogen.place_hydrogens(model = model)
  obj.run()
  return obj

# ------------------------------------------------------------------------------

def test_001():
  '''
    A bond the restraints already know about is not a missing link.

    workarounds_00345 removes a pair of H closer than 1.0 A whose parent heavy
    atoms are within 1.6 A, on the premise that the parents carry an unperceived
    bond and so should never have been hydrogenated. A real bond satisfies that
    distance test too. Residues are hydrogenated with methyls at an arbitrary
    torsion and rotated to staggered later, so an as-placed methyl can sit
    eclipsed against its own parent's H: here THR 5 CB-CG2 is an ordinary 1.52 A
    bond and HG22 lands 0.84 A from HB. Both were deleted, before the rotation
    that separates them.

    From 4RXN, where THR 5 lost HB and HG22 while THR 7 and THR 28 kept theirs
    only because their starting torsion happened to land further away.
  '''
  h_atoms = {}
  for ag in get_model_with_h(pdb_str_001).get_model().get_hierarchy().atom_groups():
    h_atoms[ag.parent().resseq.strip()] = set(
      a.name.strip() for a in ag.atoms() if a.element.strip() in ('H', 'D'))
  assert 'HB' in h_atoms['5'], h_atoms['5']
  assert 'HG22' in h_atoms['5'], h_atoms['5']
  assert h_atoms['5'] == set(
    ['H', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23']), h_atoms['5']

# ------------------------------------------------------------------------------

def test_002():
  '''
    The reported H count is the count in the model that is returned.

    n_H_final was taken before workarounds_00345 ran, so any H the workaround
    removed was still counted: 4RXN reported 372 added and wrote 370. The number
    is only printed, never asserted on, so a silent deletion stayed invisible
    from both ends.
  '''
  obj = get_model_with_h(pdb_str_001)
  n_reported = obj.get_counts().number_h_final
  n_in_model = obj.get_model().get_hd_selection().count(True)
  assert n_reported == n_in_model, (n_reported, n_in_model)

# ------------------------------------------------------------------------------

if __name__ == '__main__':
  run()
  print("OK")
