"""
A rotatable H on a shared atom whose bonded partner carries altlocs: TYR
with OH in the blank atom group and CZ under A and B (reduce2 qmready
makes this shape when it borrows OH between the conformers, mmtbx.reduce2
with 1PWM or 1R6J; deposited models have it too). place_hydrogens adds
HH to A and to B at a bogus position, then exclude_H_on_esterified_O sees
the blank OH bonded to two heavy atoms (CZ A and CZ B), takes it for an
ester O and removes both. Expected: one HH per conformer, bonded to OH,
each with the occupancy of its CZ, as HB2/HB3 on a shared CB with CG A/B.
"""
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen

pdb_str = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM     18  N   TYR A   3      23.329   7.636  13.128  1.00 14.70           N
ATOM     19  CA  TYR A   3      24.196   7.963  14.280  1.00 15.18           C
ATOM     20  C   TYR A   3      25.640   8.150  13.866  1.00 15.91           C
ATOM     21  O   TYR A   3      26.078   7.630  12.836  1.00 15.76           O
ATOM     22  CB  TYR A   3      24.098   6.884  15.350  1.00 15.35           C
ATOM     23  CG  TYR A   3      22.702   6.748  15.883  1.00 14.45           C
ATOM     24  CD1 TYR A   3      22.247   7.575  16.901  1.00 14.80           C
ATOM     25  CD2 TYR A   3      21.808   5.840  15.308  1.00 15.68           C
ATOM     26  CE1 TYR A   3      20.941   7.468  17.397  1.00 14.33           C
ATOM     27  CE2 TYR A   3      20.517   5.725  15.777  1.00 13.46           C
ATOM     28  CZ ATYR A   3      20.084   6.548  16.812  0.60 15.09           C
ATOM     29  CZ BTYR A   3      20.084   6.548  16.812  0.40 15.09           C
ATOM     30  OH  TYR A   3      18.803   6.408  17.272  1.00 14.39           O
END
"""

def run():
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info=None, lines=pdb_str),
    log = null_out())
  obj = reduce_hydrogen.place_hydrogens(model = model,
    use_neutron_distances = False)
  obj.run()
  atoms = obj.get_model().get_hierarchy().atoms()
  oh = [a for a in atoms if a.name.strip() == "OH"][0]
  hh = dict((a.parent().altloc, a) for a in atoms if a.name.strip() == "HH")
  assert sorted(hh.keys()) == ["A", "B"], sorted(hh.keys())
  for altloc, occ in (("A", 0.6), ("B", 0.4)):
    a = hh[altloc]
    assert a.distance(oh) < 1.1, (altloc, a.distance(oh))
    assert abs(a.occ - occ) < 1e-6, (altloc, a.occ)

if (__name__ == "__main__"):
  run()
  print("OK")
