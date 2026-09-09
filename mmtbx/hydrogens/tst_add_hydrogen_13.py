from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from mmtbx.hydrogens import reduce_hydrogen

def run():
  test_001()

# ------------------------------------------------------------------------------

def test_001():
  '''
    Occupancy of H atoms placed into alternate conformations.

    Only PART of this GOL is split: C1, C2, O1 and O2 are shared (blank altloc,
    occupancy 1.00) while C3 and O3 exist as conformers A (0.60) and B (0.40).

    reduce2 places H into the alt-conf atom groups only, so a hydrogen whose
    parent is a shared atom - H2 on C2, H11/H12 on C1, HO1 on O1, HO2 on O2 -
    is DUPLICATED, one copy per conformer. Each copy must carry the occupancy
    of the conformer it sits in, so that the copies sum back to the parent's
    occupancy. Giving each copy the parent's occupancy instead makes a single
    hydrogen exist at total occupancy 2.0.

    Seen on 9owa (A1CEW A 1002) and 9udu (A1L81 A 401), where the resulting
    conformer occupancies of 0.85/0.78 summed to 1.63 and validate_ligands
    flagged a correctly modelled ligand.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_001.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model=model)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()

  # occupancy of every atom, keyed by (altloc, name)
  occ = {}
  for atom in model_h_added.get_hierarchy().atoms():
    occ[(atom.parent().altloc.strip(), atom.name.strip())] = atom.occ

  # the conformers themselves are unchanged
  for name in ('C3', 'O3'):
    assert approx_equal(occ[('A', name)], 0.60, eps=1.e-4), (name, occ[('A', name)])
    assert approx_equal(occ[('B', name)], 0.40, eps=1.e-4), (name, occ[('B', name)])

  # H on a split parent: follows that parent. This already worked.
  for name in ('H31', 'H32', 'HO3'):
    assert approx_equal(occ[('A', name)], 0.60, eps=1.e-4), (name, occ[('A', name)])
    assert approx_equal(occ[('B', name)], 0.40, eps=1.e-4), (name, occ[('B', name)])

  # H on a SHARED parent, duplicated into both conformers: each copy must take
  # the occupancy of its own conformer, not the parent's 1.00.
  for name in ('H2', 'H11', 'H12', 'HO1', 'HO2'):
    assert approx_equal(occ[('A', name)], 0.60, eps=1.e-4), (name, occ[('A', name)])
    assert approx_equal(occ[('B', name)], 0.40, eps=1.e-4), (name, occ[('B', name)])

  # and the invariant that motivates it: every duplicated H sums to the
  # occupancy of the parent it belongs to
  for name in ('H2', 'H11', 'H12', 'HO1', 'HO2', 'H31', 'H32', 'HO3'):
    total = occ[('A', name)] + occ[('B', name)]
    assert approx_equal(total, 1.0, eps=1.e-4), (name, total)

  # no H was left at the placement default of zero
  for (altloc, name), value in occ.items():
    if name.startswith('H'):
      assert value > 0.0, (altloc, name, value)

# ------------------------------------------------------------------------------

pdb_str_001 = '''
CRYST1   40.000   30.000   30.000  90.00  90.00  90.00 P 1
SCALE1      0.025000  0.000000  0.000000        0.00000
SCALE2      0.000000  0.033333  0.000000        0.00000
SCALE3      0.000000  0.000000  0.033333        0.00000
HETATM    1  C1  GOL A   1       5.000   5.000   5.000  1.00 20.00           C
HETATM    2  C2  GOL A   1       6.520   5.000   5.000  1.00 20.00           C
HETATM    3  O1  GOL A   1       4.300   3.850   5.400  1.00 20.00           O
HETATM    4  O2  GOL A   1       7.220   4.050   5.700  1.00 20.00           O
HETATM    5  C3 AGOL A   1       7.100   6.400   5.000  0.60 20.00           C
HETATM    6  O3 AGOL A   1       8.500   6.350   4.700  0.60 20.00           O
HETATM    7  C3 BGOL A   1       7.000   6.500   5.200  0.40 20.00           C
HETATM    8  O3 BGOL A   1       8.400   6.500   4.900  0.40 20.00           O
END
'''

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
