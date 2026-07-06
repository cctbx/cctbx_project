from __future__ import absolute_import, division, print_function
import time
from libtbx.utils import null_out
import mmtbx.model
import iotbx.pdb
from libtbx.test_utils import approx_equal
from mmtbx.regression.tst_validate_ligands import find_lr

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL) # Only show critical errors

# ------------------------------------------------------------------------------

def run():
  run_test19()
  run_test20()
  run_test21()
  run_test22()
  run_test23()
  run_test24()
  run_test25()
  run_test26()

# ------------------------------------------------------------------------------

def run_test19():
  print('test19')
  from mmtbx.validation import validate_ligands as vlmod
  p = vlmod.master_params().extract().validate_ligands.alt_conf
  assert approx_equal(p.overlap_dist, 2.0)
  assert approx_equal(p.sym_overlap_dist, 1.5)
  assert approx_equal(p.occ_tol, 0.02)

# ------------------------------------------------------------------------------

_altconf_state_pdb_str = '''
CRYST1   70.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1 C1   GOL A   1       5.000   5.000   5.000  1.00 20.00           C
HETATM    2 C2   GOL A   1       6.500   5.000   5.000  1.00 20.00           C
HETATM    3 C3   GOL A   1       7.100   6.400   5.000  1.00 20.00           C
HETATM    4 O1   GOL A   1       4.300   3.850   5.000  1.00 20.00           O
HETATM    5 O2   GOL A   1       7.200   4.050   5.000  1.00 20.00           O
HETATM   11 C1  AGOL A   2      15.000   5.000   5.000  0.50 20.00           C
HETATM   12 C2  AGOL A   2      16.500   5.000   5.000  0.50 20.00           C
HETATM   13 C3  AGOL A   2      17.100   6.400   5.000  0.50 20.00           C
HETATM   14 O1  AGOL A   2      14.300   3.850   5.000  0.50 20.00           O
HETATM   15 O2  AGOL A   2      17.200   4.050   5.000  0.50 20.00           O
HETATM   21 C1  BGOL A   2      15.200   5.000   5.000  0.50 20.00           C
HETATM   22 C2  BGOL A   2      16.700   5.000   5.000  0.50 20.00           C
HETATM   23 C3  BGOL A   2      17.300   6.400   5.000  0.50 20.00           C
HETATM   24 O1  BGOL A   2      14.500   3.850   5.000  0.50 20.00           O
HETATM   25 O2  BGOL A   2      17.400   4.050   5.000  0.50 20.00           O
HETATM   31 C1  AGOL A   3      25.000   5.000   5.000  0.50 20.00           C
HETATM   32 C2  AGOL A   3      26.500   5.000   5.000  0.50 20.00           C
HETATM   33 C3  AGOL A   3      27.100   6.400   5.000  0.50 20.00           C
HETATM   34 O1  AGOL A   3      24.300   3.850   5.000  0.50 20.00           O
HETATM   35 O2  AGOL A   3      27.200   4.050   5.000  0.50 20.00           O
HETATM   41 C1  BGOL A   3      25.200   5.000   5.000  0.20 20.00           C
HETATM   42 C2  BGOL A   3      26.700   5.000   5.000  0.20 20.00           C
HETATM   43 C3  BGOL A   3      27.300   6.400   5.000  0.20 20.00           C
HETATM   44 O1  BGOL A   3      24.500   3.850   5.000  0.20 20.00           O
HETATM   45 O2  BGOL A   3      27.400   4.050   5.000  0.20 20.00           O
HETATM   51 C1  AGOL A   4      35.000   5.000   5.000  0.50 20.00           C
HETATM   52 C2  AGOL A   4      36.500   5.000   5.000  0.50 20.00           C
HETATM   53 C3  AGOL A   4      37.100   6.400   5.000  0.50 20.00           C
HETATM   54 O1  AGOL A   4      34.300   3.850   5.000  0.50 20.00           O
HETATM   55 O2  AGOL A   4      37.200   4.050   5.000  0.50 20.00           O
HETATM   61 C1  AGOL A   5      45.000   5.000   5.000  0.50 20.00           C
HETATM   62 C2  AGOL A   5      46.500   5.000   5.000  0.50 20.00           C
HETATM   63 C3  AGOL A   5      47.100   6.400   5.000  0.50 20.00           C
HETATM   64 O1  AGOL A   5      44.300   3.850   5.000  0.50 20.00           O
HETATM   65 O2  AGOL A   5      47.200   4.050   5.000  0.50 20.00           O
HETATM   71 C1  BGOL A   6      45.300   5.000   5.000  0.50 20.00           C
HETATM   72 C2  BGOL A   6      46.800   5.000   5.000  0.50 20.00           C
HETATM   73 C3  BGOL A   6      47.400   6.400   5.000  0.50 20.00           C
HETATM   74 O1  BGOL A   6      44.600   3.850   5.000  0.50 20.00           O
HETATM   75 O2  BGOL A   6      47.500   4.050   5.000  0.50 20.00           O
'''

def _altconf_manager(pdb_str):
  from mmtbx.validation import validate_ligands as vlmod
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  params = vlmod.master_params().extract().validate_ligands
  params.ligand_code = []
  vl_manager = vlmod.manager(
    model=model, fmodel=None, map_manager=None, params=params, log=null_out())
  vl_manager.run()
  return vl_manager

def _altconf_lr(vl_manager, resseq, altloc=''):
  sel = 'chain A and resseq %d' % resseq
  if altloc:
    sel += ' and altloc %s' % altloc
  return find_lr(vl_manager, sel)

def run_test20():
  print('test20')
  vl = _altconf_manager(_altconf_state_pdb_str)
  ac1 = _altconf_lr(vl, 1).get_alt_conf()
  assert ac1.state == 'single', ac1.state
  assert ac1.flag == 'ok', ac1.flag

  ac2 = _altconf_lr(vl, 2, 'A').get_alt_conf()
  assert ac2.state == 'alt_conf', ac2.state
  assert ac2.altlocs == ['A', 'B'], ac2.altlocs

  ac4 = _altconf_lr(vl, 4, 'A').get_alt_conf()
  assert ac4.state == 'lone_altloc', ac4.state
  assert ac4.flag == 'inspect', ac4.flag
  assert ac4.partner is None

  ac5 = _altconf_lr(vl, 5, 'A').get_alt_conf()
  assert ac5.state == 'split_residue', ac5.state
  assert ac5.flag == 'inspect', ac5.flag
  assert ac5.partner.resseq == 6, ac5.partner.resseq
  assert ac5.partner.altloc.strip() == 'B', ac5.partner.altloc

  # partner carries chain id and a fully-qualified id string
  from mmtbx.validation.validate_ligands import _partner_id_str
  assert ac5.partner.chain == 'A', ac5.partner.chain
  assert _partner_id_str(ac5.partner) == 'GOL A 6 (altloc B)', _partner_id_str(ac5.partner)

  assert 'GOL A 6 (altloc B)' in ac5.reason, ac5.reason
  from mmtbx.validation.validate_ligands import _alt_conf_short
  assert _alt_conf_short(ac5) == 'split A -> A 6 B', _alt_conf_short(ac5)

  ac6 = _altconf_lr(vl, 6, 'B').get_alt_conf()
  assert ac6.state == 'split_residue', ac6.state
  assert ac6.partner.resseq == 5, ac6.partner.resseq

def run_test21():
  print('test21')
  vl = _altconf_manager(_altconf_state_pdb_str)

  ac2 = _altconf_lr(vl, 2, 'A').get_alt_conf()      # A/B occ 0.5/0.5
  assert ac2.occupancy is not None
  assert approx_equal(ac2.occupancy.occ_sum, 1.0)
  assert ac2.occupancy.sum_ok is True
  assert ac2.flag == 'ok', ac2.flag

  ac3 = _altconf_lr(vl, 3, 'A').get_alt_conf()      # A/B occ 0.5/0.2
  assert approx_equal(ac3.occupancy.occ_sum, 0.7)
  assert ac3.occupancy.sum_ok is False
  assert ac3.flag == 'inspect', ac3.flag

  ac4 = _altconf_lr(vl, 4, 'A').get_alt_conf()      # lone occ 0.5
  assert approx_equal(ac4.occupancy.occ_sum, 0.5)
  assert 'occupancies sum to 0.50' in ac4.reason, ac4.reason

  ac1 = _altconf_lr(vl, 1).get_alt_conf()           # single
  assert ac1.occupancy is None

_altconf_sym_pdb_str = '''
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 2
HETATM    1 C1   GOL A   1       0.300   3.000   0.300  0.50 20.00           C
HETATM    2 C2   GOL A   1       0.300   4.500   0.300  0.50 20.00           C
HETATM    3 C3   GOL A   1       0.300   6.000   0.300  0.50 20.00           C
HETATM    4 O1   GOL A   1       0.300   2.000   0.300  0.50 20.00           O
HETATM    5 O2   GOL A   1       0.300   7.000   0.300  0.50 20.00           O
HETATM   11 C1   GOL A   2       8.000   3.000   8.000  1.00 20.00           C
HETATM   12 C2   GOL A   2       8.000   4.500   8.000  1.00 20.00           C
HETATM   13 C3   GOL A   2       8.000   6.000   8.000  1.00 20.00           C
HETATM   14 O1   GOL A   2       8.000   2.000   8.000  1.00 20.00           O
HETATM   15 O2   GOL A   2       8.000   7.000   8.000  1.00 20.00           O
'''

def run_test22():
  print('test22')
  vl = _altconf_manager(_altconf_sym_pdb_str)

  ac1 = _altconf_lr(vl, 1).get_alt_conf()   # near the 2-fold axis
  assert ac1.symmetry is not None
  assert ac1.symmetry.self_overlap is True
  assert ac1.symmetry.min_dist < 1.5, ac1.symmetry.min_dist
  assert ac1.flag == 'inspect', ac1.flag

  ac2 = _altconf_lr(vl, 2).get_alt_conf()   # far from any axis
  assert ac2.symmetry is None
  assert ac2.flag == 'ok', ac2.flag

_altconf_hetero_pdb_str = '''
CRYST1   50.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  C1 AGOL A   1       5.000   5.000   5.000  0.50 20.00           C
HETATM    2  C2 AGOL A   1       6.500   5.000   5.000  0.50 20.00           C
HETATM    3  O1 AGOL A   1       4.300   3.850   5.000  0.50 20.00           O
HETATM    4  C1 BEDO A   1       5.100   5.100   5.000  0.50 20.00           C
HETATM    5  C2 BEDO A   1       6.600   5.100   5.000  0.50 20.00           C
HETATM    6  O1 BEDO A   1       4.400   3.950   5.000  0.50 20.00           O
HETATM   11  C1 AGOL A 501      25.000   5.000   5.000  0.50 20.00           C
HETATM   12  C2 AGOL A 501      26.500   5.000   5.000  0.50 20.00           C
HETATM   13  O1 AGOL A 501      24.300   3.850   5.000  0.50 20.00           O
HETATM   21  C1 BEDO A 502      25.200   5.100   5.000  0.50 20.00           C
HETATM   22  C2 BEDO A 502      26.700   5.100   5.000  0.50 20.00           C
HETATM   23  O1 BEDO A 502      24.500   3.950   5.000  0.50 20.00           O
ATOM     31  N  ALYS A 700      40.150   6.800   5.000  0.50 20.00           N
ATOM     32  CA ALYS A 700      40.150   7.900   5.000  0.50 20.00           C
HETATM   41  C1 AGOL A 701      40.000   5.000   5.000  0.50 20.00           C
HETATM   42  C2 AGOL A 701      41.500   5.000   5.000  0.50 20.00           C
HETATM   43  O1 AGOL A 701      39.300   3.850   5.000  0.50 20.00           O
'''

def run_test23():
  print('test23')
  vl = _altconf_manager(_altconf_hetero_pdb_str)

  # same residue group, altloc A = GOL, altloc B = EDO
  acA = _altconf_lr(vl, 1, 'A').get_alt_conf()
  assert acA.state == 'alt_conf', acA.state
  assert acA.hetero is True, acA.hetero
  assert acA.resnames == ['EDO', 'GOL'], acA.resnames
  assert acA.flag == 'inspect', acA.flag
  assert 'different chemical entities' in acA.reason, acA.reason

  # split across resseq, different compounds
  acS = _altconf_lr(vl, 501, 'A').get_alt_conf()
  assert acS.state == 'split_residue', acS.state
  assert acS.hetero is True
  assert acS.partner.resname == 'EDO', acS.partner.resname

  # class-filter control: lone GOL with a LYS altloc sidechain within 2 A
  acL = _altconf_lr(vl, 701, 'A').get_alt_conf()
  assert acL.state == 'lone_altloc', acL.state
  assert acL.hetero is False, acL.hetero

def run_test24():
  print('test24')
  from six.moves import cStringIO as StringIO
  vl = _altconf_manager(_altconf_state_pdb_str)

  snap5 = _altconf_lr(vl, 5, 'A').as_picklable_snapshot()
  assert snap5.alt_conf is not None
  assert snap5.alt_conf.state == 'split_residue'
  assert snap5.alt_conf.flag == 'inspect'
  assert snap5.alt_conf.partner.resseq == 6
  assert snap5.alt_conf.partner.chain == 'A', snap5.alt_conf.partner.chain
  assert snap5.alt_conf.hetero is False
  assert snap5.alt_conf.resnames == ['GOL']
  assert snap5.alt_conf.occupancy.sum_ok is True
  assert snap5.alt_conf.symmetry is None

  snap1 = _altconf_lr(vl, 1).as_picklable_snapshot()
  assert snap1.alt_conf.state == 'single'
  assert snap1.alt_conf.partner is None
  assert snap1.alt_conf.occupancy is None

  sio = StringIO()
  vl.show_table(out=sio)
  text = sio.getvalue()
  assert 'alt' in text and 'conf' in text, text
  assert 'split' in text, text

# ------------------------------------------------------------------------------

# A bonded ligand at a general position in a screw-axis space group. Its own
# intramolecular bonds must NOT be reported as a symmetry mate overlapping it
# (regression: identity-operation pairs were slipping through when filtering on
# pair.j_sym instead of the actual rt_mx).
_altconf_symfp_pdb_str = '''
CRYST1   20.000   24.000   28.000  90.00  90.00  90.00 P 21 21 21
HETATM    1 C1   GOL A   1      10.000  12.000  14.000  1.00 20.00           C
HETATM    2 C2   GOL A   1      11.500  12.000  14.000  1.00 20.00           C
HETATM    3 O1   GOL A   1      10.700  13.300  14.000  1.00 20.00           O
'''

def run_test25():
  print('test25')
  vl = _altconf_manager(_altconf_symfp_pdb_str)
  ac = _altconf_lr(vl, 1).get_alt_conf()
  # No crystallographic symmetry mate actually overlaps this general-position
  # ligand; only its own bonds are within the cutoff, and those are identity
  # operations that must be filtered out.
  assert ac.symmetry is None, ac.symmetry
  assert ac.flag == 'ok', ac.flag

# ------------------------------------------------------------------------------

_altconf_order_pdb_str = '''
CRYST1   70.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1 C1  AGOL A   1       5.000   5.000   5.000  0.50 20.00           C
HETATM    2 C2  AGOL A   1       6.500   5.000   5.000  0.50 20.00           C
HETATM    3 C3  AGOL A   1       7.100   6.400   5.000  0.50 20.00           C
HETATM    4 O1  AGOL A   1       4.300   3.850   5.000  0.50 20.00           O
HETATM    5 O2  AGOL A   1       7.200   4.050   5.000  0.50 20.00           O
HETATM   11 C1   GOL A   2      30.000   5.000   5.000  1.00 20.00           C
HETATM   12 C2   GOL A   2      31.500   5.000   5.000  1.00 20.00           C
HETATM   13 C3   GOL A   2      32.100   6.400   5.000  1.00 20.00           C
HETATM   14 O1   GOL A   2      29.300   3.850   5.000  1.00 20.00           O
HETATM   15 O2   GOL A   2      32.200   4.050   5.000  1.00 20.00           O
HETATM   21 C1  BGOL A   3       5.200   5.000   5.000  0.50 20.00           C
HETATM   22 C2  BGOL A   3       6.700   5.000   5.000  0.50 20.00           C
HETATM   23 C3  BGOL A   3       7.300   6.400   5.000  0.50 20.00           C
HETATM   24 O1  BGOL A   3       4.500   3.850   5.000  0.50 20.00           O
HETATM   25 O2  BGOL A   3       7.400   4.050   5.000  0.50 20.00           O
'''

def run_test26():
  print('test26')
  vl = _altconf_manager(_altconf_order_pdb_str)
  order = vl._ordered_for_display()
  # id_str glues the altloc letter onto the resname for altloc-bearing
  # ligands (e.g. 'AGOL A   1', see tst_validate_ligands.py::run_test02),
  # so identify rows by (resname, chain, resseq, altloc) instead.
  def _key(lr):
    rg = lr._atoms_ligand[0].parent().parent()
    return (lr.resname, rg.parent().id.strip(), rg.resseq_as_int(),
             lr.altloc.strip())
  keys = [_key(lr) for lr in order]
  # every ligand appears exactly once
  assert len(order) == len(vl), (len(order), len(vl))
  # partner (resseq 3, altloc B) immediately follows resseq 1 altloc A,
  # ahead of the unrelated single GOL resseq 2
  assert keys[0] == ('GOL', 'A', 1, 'A'), keys
  assert keys[1] == ('GOL', 'A', 3, 'B'), keys
  assert keys[2] == ('GOL', 'A', 2, ''), keys

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
