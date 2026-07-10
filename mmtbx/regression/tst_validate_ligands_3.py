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
  run_test27()

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
HETATM    1  C1  GOL A   1       5.578   9.079   8.959  1.00 20.00           C
HETATM    2  C2  GOL A   1       5.404  10.193   9.989  1.00 20.00           C
HETATM    3  C3  GOL A   1       4.003  10.183  10.608  1.00 20.00           C
HETATM    4  O1  GOL A   1       5.482   7.794   9.563  1.00 20.00           O
HETATM    5  O2  GOL A   1       5.628  11.463   9.370  1.00 20.00           O
HETATM    6  O3  GOL A   1       3.905  11.289  11.512  1.00 20.00           O
HETATM    7  C1 AGOL A   2      15.578   9.079   8.959  0.50 20.00           C
HETATM    8  C2 AGOL A   2      15.404  10.193   9.989  0.50 20.00           C
HETATM    9  C3 AGOL A   2      14.003  10.183  10.608  0.50 20.00           C
HETATM   10  O1 AGOL A   2      15.482   7.794   9.563  0.50 20.00           O
HETATM   11  O2 AGOL A   2      15.628  11.463   9.370  0.50 20.00           O
HETATM   12  O3 AGOL A   2      13.905  11.289  11.512  0.50 20.00           O
HETATM   13  C1 BGOL A   2      16.035   9.666   8.959  0.50 20.00           C
HETATM   14  C2 BGOL A   2      15.185  10.407   9.989  0.50 20.00           C
HETATM   15  C3 BGOL A   2      14.119   9.499  10.608  0.50 20.00           C
HETATM   16  O1 BGOL A   2      16.788   8.620   9.563  0.50 20.00           O
HETATM   17  O2 BGOL A   2      14.540  11.525   9.370  0.50 20.00           O
HETATM   18  O3 BGOL A   2      13.333  10.283  11.512  0.50 20.00           O
HETATM   19  C1 AGOL A   3      25.578   9.079   8.959  0.50 20.00           C
HETATM   20  C2 AGOL A   3      25.404  10.193   9.989  0.50 20.00           C
HETATM   21  C3 AGOL A   3      24.003  10.183  10.608  0.50 20.00           C
HETATM   22  O1 AGOL A   3      25.482   7.794   9.563  0.50 20.00           O
HETATM   23  O2 AGOL A   3      25.628  11.463   9.370  0.50 20.00           O
HETATM   24  O3 AGOL A   3      23.905  11.289  11.512  0.50 20.00           O
HETATM   25  C1 BGOL A   3      26.035   9.666   8.959  0.20 20.00           C
HETATM   26  C2 BGOL A   3      25.185  10.407   9.989  0.20 20.00           C
HETATM   27  C3 BGOL A   3      24.119   9.499  10.608  0.20 20.00           C
HETATM   28  O1 BGOL A   3      26.788   8.620   9.563  0.20 20.00           O
HETATM   29  O2 BGOL A   3      24.540  11.525   9.370  0.20 20.00           O
HETATM   30  O3 BGOL A   3      23.333  10.283  11.512  0.20 20.00           O
HETATM   31  C1 AGOL A   4      35.578   9.079   8.959  0.50 20.00           C
HETATM   32  C2 AGOL A   4      35.404  10.193   9.989  0.50 20.00           C
HETATM   33  C3 AGOL A   4      34.003  10.183  10.608  0.50 20.00           C
HETATM   34  O1 AGOL A   4      35.482   7.794   9.563  0.50 20.00           O
HETATM   35  O2 AGOL A   4      35.628  11.463   9.370  0.50 20.00           O
HETATM   36  O3 AGOL A   4      33.905  11.289  11.512  0.50 20.00           O
HETATM   37  C1 AGOL A   5      45.578   9.079   8.959  0.50 20.00           C
HETATM   38  C2 AGOL A   5      45.404  10.193   9.989  0.50 20.00           C
HETATM   39  C3 AGOL A   5      44.003  10.183  10.608  0.50 20.00           C
HETATM   40  O1 AGOL A   5      45.482   7.794   9.563  0.50 20.00           O
HETATM   41  O2 AGOL A   5      45.628  11.463   9.370  0.50 20.00           O
HETATM   42  O3 AGOL A   5      43.905  11.289  11.512  0.50 20.00           O
HETATM   43  C1 BGOL A   6      46.035   9.666   8.959  0.50 20.00           C
HETATM   44  C2 BGOL A   6      45.185  10.407   9.989  0.50 20.00           C
HETATM   45  C3 BGOL A   6      44.119   9.499  10.608  0.50 20.00           C
HETATM   46  O1 BGOL A   6      46.788   8.620   9.563  0.50 20.00           O
HETATM   47  O2 BGOL A   6      44.540  11.525   9.370  0.50 20.00           O
HETATM   48  O3 BGOL A   6      43.333  10.283  11.512  0.50 20.00           O
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
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1 2 1
HETATM    1  C1  GOL A   1       0.475   8.886  -0.730  0.50 20.00           C
HETATM    2  C2  GOL A   1       0.300  10.000   0.300  0.50 20.00           C
HETATM    3  C3  GOL A   1      -1.101   9.990   0.919  0.50 20.00           C
HETATM    4  O1  GOL A   1       0.379   7.601  -0.127  0.50 20.00           O
HETATM    5  O2  GOL A   1       0.525  11.271  -0.319  0.50 20.00           O
HETATM    6  O3  GOL A   1      -1.199  11.096   1.823  0.50 20.00           O
HETATM    7  C1  GOL A   2       6.578   9.079   4.959  1.00 20.00           C
HETATM    8  C2  GOL A   2       6.404  10.193   5.989  1.00 20.00           C
HETATM    9  C3  GOL A   2       5.003  10.183   6.608  1.00 20.00           C
HETATM   10  O1  GOL A   2       6.482   7.794   5.563  1.00 20.00           O
HETATM   11  O2  GOL A   2       6.628  11.463   5.370  1.00 20.00           O
HETATM   12  O3  GOL A   2       4.905  11.289   7.512  1.00 20.00           O
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
HETATM    1  C1 AGOL A   1       5.578   9.079   8.959  0.50 20.00           C
HETATM    2  C2 AGOL A   1       5.404  10.193   9.989  0.50 20.00           C
HETATM    3  C3 AGOL A   1       4.003  10.183  10.608  0.50 20.00           C
HETATM    4  O1 AGOL A   1       5.482   7.794   9.563  0.50 20.00           O
HETATM    5  O2 AGOL A   1       5.628  11.463   9.370  0.50 20.00           O
HETATM    6  O3 AGOL A   1       3.905  11.289  11.512  0.50 20.00           O
HETATM    7  C1 BEDO A   1       4.844  10.913   9.726  0.50 20.00           C
HETATM    8  C2 BEDO A   1       4.502   9.470   9.376  0.50 20.00           C
HETATM    9  O1 BEDO A   1       6.135  10.939  10.336  0.50 20.00           O
HETATM   10  O2 BEDO A   1       4.519   8.679  10.563  0.50 20.00           O
HETATM   11  C1 AGOL A 501      25.578   9.079   8.959  0.50 20.00           C
HETATM   12  C2 AGOL A 501      25.404  10.193   9.989  0.50 20.00           C
HETATM   13  C3 AGOL A 501      24.003  10.183  10.608  0.50 20.00           C
HETATM   14  O1 AGOL A 501      25.482   7.794   9.563  0.50 20.00           O
HETATM   15  O2 AGOL A 501      25.628  11.463   9.370  0.50 20.00           O
HETATM   16  O3 AGOL A 501      23.905  11.289  11.512  0.50 20.00           O
HETATM   17  C1 BEDO A 502      25.467  10.800   9.726  0.50 20.00           C
HETATM   18  C2 BEDO A 502      24.278   9.914   9.376  0.50 20.00           C
HETATM   19  O1 BEDO A 502      26.473   9.989  10.336  0.50 20.00           O
HETATM   20  O2 BEDO A 502      23.782   9.297  10.563  0.50 20.00           O
HETATM   21  C1 AGOL A 701      40.578   9.079   8.959  0.50 20.00           C
HETATM   22  C2 AGOL A 701      40.404  10.193   9.989  0.50 20.00           C
HETATM   23  C3 AGOL A 701      39.003  10.183  10.608  0.50 20.00           C
HETATM   24  O1 AGOL A 701      40.482   7.794   9.563  0.50 20.00           O
HETATM   25  O2 AGOL A 701      40.628  11.463   9.370  0.50 20.00           O
HETATM   26  O3 AGOL A 701      38.905  11.289  11.512  0.50 20.00           O
ATOM     27  N  ALYS A 700      38.239  14.866  10.239  0.50 20.00           N
ATOM     28  CA ALYS A 700      38.335  13.591   9.489  0.50 20.00           C
ATOM     29  C  ALYS A 700      37.032  12.796   9.649  0.50 20.00           C
ATOM     30  O  ALYS A 700      36.086  13.093  10.363  0.50 20.00           O
ATOM     31  CB ALYS A 700      39.529  12.763   9.981  0.50 20.00           C
ATOM     32  CG ALYS A 700      40.876  13.457   9.747  0.50 20.00           C
ATOM     33  CD ALYS A 700      42.042  12.570  10.193  0.50 20.00           C
ATOM     34  CE ALYS A 700      43.383  13.263   9.958  0.50 20.00           C
ATOM     35  NZ ALYS A 700      44.479  12.399  10.381  0.50 20.00           N
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
  assert acS.partner.is_ligand is True, acS.partner.is_ligand

  # same-altloc control: the only neighbour within 2 A is a LYS sidechain in
  # GOL's *own* altloc A, so it is not a partner -> still lone. (The partner
  # search no longer excludes residues by class; the altloc-match rule keeps
  # this out.)
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
HETATM    1  C1  GOL A   1      10.578  11.079  12.959  1.00 20.00           C
HETATM    2  C2  GOL A   1      10.404  12.193  13.989  1.00 20.00           C
HETATM    3  C3  GOL A   1       9.003  12.183  14.608  1.00 20.00           C
HETATM    4  O1  GOL A   1      10.482   9.794  13.563  1.00 20.00           O
HETATM    5  O2  GOL A   1      10.628  13.463  13.370  1.00 20.00           O
HETATM    6  O3  GOL A   1       8.905  13.289  15.512  1.00 20.00           O
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
HETATM    1  C1 AGOL A   1       5.578   9.079   8.959  0.50 20.00           C
HETATM    2  C2 AGOL A   1       5.404  10.193   9.989  0.50 20.00           C
HETATM    3  C3 AGOL A   1       4.003  10.183  10.608  0.50 20.00           C
HETATM    4  O1 AGOL A   1       5.482   7.794   9.563  0.50 20.00           O
HETATM    5  O2 AGOL A   1       5.628  11.463   9.370  0.50 20.00           O
HETATM    6  O3 AGOL A   1       3.905  11.289  11.512  0.50 20.00           O
HETATM    7  C1  GOL A   2      30.578   9.079   8.959  1.00 20.00           C
HETATM    8  C2  GOL A   2      30.404  10.193   9.989  1.00 20.00           C
HETATM    9  C3  GOL A   2      29.003  10.183  10.608  1.00 20.00           C
HETATM   10  O1  GOL A   2      30.482   7.794   9.563  1.00 20.00           O
HETATM   11  O2  GOL A   2      30.628  11.463   9.370  1.00 20.00           O
HETATM   12  O3  GOL A   2      28.905  11.289  11.512  1.00 20.00           O
HETATM   13  C1 BGOL A   3       6.035   9.666   8.959  0.50 20.00           C
HETATM   14  C2 BGOL A   3       5.185  10.407   9.989  0.50 20.00           C
HETATM   15  C3 BGOL A   3       4.119   9.499  10.608  0.50 20.00           C
HETATM   16  O1 BGOL A   3       6.788   8.620   9.563  0.50 20.00           O
HETATM   17  O2 BGOL A   3       4.540  11.525   9.370  0.50 20.00           O
HETATM   18  O3 BGOL A   3       3.333  10.283  11.512  0.50 20.00           O
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

_altconf_nonligand_pdb_str = '''
CRYST1   90.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  C1 AGOL A   1       5.578   9.079   8.959  0.60 20.00           C
HETATM    2  C2 AGOL A   1       5.404  10.193   9.989  0.60 20.00           C
HETATM    3  C3 AGOL A   1       4.003  10.183  10.608  0.60 20.00           C
HETATM    4  O1 AGOL A   1       5.482   7.794   9.563  0.60 20.00           O
HETATM    5  O2 AGOL A   1       5.628  11.463   9.370  0.60 20.00           O
HETATM    6  O3 AGOL A   1       3.905  11.289  11.512  0.60 20.00           O
HETATM    7  O  BHOH A 101       5.000  10.000  10.000  0.40 20.00           O
HETATM    8  C1 AGOL A   2      25.578   9.079   8.959  0.30 20.00           C
HETATM    9  C2 AGOL A   2      25.404  10.193   9.989  0.30 20.00           C
HETATM   10  C3 AGOL A   2      24.003  10.183  10.608  0.30 20.00           C
HETATM   11  O1 AGOL A   2      25.482   7.794   9.563  0.30 20.00           O
HETATM   12  O2 AGOL A   2      25.628  11.463   9.370  0.30 20.00           O
HETATM   13  O3 AGOL A   2      23.905  11.289  11.512  0.30 20.00           O
HETATM   14  O  BHOH A 102      25.000  10.000  10.000  0.30 20.00           O
HETATM   15  C1 AGOL A   3      45.578   9.079   8.959  0.50 20.00           C
HETATM   16  C2 AGOL A   3      45.404  10.193   9.989  0.50 20.00           C
HETATM   17  C3 AGOL A   3      44.003  10.183  10.608  0.50 20.00           C
HETATM   18  O1 AGOL A   3      45.482   7.794   9.563  0.50 20.00           O
HETATM   19  O2 AGOL A   3      45.628  11.463   9.370  0.50 20.00           O
HETATM   20  O3 AGOL A   3      43.905  11.289  11.512  0.50 20.00           O
HETATM   21  C1 AGOL A   4      65.578   9.079   8.959  0.50 20.00           C
HETATM   22  C2 AGOL A   4      65.404  10.193   9.989  0.50 20.00           C
HETATM   23  C3 AGOL A   4      64.003  10.183  10.608  0.50 20.00           C
HETATM   24  O1 AGOL A   4      65.482   7.794   9.563  0.50 20.00           O
HETATM   25  O2 AGOL A   4      65.628  11.463   9.370  0.50 20.00           O
HETATM   26  O3 AGOL A   4      63.905  11.289  11.512  0.50 20.00           O
ATOM     27  N  BLYS A 400      58.760  12.467   9.858  0.50 20.00           N
ATOM     28  CA BLYS A 400      58.856  11.192   9.109  0.50 20.00           C
ATOM     29  C  BLYS A 400      57.553  10.397   9.268  0.50 20.00           C
ATOM     30  O  BLYS A 400      56.607  10.694   9.983  0.50 20.00           O
ATOM     31  CB BLYS A 400      60.050  10.364   9.600  0.50 20.00           C
ATOM     32  CG BLYS A 400      61.397  11.058   9.366  0.50 20.00           C
ATOM     33  CD BLYS A 400      62.563  10.171   9.812  0.50 20.00           C
ATOM     34  CE BLYS A 400      63.904  10.864   9.577  0.50 20.00           C
ATOM     35  NZ BLYS A 400      65.000  10.000  10.000  0.50 20.00           N
'''

def run_test27():
  print('test27')
  vl = _altconf_manager(_altconf_nonligand_pdb_str)

  # (1) ligand alt conf A shares its site with a single-altloc-B water;
  #     occupancies sum to 1.0 -> recognized as split (not lone) and OK.
  ac1 = _altconf_lr(vl, 1, 'A').get_alt_conf()
  assert ac1.state == 'split_residue', ac1.state
  assert ac1.partner.resname == 'HOH', ac1.partner.resname
  assert ac1.partner.is_ligand is False, ac1.partner.is_ligand
  assert ac1.hetero is False, ac1.hetero
  assert approx_equal(ac1.occupancy.occ_sum, 1.0)
  assert ac1.flag == 'ok', ac1.flag

  # (2) same, but occupancies sum to 0.6 -> still split, flagged for occupancy.
  ac2 = _altconf_lr(vl, 2, 'A').get_alt_conf()
  assert ac2.state == 'split_residue', ac2.state
  assert ac2.partner.resname == 'HOH', ac2.partner.resname
  assert ac2.flag == 'inspect', ac2.flag
  assert 'occupancies sum to 0.60' in ac2.reason, ac2.reason

  # (3) isolated partial-occupancy ligand: nothing at the site -> genuinely lone.
  ac3 = _altconf_lr(vl, 3, 'A').get_alt_conf()
  assert ac3.state == 'lone_altloc', ac3.state
  assert ac3.partner is None
  assert ac3.flag == 'inspect', ac3.flag

  # (4) residue (LYS) partner in the complementary altloc; occ sums to 1.0 -> OK.
  ac4 = _altconf_lr(vl, 4, 'A').get_alt_conf()
  assert ac4.state == 'split_residue', ac4.state
  assert ac4.partner.resname == 'LYS', ac4.partner.resname
  assert ac4.partner.is_ligand is False, ac4.partner.is_ligand
  assert ac4.flag == 'ok', ac4.flag

  # snapshot carries the partner classification
  snap1 = _altconf_lr(vl, 1, 'A').as_picklable_snapshot()
  assert snap1.alt_conf.partner.resname == 'HOH'
  assert snap1.alt_conf.partner.is_ligand is False
  assert snap1.alt_conf.flag == 'ok'

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
