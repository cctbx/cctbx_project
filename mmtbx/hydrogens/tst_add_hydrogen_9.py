from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.model
from scitbx import matrix
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models

def run():
  test_001()
  test_002()
  test_003()
  test_004()
  test_005()
  test_006()

# ------------------------------------------------------------------------------

def test_001():
  '''
    CD1 is missing --> H, HD1, HE1 can't be parameterized --> should not be added.
  '''
  compare_models(pdb_str      = pdb_str_001,
                 not_contains = ' HE1')

# ------------------------------------------------------------------------------

def test_002():
  '''
    N-terminal residue with two altloc conformers --> both conformers receive
    NH3 propeller (H1, H2, H3); neither retains the standard backbone 'H'.
    Regression test for shadowed nested loop in place_n_terminal_propeller
    that only removed the placeholder 'H' from the first altloc, leaving the
    other altloc with [H, H2, H3] -- a combination pdb_interpretation cannot
    parameterize, so H2 and H3 were silently dropped.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_002.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model = model.select(~model.get_hd_selection())
  add_h = reduce_hydrogen.place_hydrogens(model=model)
  add_h.run()
  model_h = add_h.get_model()

  altloc_h_names = {}
  rg = model_h.get_hierarchy().models()[0].chains()[0].residue_groups()[0]
  assert rg.resseq_as_int() == 1
  for ag in rg.atom_groups():
    if ag.altloc == '': continue
    altloc_h_names[ag.altloc] = set(
      a.name.strip() for a in ag.atoms() if a.element.strip() in ('H', 'D'))

  assert set(altloc_h_names) == {'A', 'B'}, \
    "expected altlocs A and B, got %s" % sorted(altloc_h_names)
  for alt in ('A', 'B'):
    for h in ('H1', 'H2', 'H3'):
      assert h in altloc_h_names[alt], \
        "altloc %s SER 1 is missing %s; got %s" % (
          alt, h, sorted(altloc_h_names[alt]))
    assert 'H' not in altloc_h_names[alt], \
      "altloc %s SER 1 retained plain backbone 'H' instead of NH3 propeller" % alt

# ------------------------------------------------------------------------------

def test_003():
  '''
    A nucleotide with a missing dihedral.
  '''
  compare_models(pdb_str = pdb_str_003)

# ------------------------------------------------------------------------------

def test_004():
  '''
    Glycine cross-linked to an adamantane scaffold via CA-C1 (auto-detected
    LINK at ~1.53 A). exclude_H_on_links must drop the HA atom occupying the
    bonded position; the other HA must survive. resseq=10 -> no NH3 propeller.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_004.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model = model.select(~model.get_hd_selection())
  add_h = reduce_hydrogen.place_hydrogens(model=model)
  add_h.run()
  names = set(a.name.strip() for a in add_h.get_model().get_hierarchy().atoms()
              if a.parent().resname.strip() == 'GLY')
  assert 'H'   in names, "missing backbone amide H on GLY"
  assert 'HA2' in names, "GLY HA2 should survive (opposite face from C1[ADM])"
  assert 'HA3' not in names, "GLY HA3 should be removed by exclude_H_on_links"
  assert 'H1'  not in names, "no NH3 propeller expected at resseq=10"

# ------------------------------------------------------------------------------

def test_005():
  '''
    Cyclohexyl ligand (CHX) bridges CYS-SG and LYS-NZ via auto-detected
    C5-SG (~1.77 A) and C6-NZ (~1.43 A) links. exclude_H_on_links must drop
    CHX H52 (C5 side facing SG) and H62 (C6 side facing NZ); H51 and H61
    must survive on the opposite faces.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_005.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model = model.select(~model.get_hd_selection())
  add_h = reduce_hydrogen.place_hydrogens(model=model)
  add_h.run()
  chx_names = set(
    a.name.strip() for a in add_h.get_model().get_hierarchy().atoms()
    if a.parent().resname.strip() == 'CHX')
  assert 'H51' in chx_names, "CHX H51 should survive (opposite face of C5-SG link)"
  assert 'H61' in chx_names, "CHX H61 should survive (opposite face of C6-NZ link)"
  assert 'H52' not in chx_names, "CHX H52 should be removed (C5-SG link)"
  assert 'H62' not in chx_names, "CHX H62 should be removed (C6-NZ link)"

# ------------------------------------------------------------------------------

def test_006():
  '''
    B0I makes two thioether links (6b17): H00 SD - B0I CB1 and CYS SG - B0I CB2
    (both auto-detected at ~1.8 A). CB1 and CB2 are methyls in the free ligand;
    each link removes one H, leaving a CH2. exclude_H_on_links must leave the
    two surviving H at distinct tetrahedral positions -- not collapse both onto
    a single coincident site (the bug this guards: the kept-H repositioning
    placed every kept H of a 2-heavy-neighbour atom at the same in-plane
    bisector position).
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_006.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.set_stop_for_unknowns(False)
  model = model.select(~model.get_hd_selection())
  add_h = reduce_hydrogen.place_hydrogens(model=model, stop_for_unknowns=False)
  add_h.run()
  mh = add_h.get_model()
  atoms = mh.get_hierarchy().atoms()
  grm = mh.get_restraints_manager().geometry
  bond_proxies, asu = grm.get_all_bond_proxies(sites_cart=mh.get_sites_cart())

  def neighbors(resn, nm):
    iseq = None
    for a in atoms:
      if a.parent().resname.strip() == resn and a.name.strip() == nm:
        iseq = a.i_seq
        break
    assert iseq is not None, '%s %s not found' % (resn, nm)
    out = set()
    for p in bond_proxies:
      i, j = p.i_seqs
      if   i == iseq: out.add(j)
      elif j == iseq: out.add(i)
    return out

  def angle(i, j, k):
    return (matrix.col(atoms[i].xyz) - matrix.col(atoms[j].xyz)).angle(
            matrix.col(atoms[k].xyz) - matrix.col(atoms[j].xyz), deg=True)

  for cb, link_partner in (('CB1', 'SD'), ('CB2', 'SG')):
    c_iseq = [a.i_seq for a in atoms
              if a.parent().resname.strip() == 'B0I' and a.name.strip() == cb][0]
    nb = neighbors('B0I', cb)
    # the thioether link must be present, otherwise the scenario isn't exercised
    assert link_partner in set(atoms[k].name.strip() for k in nb), \
      'expected %s-%s thioether link' % (cb, link_partner)
    h_iseqs   = [k for k in nb if atoms[k].element.strip() in ('H', 'D')]
    heavy_nb  = [k for k in nb if atoms[k].element.strip() not in ('H', 'D')]
    assert len(h_iseqs) == 2, \
      'linked %s (a CH2) should keep 2 H, got %s' % (
        cb, sorted(atoms[k].name.strip() for k in h_iseqs))
    d = (matrix.col(atoms[h_iseqs[0]].xyz) -
         matrix.col(atoms[h_iseqs[1]].xyz)).length()
    assert d > 0.5, \
      'the two H on linked %s are coincident (d=%.3f A)' % (cb, d)
    # The two H must sit tetrahedrally about the two real heavy neighbours
    # (the C-C and C-S directions), not eclipse one of them. Every
    # heavy-C-H angle should be near tetrahedral, never collapsed toward 0.
    for h in h_iseqs:
      for hv in heavy_nb:
        a = angle(hv, c_iseq, h)
        assert 95.0 < a < 125.0, \
          '%s: %s-%s-%s angle = %.1f (not tetrahedral)' % (
            cb, atoms[hv].name.strip(), cb, atoms[h].name.strip(), a)

# ------------------------------------------------------------------------------

pdb_str_001 = """
REMARK CD1 is missing --> HD1, HE1 cannot not be placed; H is not placed
CRYST1   17.955   13.272   13.095  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A 139      10.241   7.920   5.000  1.00 10.00           N
ATOM      2  CA  TYR A 139      10.853   7.555   6.271  1.00 10.00           C
ATOM      3  C   TYR A 139      12.362   7.771   6.227  1.00 10.00           C
ATOM      4  O   TYR A 139      12.955   8.272   7.181  1.00 10.00           O
ATOM      5  CB  TYR A 139      10.540   6.098   6.617  1.00 10.00           C
ATOM      6  CG  TYR A 139       9.063   5.805   6.749  1.00 10.00           C
ATOM      7  CD2 TYR A 139       8.414   5.943   7.969  1.00 10.00           C
ATOM      8  CE1 TYR A 139       6.966   5.122   5.770  1.00 10.00           C
ATOM      9  CE2 TYR A 139       7.064   5.676   8.095  1.00 10.00           C
ATOM     10  CZ  TYR A 139       6.345   5.266   6.993  1.00 10.00           C
ATOM     11  OH  TYR A 139       5.000   5.000   7.113  1.00 10.00           O
ATOM     12  H   TYR A 139       9.382   7.879   5.001  1.00 10.00           H
ATOM     13  HA  TYR A 139      10.487   8.115   6.973  1.00 10.00           H
ATOM     14  HB2 TYR A 139      10.961   5.881   7.464  1.00 10.00           H
ATOM     15  HB3 TYR A 139      10.893   5.529   5.916  1.00 10.00           H
ATOM     16  HD2 TYR A 139       8.896   6.220   8.714  1.00 10.00           H
ATOM     17  HE2 TYR A 139       6.643   5.772   8.919  1.00 10.00           H
ATOM     18  HH  TYR A 139       4.752   5.127   7.905  1.00 10.00           H
TER
"""

pdb_str_002 = """
REMARK N-terminal SER in double conformation (altlocs A/B) -- both should get NH3
CRYST1   46.893   46.017   95.993  90.00  98.40  90.00 P 1 21 1
ATOM      1  C   SER A   1     -12.619  34.880  11.188  1.00 18.10           C
ATOM      2  O   SER A   1     -12.707  34.745   9.948  1.00 22.52           O
ATOM      3  N  ASER A   1     -14.822  36.058  11.168  0.55 24.23           N
ATOM      4  CA ASER A   1     -13.798  35.463  11.990  0.55 23.39           C
ATOM      5  CB ASER A   1     -14.561  34.377  12.750  0.55 24.81           C
ATOM      6  OG ASER A   1     -15.918  34.652  13.022  0.55 34.66           O
ATOM      7  N  BSER A   1     -14.115  36.726  11.048  0.46 26.00           N
ATOM      8  CA BSER A   1     -13.579  35.763  11.997  0.46 23.88           C
ATOM      9  CB BSER A   1     -14.690  34.985  12.696  0.46 23.97           C
ATOM     10  OG BSER A   1     -15.372  34.116  11.815  0.46 29.41           O
ATOM     11  N   TYR A   2     -11.649  34.403  11.891  1.00 10.98           N
ATOM     12  CA  TYR A   2     -10.600  33.615  11.277  1.00  9.36           C
ATOM     13  C   TYR A   2     -11.075  32.217  10.979  1.00  8.66           C
ATOM     14  O   TYR A   2     -11.965  31.658  11.655  1.00  9.72           O
ATOM     15  CB  TYR A   2      -9.378  33.538  12.232  1.00  8.96           C
ATOM     16  CG  TYR A   2      -8.734  34.820  12.606  1.00  8.86           C
ATOM     17  CD1 TYR A   2      -8.761  35.972  11.813  1.00  9.47           C
ATOM     18  CD2 TYR A   2      -7.981  34.882  13.783  1.00  9.74           C
ATOM     19  CE1 TYR A   2      -8.113  37.147  12.185  1.00  9.54           C
ATOM     20  CE2 TYR A   2      -7.335  36.044  14.145  1.00  9.97           C
ATOM     21  CZ  TYR A   2      -7.406  37.182  13.362  1.00  8.39           C
ATOM     22  OH  TYR A   2      -6.760  38.306  13.784  1.00  9.22           O
ATOM     23  N   THR A   3     -10.473  31.593   9.961  1.00  8.54           N
ATOM     24  CA  THR A   3     -10.697  30.229   9.657  1.00  8.82           C
ATOM     25  C   THR A   3      -9.356  29.515   9.445  1.00  7.80           C
ATOM     26  O   THR A   3      -8.339  30.115   9.091  1.00  8.51           O
ATOM     27  CB  THR A   3     -11.537  30.085   8.363  1.00 11.47           C
ATOM     28  OG1 THR A   3     -10.781  30.625   7.267  1.00 13.80           O
ATOM     29  CG2 THR A   3     -12.878  30.751   8.531  1.00 16.31           C
TER
"""

pdb_str_003 = '''
CRYST1   60.683   61.851   76.893  90.00  90.00  90.00 P 21 21 21
SCALE1      0.016479  0.000000  0.000000        0.00000
SCALE2      0.000000  0.016168  0.000000        0.00000
SCALE3      0.000000  0.000000  0.013005        0.00000
HETATM    1  C1' ADP A1311      -7.459  14.326  10.821  0.60 14.02           C
HETATM    2  C2  ADP A1311      -5.449  12.545   7.224  0.60 15.85           C
HETATM    3  C2' ADP A1311      -7.594  15.611  11.604  0.60 14.04           C
HETATM    4  C3' ADP A1311      -8.097  15.126  12.927  0.60 13.75           C
HETATM    5  C4  ADP A1311      -6.903  13.885   8.453  0.60 14.26           C
HETATM    6  C4' ADP A1311      -8.964  13.965  12.513  0.60 14.88           C
HETATM    7  C5  ADP A1311      -7.414  14.430   7.198  0.60 14.79           C
HETATM    8  C5' ADP A1311     -10.403  14.410  12.351  0.60 12.89           C
HETATM    9  C6  ADP A1311      -6.877  13.954   5.942  0.60 15.37           C
HETATM   10  C8  ADP A1311      -8.458  15.356   8.818  0.60 13.69           C
HETATM   11  N1  ADP A1311      -5.907  13.025   6.029  0.60 16.08           N
HETATM   12  N3  ADP A1311      -5.925  12.949   8.425  0.60 15.93           N
HETATM   13  N6  ADP A1311      -7.358  14.471   4.786  0.60 16.26           N
HETATM   14  N7  ADP A1311      -8.348  15.322   7.484  0.60 15.22           N
HETATM   15  N9  ADP A1311      -7.603  14.521   9.383  0.60 14.10           N
HETATM   16  O1A ADP A1311     -12.593  15.064  10.083  0.60  6.78           O
HETATM   17  O1B ADP A1311     -11.480  19.808  11.509  0.30 13.50           O
HETATM   18  O2' ADP A1311      -6.345  16.239  11.737  0.60 12.16           O
HETATM   19  O2A ADP A1311     -12.829  15.834  12.325  0.60  8.51           O
HETATM   20  O2B ADP A1311     -10.576  17.778  12.692  0.30 17.96           O
HETATM   21  O3' ADP A1311      -7.034  14.659  13.737  0.60 15.01           O
HETATM   22  O3A ADP A1311     -11.739  17.443  10.617  0.60 12.19           O
HETATM   23  O3B ADP A1311     -13.078  18.233  12.535  0.20 17.22           O
HETATM   24  O4' ADP A1311      -8.522  13.493  11.245  0.60 15.01           O
HETATM   25  O5' ADP A1311     -10.557  15.534  11.511  0.60 12.05           O
HETATM   26  PA  ADP A1311     -12.004  15.966  11.098  0.60 10.28           P
HETATM   27  PB  ADP A1311     -11.709  18.390  11.914  0.30 15.89           P
HETATM   28  H2  ADP A1311      -4.662  11.801   7.209  0.60 15.85           H
HETATM   29  H1' ADP A1311      -6.457  13.916  11.009  0.60 14.02           H
HETATM   30  H2' ADP A1311      -8.248  16.358  11.133  0.60 14.04           H
HETATM   31  H3' ADP A1311      -8.605  15.901  13.517  0.60 13.75           H
HETATM   32  H4' ADP A1311      -8.894  13.181  13.280  0.60 14.88           H
HETATM   33  H8  ADP A1311      -9.154  15.983   9.361  0.60 13.69           H
HETATM   34 H5'1 ADP A1311     -10.983  13.581  11.940  0.60 12.89           H
HETATM   35 H5'2 ADP A1311     -10.810  14.648  13.336  0.60 12.89           H
HETATM   36 HN61 ADP A1311      -6.985  14.160   3.900  0.60 16.26           H
HETATM   37 HN62 ADP A1311      -8.089  15.167   4.810  0.60 16.26           H
HETATM   38 HO2' ADP A1311      -6.432  17.010  12.315  0.60 12.16           H
HETATM   39 HO3' ADP A1311      -6.479  15.405  14.002  0.60 15.01           H
'''

pdb_str_004 = """
REMARK GLY 10 with CA cross-linked to adamantane C1 (auto-link, ~1.53 A)
CRYST1   30.746   56.560   99.515  90.00  90.00  90.00 P 21 21 21
SCALE1      0.032525  0.000000  0.000000        0.00000
SCALE2      0.000000  0.017680  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010049        0.00000
ATOM      1  N   GLY B  10       2.957 -10.749  -1.364  1.00 22.36           N
ATOM      2  CA  GLY B  10       2.830 -10.787  -2.810  1.00 23.42           C
ATOM      3  C   GLY B  10       3.228 -12.114  -3.425  1.00 21.76           C
ATOM      4  O   GLY B  10       2.672 -12.516  -4.453  1.00 20.78           O
TER
HETATM    5  C1  ADM B 101       3.506  -9.585  -3.485  1.00 29.99           C
HETATM    6  C10 ADM B 101       3.272  -6.994  -4.853  1.00 33.86           C
HETATM    7  C2  ADM B 101       2.997  -9.452  -4.919  1.00 37.40           C
HETATM    8  C3  ADM B 101       3.647  -8.260  -5.608  1.00 30.02           C
HETATM    9  C4  ADM B 101       5.159  -8.444  -5.610  1.00 40.96           C
HETATM   10  C5  ADM B 101       5.668  -8.541  -4.180  1.00 34.29           C
HETATM   11  C6  ADM B 101       5.301  -7.254  -3.457  1.00 29.01           C
HETATM   12  C7  ADM B 101       3.789  -7.091  -3.424  1.00 35.16           C
HETATM   13  C8  ADM B 101       3.171  -8.301  -2.730  1.00 32.12           C
HETATM   14  C9  ADM B 101       5.021  -9.747  -3.503  1.00 25.37           C
HETATM   15  H3  ADM B 101       3.284  -8.191  -6.643  1.00 36.03           H
HETATM   16  H21 ADM B 101       1.912  -9.324  -4.912  1.00 44.89           H
HETATM   17  H22 ADM B 101       3.225 -10.365  -5.473  1.00 44.89           H
HETATM   18  H41 ADM B 101       5.419  -9.353  -6.156  1.00 49.15           H
HETATM   19  H42 ADM B 101       5.633  -7.598  -6.112  1.00 49.15           H
HETATM   20  H5  ADM B 101       6.760  -8.663  -4.185  1.00 41.15           H
HETATM   21  H61 ADM B 101       5.752  -6.402  -3.971  1.00 34.81           H
HETATM   22  H62 ADM B 101       5.690  -7.280  -2.437  1.00 34.81           H
HETATM   23  H7  ADM B 101       3.528  -6.176  -2.874  1.00 42.19           H
HETATM   24  H81 ADM B 101       3.549  -8.368  -1.708  1.00 38.54           H
HETATM   25  H82 ADM B 101       2.087  -8.178  -2.681  1.00 38.54           H
HETATM   26  H91 ADM B 101       5.287 -10.658  -4.043  1.00 30.45           H
HETATM   27  H92 ADM B 101       5.394  -9.838  -2.480  1.00 30.45           H
HETATM   28 H101 ADM B 101       2.187  -6.873  -4.847  1.00 40.63           H
HETATM   29 H102 ADM B 101       3.711  -6.124  -5.347  1.00 40.63           H
END
"""

pdb_str_005 = """
REMARK CHX cyclohexyl ligand bridges CYS-SG (via C5) and LYS-NZ (via C6)
CRYST1   69.873   71.696  120.736  90.00  90.00  90.00 P 21 21 21
SCALE1      0.014312  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013948  0.000000        0.00000
SCALE3      0.000000  0.000000  0.008283        0.00000
ATOM      1  N   CYS W   0      20.663  29.193 -38.814  1.00 11.26           N
ATOM      2  CA  CYS W   0      21.799  29.059 -37.916  1.00 10.56           C
ATOM      3  C   CYS W   0      23.068  28.775 -38.699  1.00 11.18           C
ATOM      4  O   CYS W   0      23.053  28.093 -39.724  1.00 10.81           O
ATOM      5  CB  CYS W   0      21.598  27.851 -37.004  1.00  9.35           C
ATOM      6  SG  CYS W   0      20.153  27.941 -35.989  1.00 14.74           S
ATOM      7  N   LYS W   2      26.671  27.556 -37.018  1.00 13.10           N
ATOM      8  CA  LYS W   2      27.479  27.245 -35.835  1.00 15.84           C
ATOM      9  C   LYS W   2      28.394  26.064 -36.128  1.00 14.04           C
ATOM     10  O   LYS W   2      27.924  25.027 -36.678  1.00 14.05           O
ATOM     11  CB  LYS W   2      26.565  26.937 -34.640  1.00 11.96           C
ATOM     12  CG  LYS W   2      25.577  28.094 -34.319  1.00 12.87           C
ATOM     13  CD  LYS W   2      24.634  27.641 -33.213  1.00 17.44           C
ATOM     14  CE  LYS W   2      23.617  28.719 -32.809  1.00 13.32           C
ATOM     15  NZ  LYS W   2      22.583  28.813 -33.871  1.00 25.21           N
TER
HETATM   16  C1  CHX W 101      20.190  28.826 -32.556  1.00 27.31           C
HETATM   17  C2  CHX W 101      18.722  28.775 -32.978  1.00 31.27           C
HETATM   18  C3  CHX W 101      18.316  29.955 -33.845  1.00 40.64           C
HETATM   19  C4  CHX W 101      19.302  30.261 -34.967  1.00 24.46           C
HETATM   20  C5  CHX W 101      20.557  29.381 -35.046  1.00 15.31           C
HETATM   21  C6  CHX W 101      21.170  28.958 -33.715  1.00 27.48           C
END
"""

# 6b17 chain E: H00 (modified Met, SD) and a CYS each thioether-link to the two
# methyl carbons CB1/CB2 of the bridging ligand B0I.
pdb_str_006 = """
CRYST1   36.758   35.583   38.076  90.00 118.88  90.00 P 1 21 1
HETATM  209  CA  H00 E   1       6.717 -12.704   5.860  1.00 16.88           C
HETATM  210  CB  H00 E   1       6.318 -14.063   5.292  1.00 16.74           C
HETATM  211  CG  H00 E   1       6.017 -15.116   6.354  1.00 16.68           C
HETATM  212  SD  H00 E   1       5.293 -16.601   5.615  1.00 16.31           S
HETATM  213  C   H00 E   1       7.931 -12.777   6.784  1.00 18.18           C
HETATM  214  O   H00 E   1       9.010 -13.154   6.318  1.00 20.57           O
ATOM    309  N   CYS E  14      -0.395 -23.942   7.399  1.00 13.33           N
ATOM    310  CA  CYS E  14      -0.510 -23.103   6.204  1.00 14.58           C
ATOM    311  C   CYS E  14       0.780 -23.163   5.401  1.00 16.24           C
ATOM    312  O   CYS E  14       0.762 -23.109   4.165  1.00 17.78           O
ATOM    313  CB  CYS E  14      -0.779 -21.656   6.597  1.00 14.48           C
ATOM    314  SG  CYS E  14      -2.347 -21.437   7.427  1.00 14.92           S
HETATM  645  CB1 B0I E 101       3.656 -16.038   5.070  1.00 13.50           C
HETATM  646  CG1 B0I E 101       2.744 -15.629   6.189  1.00 12.92           C
HETATM  647  CD1 B0I E 101       1.936 -16.583   6.802  1.00 12.07           C
HETATM  648  CD2 B0I E 101       2.687 -14.302   6.601  1.00 12.54           C
HETATM  649  CE1 B0I E 101       1.087 -16.213   7.847  1.00 12.06           C
HETATM  650  CE2 B0I E 101       1.830 -13.932   7.638  1.00 12.24           C
HETATM  651  CZ1 B0I E 101       1.031 -14.885   8.264  1.00 12.17           C
HETATM  652  CB2 B0I E 101      -2.466 -19.639   7.341  1.00 13.09           C
HETATM  653  CG2 B0I E 101      -1.533 -18.870   8.236  1.00 12.42           C
HETATM  654  CD3 B0I E 101      -0.661 -17.942   7.685  1.00 11.87           C
HETATM  655  CD4 B0I E 101      -1.557 -19.074   9.612  1.00 12.57           C
HETATM  656  CE3 B0I E 101       0.198 -17.211   8.495  1.00 11.38           C
HETATM  657  CE4 B0I E 101      -0.697 -18.346  10.435  1.00 12.86           C
HETATM  658  CZ2 B0I E 101       0.174 -17.408   9.876  1.00 12.11           C
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
