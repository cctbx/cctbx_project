from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from libtbx.utils import null_out

# ------------------------------------------------------------------------------

def run():
  test_000()
  test_001()
  test_002()

# ------------------------------------------------------------------------------

def test_000():
  '''
    An unknown ligand (FCO, no CCD/GeoStd/user restraints) gets a throwaway
    restraint dictionary auto-generated during H placement so that pdb
    interpretation can build it. That dictionary must not leak out as a
    permanent restraint object on the returned model: it is purpose-built for
    placement (bonds shortened to 0.9x, idealized esd=1/period=1 torsions from
    the CCD conformer) and would produce bogus geometry outliers if re-used by
    downstream validation (e.g. mmtbx.development.validate_ligands).
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_000.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  assert(model_initial.get_hd_selection().count(True) == 0)

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    stop_for_unknowns = False)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()

  # Sanity: H placement ran fully (H atoms placed on the LEU residue).
  assert(model_h_added.get_hd_selection().count(True) > 0)

  # The auto-generated FCO placement restraints must not persist on the model.
  ro = model_h_added.get_restraint_objects()
  names = [] if ro is None else [name for name, _ in ro]
  assert('auto_FCO' not in names), \
    'auto-generated placement restraints leaked: %s' % names

  # Behavioral check: a downstream consumer that re-derives restraints from the
  # returned model the way validate_ligands does (harvest restraint objects,
  # set them, re-process) must not synthesize bogus FCO geometry. Re-processing
  # rebuilds mon_lib_srv from the restraint objects only, so with the leak fixed
  # the unknown ligand stays unknown and gets no restraints.
  ro = model_h_added.get_restraint_objects()
  model_h_added.set_restraint_objects([] if ro is None else list(ro))
  model_h_added.set_stop_for_unknowns(False)
  model_h_added.process(make_restraints=True)
  isel = model_h_added.iselection('resname FCO and not (element H or element D)')
  fco = model_h_added.select(isel)
  grm = fco.get_restraints_manager().geometry
  assert(grm.dihedral_proxies.size() == 0), \
    'bogus FCO dihedral restraints present: %d' % grm.dihedral_proxies.size()

# ------------------------------------------------------------------------------

def test_001():
  '''
    N-terminal propeller H atoms (H1/H2/H3) must inherit the hetero flag of
    their parent residue. For a modified amino acid written as HETATM (e.g.
    N-terminal MSE), the terminal hydrogens were emitted as ATOM while the rest
    of the residue stayed HETATM, producing a residue with mixed ATOM/HETATM
    records that breaks downstream parsers (e.g. Biopython mmCIF).
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_001.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model_initial.set_stop_for_unknowns(False)

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    stop_for_unknowns = False)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  ph = model_h_added.get_hierarchy()

  # Collect the MSE residue atoms.
  mse_atoms = {}
  for a in ph.atoms():
    if a.parent().resname.strip() == 'MSE':
      mse_atoms[a.name.strip()] = a

  # Sanity: propeller H were placed.
  for name in ('H1', 'H2', 'H3'):
    assert(name in mse_atoms), 'expected propeller H %s not placed' % name

  # The parent residue is HETATM, so every atom (including the propeller H)
  # must carry hetero=True. No mixing of ATOM/HETATM within one residue.
  heteros = set(a.hetero for a in mse_atoms.values())
  assert(heteros == set([True])), \
    'mixed ATOM/HETATM within MSE: ' + \
    ', '.join('%s=%s' % (n, a.hetero) for n, a in sorted(mse_atoms.items()))

# ------------------------------------------------------------------------------

def test_002():
  '''
    Cyclic peptide (3njw): the N-terminal GLY 1 is linked (lactam) to the CG of
    ASP 9, so its N is a secondary amide, not a free N-terminus. It must get a
    single backbone H, not an N-terminal H1/H2/H3 propeller. (Without the link
    to ASP 9 the same GLY would get the full propeller.)
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_002.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model_initial.set_stop_for_unknowns(False)

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    stop_for_unknowns = False)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  atoms = model_h_added.get_hierarchy().atoms()

  # Locate the N atom of the N-terminal GLY 1.
  n_iseq = None
  for a in atoms:
    ag = a.parent()
    if (ag.resname.strip() == 'GLY' and
        ag.parent().resseq_as_int() == 1 and
        a.name.strip() == 'N'):
      n_iseq = a.i_seq
      break
  assert(n_iseq is not None), 'GLY 1 N not found'

  # Atoms bonded to that N, via the geometry restraints.
  grm = model_h_added.get_restraints_manager().geometry
  bond_proxies, asu = grm.get_all_bond_proxies(
    sites_cart = model_h_added.get_sites_cart())
  bonded = set()
  for p in bond_proxies:
    i, j = p.i_seqs
    if   i == n_iseq: bonded.add(j)
    elif j == n_iseq: bonded.add(i)
  partners = set(atoms[k].name.strip() for k in bonded)

  # The lactam link to ASP 9 must be present (CG), otherwise the scenario under
  # test (linked N-terminus) is not being exercised.
  assert('CG' in partners), 'expected lactam link N(GLY1)-CG(ASP9), got %s' % \
    sorted(partners)

  # Exactly one H on N: a single backbone H, no NH3 propeller.
  n_hydrogens = sorted(atoms[k].name.strip() for k in bonded
                       if atoms[k].element.strip() in ('H', 'D'))
  assert(len(n_hydrogens) == 1), \
    'N-terminal linked GLY 1 should have a single N-H, got %s' % n_hydrogens

# ------------------------------------------------------------------------------

pdb_str_000 = """
CRYST1  100.667  101.210  170.826  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   LEU A 482     114.924  99.962 -27.431  1.00 10.00           N
ATOM      2  CA  LEU A 482     114.188 101.177 -27.843  1.00 10.00           C
ATOM      3  C   LEU A 482     115.223 102.183 -28.350  1.00 10.00           C
ATOM      4  O   LEU A 482     116.313 102.353 -27.708  1.00 10.00           O
ATOM      5  CB  LEU A 482     113.411 101.726 -26.625  1.00 10.00           C
ATOM      6  CG  LEU A 482     112.764 103.127 -26.762  1.00 10.00           C
ATOM      7  CD1 LEU A 482     111.736 103.154 -27.867  1.00 10.00           C
ATOM      8  CD2 LEU A 482     112.154 103.534 -25.429  1.00 10.00           C
TER
HETATM    9  C1  FCO A 601     107.099  98.547 -26.445  1.00 10.00           C
HETATM   10  C2  FCO A 601     107.980  98.742 -23.921  1.00 10.00           C
HETATM   11  C3  FCO A 601     108.613 100.597 -25.504  1.00 10.00           C
HETATM   12  N1  FCO A 601     107.159  97.752 -27.314  1.00 10.00           N
HETATM   13  N2  FCO A 601     108.687  98.106 -23.190  1.00 10.00           N
HETATM   14  O3  FCO A 601     109.639 101.050 -25.723  1.00 10.00           O
HETATM   15 FE   FCO A 601     107.039  99.845 -25.099  1.00 10.00          FE
END
"""

# ------------------------------------------------------------------------------

pdb_str_001 = """
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1
HETATM    1  N   MSE A   1      36.554 -26.382 -54.813  1.00129.74           N
HETATM    2  CA  MSE A   1      37.413 -25.165 -54.854  1.00129.96           C
HETATM    3  C   MSE A   1      38.272 -25.033 -53.595  1.00128.03           C
HETATM    4  O   MSE A   1      37.746 -24.824 -52.499  1.00128.43           O
HETATM    5  CB  MSE A   1      36.535 -23.924 -55.017  1.00133.18           C
HETATM    6  CG  MSE A   1      37.311 -22.627 -55.097  1.00137.42           C
HETATM    7 SE   MSE A   1      36.236 -21.220 -55.847  1.00142.16          SE
HETATM    8  CE  MSE A   1      36.573 -21.565 -57.723  1.00140.74           C
ATOM      9  N   ALA A   2      39.580 -25.160 -53.700  1.00128.00           N
ATOM     10  CA  ALA A   2      40.450 -25.050 -52.540  1.00128.00           C
ATOM     11  C   ALA A   2      41.900 -25.250 -52.950  1.00128.00           C
ATOM     12  O   ALA A   2      42.250 -25.100 -54.120  1.00128.00           O
ATOM     13  CB  ALA A   2      40.050 -26.080 -51.500  1.00128.00           C
END
"""

# ------------------------------------------------------------------------------

# 3njw, residues 1-9: a cyclic peptide whose N-terminal GLY 1 is linked
# (lactam) to the CG of ASP 9.
pdb_str_002 = """
CRYST1   19.465   21.432   29.523  90.00  90.00  90.00 P 21 21 21
LINK         N   GLY A   1                 CG  ASP A   9     1555   1555  1.34
ATOM      1  N   GLY A   1       6.011  23.726   5.538  1.00  4.36           N
ATOM      2  CA  GLY A   1       7.279  24.418   5.504  1.00  4.79           C
ATOM      3  C   GLY A   1       8.370  23.751   6.291  1.00  4.66           C
ATOM      4  O   GLY A   1       9.449  24.344   6.484  1.00  5.69           O
ATOM      5  N   LEU A   2       8.134  22.534   6.793  1.00  4.84           N
ATOM      6  CA  LEU A   2       9.036  21.914   7.745  1.00  5.30           C
ATOM      7  C   LEU A   2       9.922  20.861   7.111  1.00  5.17           C
ATOM      8  O   LEU A   2       9.581  20.226   6.109  1.00  5.18           O
ATOM      9  CB  LEU A   2       8.200  21.239   8.862  1.00  5.04           C
ATOM     10  CG  LEU A   2       7.320  22.154   9.683  1.00  6.43           C
ATOM     11  CD1 LEU A   2       6.376  21.399  10.529  1.00  7.04           C
ATOM     12  CD2 LEU A   2       8.211  23.064  10.553  1.00 10.71           C
ATOM     13  N   PRO A   3      11.058  20.597   7.758  1.00  5.43           N
ATOM     14  CA  PRO A   3      11.927  19.506   7.347  1.00  6.16           C
ATOM     15  C   PRO A   3      11.539  18.165   7.931  1.00  5.34           C
ATOM     16  O   PRO A   3      12.394  17.258   8.065  1.00  6.40           O
ATOM     17  CB  PRO A   3      13.314  19.996   7.864  1.00  8.27           C
ATOM     18  CG  PRO A   3      13.002  20.748   9.093  1.00  8.11           C
ATOM     19  CD  PRO A   3      11.699  21.426   8.805  1.00  6.80           C
ATOM     20  N   TRP A   4      10.285  17.990   8.316  1.00  4.64           N
ATOM     21  CA  TRP A   4       9.755  16.714   8.712  1.00  4.53           C
ATOM     22  C   TRP A   4       8.273  16.693   8.400  1.00  4.33           C
ATOM     23  O   TRP A   4       7.623  17.716   8.189  1.00  4.81           O
ATOM     24  CB  TRP A   4      10.002  16.386  10.189  1.00  5.02           C
ATOM     25  CG  TRP A   4       9.274  17.201  11.154  1.00  5.11           C
ATOM     26  CD1 TRP A   4       8.113  16.856  11.754  1.00  5.90           C
ATOM     27  CD2 TRP A   4       9.617  18.500  11.725  1.00  5.51           C
ATOM     28  NE1 TRP A   4       7.686  17.813  12.601  1.00  6.14           N
ATOM     29  CE2 TRP A   4       8.602  18.828  12.634  1.00  5.65           C
ATOM     30  CE3 TRP A   4      10.738  19.304  11.615  1.00  6.51           C
ATOM     31  CZ2 TRP A   4       8.656  19.968  13.385  1.00  6.75           C
ATOM     32  CZ3 TRP A   4      10.784  20.456  12.384  1.00  7.66           C
ATOM     33  CH2 TRP A   4       9.736  20.801  13.263  1.00  8.03           C
ATOM     34  N   GLY A   5       7.735  15.465   8.408  1.00  4.59           N
ATOM     35  CA  GLY A   5       6.363  15.239   8.055  1.00  4.78           C
ATOM     36  C   GLY A   5       6.228  14.644   6.688  1.00  4.30           C
ATOM     37  O   GLY A   5       7.176  14.072   6.142  1.00  5.39           O
ATOM     38  N   CYS A   6       5.038  14.805   6.111  1.00  4.95           N
ATOM     39  CA  CYS A   6       4.737  14.291   4.792  1.00  5.27           C
ATOM     40  C   CYS A   6       4.677  15.413   3.789  1.00  5.02           C
ATOM     41  O   CYS A   6       4.405  16.592   4.151  1.00  4.93           O
ATOM     42  CB  CYS A   6       3.416  13.515   4.742  1.00  5.80           C
ATOM     43  SG  CYS A   6       3.392  12.095   5.898  1.00  6.51           S
ATOM     44  N   PRO A   7       4.842  15.126   2.498  1.00  6.21           N
ATOM     45  CA  PRO A   7       4.686  16.174   1.493  1.00  6.22           C
ATOM     46  C   PRO A   7       3.217  16.586   1.450  1.00  5.81           C
ATOM     47  O   PRO A   7       2.354  15.687   1.449  1.00  7.43           O
ATOM     48  CB  PRO A   7       5.106  15.456   0.205  1.00  8.78           C
ATOM     49  CG  PRO A   7       5.972  14.380   0.602  1.00 11.03           C
ATOM     50  CD  PRO A   7       5.337  13.848   1.938  1.00  8.80           C
ATOM     51  N   SER A   8       2.868  17.871   1.378  1.00  5.98           N
ATOM     52  CA  SER A   8       3.706  19.022   1.408  1.00  6.08           C
ATOM     53  C   SER A   8       3.018  20.110   2.257  1.00  5.61           C
ATOM     54  O   SER A   8       1.815  20.173   2.316  1.00  7.41           O
ATOM     55  CB  SER A   8       3.993  19.575   0.044  1.00  8.00           C
ATOM     56  OG  SER A   8       4.762  18.631  -0.709  1.00 10.23           O
ATOM     57  N   ASP A   9       3.861  20.962   2.856  1.00  4.41           N
ATOM     58  CA  ASP A   9       3.324  22.217   3.387  1.00  4.26           C
ATOM     59  C   ASP A   9       3.240  23.233   2.275  1.00  4.25           C
ATOM     60  O   ASP A   9       3.655  23.005   1.142  1.00  4.88           O
ATOM     61  CB  ASP A   9       4.061  22.659   4.637  1.00  3.84           C
ATOM     62  CG  ASP A   9       5.454  23.195   4.444  1.00  3.99           C
ATOM     63  OD1 ASP A   9       6.011  23.127   3.350  1.00  4.95           O
END
"""

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
