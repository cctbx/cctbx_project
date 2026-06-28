from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from libtbx.utils import null_out
from scitbx import matrix

# ------------------------------------------------------------------------------

def run():
  test_000()
  test_001()
  test_002()
  test_003()
  test_004()

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
    N-methylated (N-substituted) amino acids are tertiary amides when they are
    peptide-bonded: their backbone N already has three heavy neighbours (CA, the
    N-methyl carbon CN, and the previous residue's C), so it cannot carry a
    backbone amide H. The GeoStd monomers (e.g. MLE, MVA, SAR, BMT) still define
    an 'H' on N (valid only for the free/N-terminal form), so reduce used to add
    a spurious H. Example: the N-terminal stretch of chain C in 2oju (cyclosporin).

    The true N-terminus (DAL 1) must still get its propeller NH3, and ordinary
    secondary-amide residues (ABA 6, VAL 9) must still keep their single backbone H.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_002.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model_initial.set_stop_for_unknowns(False)

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    stop_for_unknowns = False)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  ph = model_h_added.get_hierarchy()

  # Map (resseq, resname) -> set of atom names placed in that residue.
  names_by_rg = {}
  for rg in ph.residue_groups():
    for ag in rg.atom_groups():
      key = (rg.resseq.strip(), ag.resname.strip())
      names_by_rg.setdefault(key, set()).update(
        a.name.strip() for a in ag.atoms())

  # Sanity: H placement actually ran.
  assert(model_h_added.get_hd_selection().count(True) > 0)

  # Tertiary-amide (N-methylated) residues: no backbone H on N.
  tertiary = [('2', 'MLE'), ('3', 'MLE'), ('4', 'MVA'),
              ('5', 'BMT'), ('7', 'SAR'), ('8', 'MLE')]
  for key in tertiary:
    names = names_by_rg[key]
    assert('CN' in names), '%s should be N-methylated (CN present)' % str(key)
    assert('H' not in names and 'D' not in names), \
      'spurious backbone H on tertiary amide N of %s: %s' % (str(key), sorted(names))

  # True N-terminus DAL 1 still gets its propeller NH3.
  dal = names_by_rg[('1', 'DAL')]
  for name in ('H1', 'H2', 'H3'):
    assert(name in dal), 'N-terminal propeller H %s missing on DAL 1' % name

  # Ordinary secondary-amide residues still keep their single backbone H.
  for key in [('6', 'ABA'), ('9', 'VAL')]:
    assert('H' in names_by_rg[key]), \
      'backbone amide H missing on %s' % str(key)

# ------------------------------------------------------------------------------

def test_003():
  '''
    The backbone N of pyroglutamate (PCA, a 5-membered lactam) is an amide
    nitrogen: it is bonded to CA, the carbonyl carbon CD (CD=OE), and H, and the
    monomer restrains N/CA/CG/CD/OE/H to be coplanar. Its H must therefore lie
    in that plane. The GeoStd ideal angles around N (H-N-CD=120.3, H-N-CA=114.0,
    CD-N-CA=114.2) sum to only ~348 deg, so the riding-H parameterization used to
    treat N as pyramidal ('2neigbs') and place H ~0.5 A out of plane, ignoring
    the planarity restraint. Example: N-terminal PCA of chain A in 4jp6.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_003.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model_initial.set_stop_for_unknowns(False)

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    stop_for_unknowns = False)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  ph = model_h_added.get_hierarchy()

  xyz = {}
  for a in ph.atoms():
    if a.parent().parent().resseq.strip() == '1' and a.parent().resname.strip() == 'PCA':
      xyz[a.name.strip()] = matrix.col(a.xyz)

  assert('H' in xyz), 'PCA backbone amide H was not placed'
  N, CA, CD, H = xyz['N'], xyz['CA'], xyz['CD'], xyz['H']

  # Bond length sanity.
  assert(0.8 < (H - N).length() < 1.1), 'unexpected N-H length: %.3f' % (H-N).length()

  # H must be in the amide plane: its distance from the plane through N, CA, CD
  # must be ~0 (was ~0.48 A with the pyramidal placement).
  normal = (CA - N).cross(CD - N).normalize()
  out_of_plane = abs((H - N).dot(normal))
  assert(out_of_plane < 0.05), \
    'PCA amide H is out of plane by %.3f A (should be ~0)' % out_of_plane

# ------------------------------------------------------------------------------

def test_004():
  '''
    The backbone N of an internal proline-type residue (e.g. HYP,
    4-hydroxyproline) is a tertiary amide: it is bonded to CA, the ring CD, and
    the previous residue's C. The monomer expects a backbone N-H that cannot be
    placed (the ring N is already fully substituted). reduce must report this as
    a tertiary amide, not as a generic "could not be parameterized (not enough
    restraints)" failure. Example: internal HYP of chain A in 5mas.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_004.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model_initial.set_stop_for_unknowns(False)

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    stop_for_unknowns = False)
  reduce_add_h_obj.run()
  ph = reduce_add_h_obj.get_model().get_hierarchy()

  # The HYP backbone N must not carry an H in the output.
  for ag in ph.atom_groups():
    if ag.resname.strip() == 'HYP':
      names = [a.name.strip() for a in ag.atoms()]
      assert('H' not in names), 'HYP must not have a backbone N-H: %s' % names

  # The omitted HYP H is reported as a tertiary amide, not as a generic
  # "could not be parameterized" failure.
  ta = reduce_add_h_obj.site_labels_tertiary_amide
  assert(any('HYP' in s for s in ta)), \
    'HYP backbone H not reported as tertiary amide: %s' % ta
  assert(not any('HYP' in s for s in reduce_add_h_obj.site_labels_no_para)), \
    'HYP H wrongly reported under could-not-parameterize: %s' % \
    reduce_add_h_obj.site_labels_no_para

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

pdb_str_002 = """
CRYST1   64.387   64.387  200.943  90.00  90.00 120.00 P 31 2 1
HETATM    1  N   DAL C   1      30.209  22.409  72.196  1.00 40.35           N
HETATM    2  CA  DAL C   1      30.152  23.547  73.106  1.00 39.57           C
HETATM    3  C   DAL C   1      28.727  23.723  73.642  1.00 38.87           C
HETATM    4  O   DAL C   1      28.069  22.750  74.010  1.00 38.69           O
HETATM    5  CB  DAL C   1      31.135  23.351  74.252  1.00 39.01           C
HETATM   11  N   MLE C   2      28.233  24.972  73.686  1.00 37.76           N
HETATM   12  CA  MLE C   2      26.870  25.138  74.205  1.00 38.15           C
HETATM   13  C   MLE C   2      26.011  25.920  73.220  1.00 36.94           C
HETATM   14  O   MLE C   2      25.704  27.088  73.435  1.00 38.29           O
HETATM   15  CB  MLE C   2      27.075  25.953  75.486  1.00 38.11           C
HETATM   16  CG  MLE C   2      27.879  25.267  76.592  1.00 41.61           C
HETATM   17  CD1 MLE C   2      28.104  26.251  77.738  1.00 40.41           C
HETATM   18  CD2 MLE C   2      27.140  24.012  77.083  1.00 39.25           C
HETATM   19  CN  MLE C   2      28.914  26.117  73.018  1.00 34.37           C
HETATM   34  N   MLE C   3      25.638  25.289  72.103  1.00 35.52           N
HETATM   35  CA  MLE C   3      24.798  26.122  71.239  1.00 35.63           C
HETATM   36  C   MLE C   3      23.386  25.537  71.138  1.00 35.10           C
HETATM   37  O   MLE C   3      23.133  24.408  71.581  1.00 34.67           O
HETATM   38  CB  MLE C   3      25.511  26.056  69.876  1.00 35.05           C
HETATM   39  CG  MLE C   3      26.868  26.749  69.694  1.00 37.22           C
HETATM   40  CD1 MLE C   3      27.399  26.459  68.311  1.00 36.57           C
HETATM   41  CD2 MLE C   3      26.733  28.258  69.882  1.00 40.46           C
HETATM   42  CN  MLE C   3      25.970  23.876  71.742  1.00 30.86           C
HETATM   57  N   MVA C   4      22.437  26.313  70.590  1.00 34.64           N
HETATM   58  CA  MVA C   4      21.120  25.686  70.454  1.00 34.92           C
HETATM   59  C   MVA C   4      20.694  25.787  68.991  1.00 35.62           C
HETATM   60  O   MVA C   4      20.108  26.780  68.568  1.00 36.24           O
HETATM   61  CB  MVA C   4      20.018  26.263  71.366  1.00 35.21           C
HETATM   62  CG1 MVA C   4      20.511  26.363  72.815  1.00 32.56           C
HETATM   63  CG2 MVA C   4      18.785  25.369  71.291  1.00 31.42           C
HETATM   64  CN  MVA C   4      22.646  27.721  70.170  1.00 29.91           C
HETATM   77  N   BMT C   5      20.987  24.755  68.195  1.00 36.70           N
HETATM   78  CA  BMT C   5      20.584  24.894  66.789  1.00 37.85           C
HETATM   79  C   BMT C   5      19.453  23.972  66.354  1.00 37.25           C
HETATM   80  O   BMT C   5      18.969  23.160  67.141  1.00 38.13           O
HETATM   81  CB  BMT C   5      21.819  24.724  65.804  1.00 37.56           C
HETATM   82  OG1 BMT C   5      22.218  23.344  65.803  1.00 38.73           O
HETATM   83  CG2 BMT C   5      23.058  25.605  66.152  1.00 36.67           C
HETATM   84  CD1 BMT C   5      22.671  27.101  66.346  1.00 31.96           C
HETATM   85  CD2 BMT C   5      24.268  25.464  65.175  1.00 40.19           C
HETATM   86  CE  BMT C   5      24.634  26.661  64.279  1.00 42.17           C
HETATM   87  CZ  BMT C   5      24.227  26.734  62.914  1.00 46.10           C
HETATM   88  CH  BMT C   5      24.602  27.952  62.060  1.00 43.62           C
HETATM   89  CN  BMT C   5      21.555  23.487  68.730  1.00 36.97           C
HETATM  108  N   ABA C   6      19.047  24.098  65.091  1.00 38.16           N
HETATM  109  CA  ABA C   6      17.952  23.298  64.545  1.00 40.43           C
HETATM  110  C   ABA C   6      18.396  22.524  63.303  1.00 41.22           C
HETATM  111  O   ABA C   6      18.013  22.864  62.180  1.00 41.72           O
HETATM  112  CB  ABA C   6      16.719  24.148  64.224  1.00 39.98           C
HETATM  113  CG  ABA C   6      16.196  24.964  65.401  1.00 39.91           C
HETATM  122  N   SAR C   7      19.194  21.458  63.490  1.00 40.14           N
HETATM  123  CA  SAR C   7      19.660  20.675  62.345  1.00 41.02           C
HETATM  124  C   SAR C   7      21.133  20.272  62.494  1.00 41.75           C
HETATM  125  O   SAR C   7      21.429  19.109  62.739  1.00 42.59           O
HETATM  126  CN  SAR C   7      19.640  20.993  64.828  1.00 35.07           C
HETATM  133  N   MLE C   8      22.070  21.230  62.332  1.00 42.04           N
HETATM  134  CA  MLE C   8      23.452  20.760  62.543  1.00 41.88           C
HETATM  135  C   MLE C   8      23.889  21.053  63.975  1.00 40.61           C
HETATM  136  O   MLE C   8      24.313  22.167  64.284  1.00 40.93           O
HETATM  137  CB  MLE C   8      24.289  21.592  61.551  1.00 42.91           C
HETATM  138  CG  MLE C   8      25.549  20.947  60.934  1.00 45.33           C
HETATM  139  CD1 MLE C   8      26.543  22.040  60.543  1.00 44.56           C
HETATM  140  CD2 MLE C   8      26.214  19.993  61.911  1.00 45.50           C
HETATM  141  CN  MLE C   8      21.764  22.672  61.984  1.00 41.22           C
ATOM    156  N   VAL C   9      23.780  20.056  64.845  1.00 39.54           N
ATOM    157  CA  VAL C   9      24.160  20.220  66.246  1.00 41.09           C
ATOM    158  C   VAL C   9      25.679  20.332  66.350  1.00 41.25           C
ATOM    159  O   VAL C   9      26.402  19.504  65.798  1.00 43.29           O
ATOM    160  CB  VAL C   9      23.694  19.009  67.102  1.00 41.32           C
ATOM    161  CG1 VAL C   9      23.802  19.340  68.588  1.00 39.65           C
ATOM    162  CG2 VAL C   9      22.267  18.632  66.739  1.00 40.83           C
TER
END
"""

# ------------------------------------------------------------------------------

pdb_str_003 = """
CRYST1   27.688   54.971   73.509  90.00  90.00 90.00 P 21 21 21
HETATM    1  N   PCA A   1     -11.673  -9.095 -16.308  1.00  8.92           N
HETATM    2  CA  PCA A   1     -10.819  -9.389 -15.195  1.00  8.33           C
HETATM    3  C   PCA A   1     -10.830  -8.243 -14.205  1.00  7.30           C
HETATM    4  O   PCA A   1     -11.068  -7.106 -14.550  1.00  7.71           O
HETATM    5  CB  PCA A   1      -9.428  -9.570 -15.897  1.00  9.07           C
HETATM    6  CG  PCA A   1      -9.711  -9.650 -17.380  1.00 10.60           C
HETATM    7  CD  PCA A   1     -11.141  -9.250 -17.510  1.00  9.08           C
HETATM    8  OE  PCA A   1     -11.716  -9.140 -18.582  1.00 11.01           O
ATOM      9  N   SER A   2     -10.385  -8.545 -13.002  1.00  8.03           N
ATOM     10  CA  SER A   2     -10.366  -7.552 -11.919  1.00  8.31           C
ATOM     11  C   SER A   2      -9.236  -7.832 -10.951  1.00  6.76           C
ATOM     12  O   SER A   2      -8.826  -8.990 -10.799  1.00  7.56           O
ATOM     13  CB  SER A   2     -11.691  -7.490 -11.164  1.00  9.75           C
ATOM     14  OG  SER A   2     -11.883  -8.665 -10.484  1.00 13.96           O
TER
END
"""

# ------------------------------------------------------------------------------

pdb_str_004 = """
CRYST1   27.688   54.971   73.509  90.00  90.00 90.00 P 21 21 21
HETATM   54  N   AIB A   9       8.193   6.231   1.668  1.00  2.65           N
HETATM   55  CA  AIB A   9       8.830   5.063   2.289  1.00  2.66           C
HETATM   56  C   AIB A   9      10.204   5.417   2.905  1.00  2.37           C
HETATM   57  O   AIB A   9      10.503   4.963   4.008  1.00  2.49           O
HETATM   58  CB1 AIB A   9       9.066   3.999   1.207  1.00  2.69           C
HETATM   59  CB2 AIB A   9       7.887   4.519   3.361  1.00  2.79           C
HETATM   60  N   HYP A  10      11.088   6.150   2.203  1.00  2.24           N
HETATM   61  CA  HYP A  10      12.443   6.314   2.751  1.00  2.55           C
HETATM   62  C   HYP A  10      12.513   7.043   4.090  1.00  2.82           C
HETATM   63  O   HYP A  10      13.527   6.972   4.759  1.00  3.62           O
HETATM   64  CB  HYP A  10      13.173   7.123   1.668  1.00  2.69           C
HETATM   65  CG  HYP A  10      12.373   6.906   0.400  1.00  2.36           C
HETATM   66  CD  HYP A  10      10.962   6.808   0.895  1.00  2.35           C
HETATM   67  OD1 HYP A  10      12.699   5.681  -0.242  1.00  2.81           O
ATOM     68  N   GLN A  11      11.449   7.823   4.417  1.00  2.57           N
ATOM     69  CA  GLN A  11      11.425   8.570   5.666  1.00  2.59           C
ATOM     70  C   GLN A  11      10.952   7.707   6.834  1.00  2.68           C
ATOM     71  O   GLN A  11      11.109   8.109   7.988  1.00  3.44           O
ATOM     72  CB  GLN A  11      10.543   9.792   5.569  1.00  2.68           C
ATOM     73  CG  GLN A  11      10.831  10.646   4.324  1.00  3.13           C
ATOM     74  CD  GLN A  11      12.289  10.903   4.096  1.00  2.99           C
ATOM     75  OE1 GLN A  11      13.051  11.164   5.025  1.00  3.82           O
ATOM     76  NE2 GLN A  11      12.692  10.845   2.804  1.00  3.56           N
TER
END
"""

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
