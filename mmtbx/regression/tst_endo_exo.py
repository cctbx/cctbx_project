"""Regression test: ``mmtbx.development.endo_exo`` produces a stable
QM region around two contrasting Fe sites:

* **1BQ8** (Pyrococcus rubredoxin, P 21 21 21).  All four Cys ligands
  live inside the ASU, so symmetry plays no role.  Exercises the
  baseline BFS + capping pipeline.

* **2C2U** (Halobacterium ferritin core, P 2 3).  The Fe sits on a
  3-fold special position; three symmetry-related copies of Asp 93 and
  HOH 2154 from the ASU map onto the metal.  Exercises:

  - symmetry-aware adjacency (the BFS picks up the Asp 93 / HOH 2154
    symmetry images even though the parent atoms are far from the
    metal in the ASU);
  - special-position deduplication (three sym_ops put Fe at the same
    physical point; only one Fe atom survives in the materialized
    region);
  - per-sym_op chain IDs (the identity copy stays in chain "A"; the
    two extra images land in single-character chain IDs).
"""

from __future__ import absolute_import, division, print_function

import io
from collections import defaultdict

import iotbx.pdb
import libtbx.phil
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times

import mmtbx.model
from iotbx.data_manager import DataManager
from mmtbx.regression import model_1yjp
from mmtbx.programs.endo_exo import Program as EndoexoProgram
from scitbx import matrix
from cctbx import geometry_restraints, sgtbx

from mmtbx.geometry_restraints.endo_exo.util import _canon_op
from mmtbx.geometry_restraints.endo_exo.capping import HydrogenCapper
from mmtbx.geometry_restraints.endo_exo.cutting import BondCutDetector
from mmtbx.geometry_restraints.endo_exo.graph import AtomGraphBuilder
from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower


# 8 A sphere around the Fe of 1BQ8 (29 residues, 154 atoms).  Slightly
# larger than the QM region endo_exo extracts (72 atoms) so the BFS has
# scaffold to cap into.
_1BQ8_FE_SPHERE_PDB = """\
CRYST1   33.823   34.705   43.204  90.00  90.00  90.00 P 21 21 21
SCALE1      0.029566  0.000000  0.000000        0.00000
SCALE2      0.000000  0.028814  0.000000        0.00000
SCALE3      0.000000  0.000000  0.023146        0.00000
ATOM     37  N   VAL A   5      16.394   5.495   0.266  1.00  3.20           N
ATOM     38  CA  VAL A   5      15.980   6.892   0.493  1.00  3.49           C
ATOM     39  C   VAL A   5      14.480   7.032   0.278  1.00  2.98           C
ATOM     40  O   VAL A   5      13.922   6.570  -0.728  1.00  3.30           O
ATOM     41  CB  VAL A   5      16.763   7.863  -0.366  1.00  3.68           C
ATOM     42  CG1 VAL A   5      16.514   7.724  -1.870  1.00  4.62           C
ATOM     43  CG2 VAL A   5      16.523   9.338   0.046  1.00  5.32           C
ATOM     44  N   CYS A   6      13.876   7.840   1.115  1.00  2.56           N
ATOM     45  CA  CYS A   6      12.515   8.342   0.910  1.00  2.77           C
ATOM     46  C   CYS A   6      12.644   9.458  -0.142  1.00  2.56           C
ATOM     47  O   CYS A   6      13.271  10.473   0.122  1.00  3.23           O
ATOM     48  CB  CYS A   6      11.937   8.897   2.187  1.00  2.83           C
ATOM     49  SG  CYS A   6      10.256   9.599   1.876  1.00  3.04           S
ATOM     50  N   LYS A   7      12.120   9.227  -1.346  1.00  2.98           N
ATOM     51  CA  LYS A   7      12.260  10.201  -2.421  1.00  3.51           C
ATOM     52  C   LYS A   7      11.541  11.504  -2.103  1.00  3.88           C
ATOM     53  O   LYS A   7      11.857  12.540  -2.672  1.00  5.65           O
ATOM     54  CB  LYS A   7      11.742   9.650  -3.758  1.00  5.11           C
ATOM     55  CG  LYS A   7      12.524   8.451  -4.265  1.00  6.96           C
ATOM     56  CD  LYS A   7      11.938   7.989  -5.583  1.00  8.84           C
ATOM     57  CE  LYS A   7      10.500   7.483  -5.529  1.00  9.73           C
ATOM     58  NZ  LYS A   7      10.079   6.961  -6.833  1.00 11.95           N
ATOM     59  N   ILE A   8      10.572  11.487  -1.162  1.00  3.53           N
ATOM     60  CA  ILE A   8       9.858  12.717  -0.799  1.00  3.44           C
ATOM     61  C   ILE A   8      10.661  13.609   0.123  1.00  3.26           C
ATOM     62  O   ILE A   8      10.797  14.826  -0.111  1.00  5.13           O
ATOM     63  CB  ILE A   8       8.487  12.375  -0.162  1.00  5.11           C
ATOM     64  CG1 ILE A   8       7.671  11.360  -0.956  1.00  5.47           C
ATOM     65  CG2 ILE A   8       7.719  13.652   0.102  1.00  5.65           C
ATOM     66  CD1 ILE A   8       7.349  11.801  -2.358  1.00  8.43           C
ATOM     67  N   CYS A   9      11.210  13.074   1.223  1.00  3.35           N
ATOM     68  CA  CYS A   9      11.769  13.851   2.272  1.00  3.49           C
ATOM     69  C   CYS A   9      13.219  13.660   2.581  1.00  3.53           C
ATOM     70  O   CYS A   9      13.796  14.406   3.390  1.00  3.48           O
ATOM     71  CB  CYS A   9      10.948  13.686   3.564  1.00  3.16           C
ATOM     72  SG  CYS A   9      11.311  12.173   4.470  1.00  3.34           S
ATOM     73  N   GLY A  10      13.879  12.653   2.014  1.00  2.95           N
ATOM     74  CA  GLY A  10      15.282  12.416   2.247  1.00  2.72           C
ATOM     75  C   GLY A  10      15.614  11.546   3.426  1.00  2.73           C
ATOM     76  O   GLY A  10      16.820  11.254   3.626  1.00  3.00           O
ATOM     77  N   TYR A  11      14.667  11.081   4.202  1.00  2.81           N
ATOM     78  CA  TYR A  11      14.915  10.112   5.274  1.00  2.91           C
ATOM     79  C   TYR A  11      15.553   8.845   4.649  1.00  2.90           C
ATOM     80  O   TYR A  11      15.159   8.433   3.551  1.00  3.10           O
ATOM     81  CB  TYR A  11      13.580   9.729   5.930  1.00  2.88           C
ATOM     82  CG  TYR A  11      13.703   8.536   6.859  1.00  2.57           C
ATOM     83  CD1 TYR A  11      14.232   8.652   8.139  1.00  3.41           C
ATOM     84  CD2 TYR A  11      13.305   7.275   6.410  1.00  3.49           C
ATOM     85  CE1 TYR A  11      14.336   7.547   8.970  1.00  4.29           C
ATOM     86  CE2 TYR A  11      13.427   6.178   7.235  1.00  3.09           C
ATOM     87  CZ  TYR A  11      13.946   6.307   8.493  1.00  3.73           C
ATOM     88  OH  TYR A  11      14.053   5.156   9.266  1.00  4.47           O
ATOM     89  N   ILE A  12      16.530   8.325   5.360  1.00  2.50           N
ATOM     90  CA  ILE A  12      17.176   7.071   4.961  1.00  2.77           C
ATOM     91  C   ILE A  12      16.631   5.931   5.778  1.00  2.84           C
ATOM     92  O   ILE A  12      16.815   5.923   6.989  1.00  3.43           O
ATOM     93  CB  ILE A  12      18.702   7.193   5.103  1.00  3.75           C
ATOM     94  CG1 ILE A  12      19.290   8.421   4.373  1.00  4.12           C
ATOM     95  CG2 ILE A  12      19.388   5.922   4.674  1.00  4.75           C
ATOM     96  CD1 ILE A  12      19.194   8.339   2.859  1.00  6.27           C
ATOM    279  N   TRP A  37       3.798  -0.020   6.436  1.00  4.44           N
ATOM    280  CA  TRP A  37       4.860   0.808   7.049  1.00  3.68           C
ATOM    281  C   TRP A  37       5.032   2.052   6.206  1.00  3.09           C
ATOM    282  O   TRP A  37       4.931   2.011   5.003  1.00  4.07           O
ATOM    283  CB  TRP A  37       6.150  -0.001   7.149  1.00  3.69           C
ATOM    284  CG  TRP A  37       7.317   0.766   7.681  1.00  3.57           C
ATOM    285  CD1 TRP A  37       7.768   0.734   8.982  1.00  3.65           C
ATOM    286  CD2 TRP A  37       8.215   1.617   6.994  1.00  3.73           C
ATOM    287  NE1 TRP A  37       8.896   1.532   9.118  1.00  3.98           N
ATOM    288  CE2 TRP A  37       9.187   2.061   7.896  1.00  3.73           C
ATOM    289  CE3 TRP A  37       8.319   2.027   5.639  1.00  3.78           C
ATOM    290  CZ2 TRP A  37      10.240   2.931   7.536  1.00  3.88           C
ATOM    291  CZ3 TRP A  37       9.358   2.853   5.292  1.00  4.00           C
ATOM    292  CH2 TRP A  37      10.308   3.297   6.238  1.00  4.03           C
ATOM    293  N   VAL A  38       5.270   3.185   6.882  1.00  2.88           N
ATOM    294  CA  VAL A  38       5.473   4.448   6.200  1.00  3.57           C
ATOM    295  C   VAL A  38       6.750   5.130   6.648  1.00  3.48           C
ATOM    296  O   VAL A  38       7.242   4.891   7.767  1.00  3.30           O
ATOM    297  CB  VAL A  38       4.286   5.394   6.306  1.00  4.43           C
ATOM    298  CG1 VAL A  38       3.017   4.799   5.695  1.00  5.23           C
ATOM    299  CG2 VAL A  38       4.031   5.795   7.745  1.00  5.51           C
ATOM    300  N   CYS A  39       7.214   6.059   5.822  1.00  3.08           N
ATOM    301  CA  CYS A  39       8.351   6.911   6.228  1.00  3.01           C
ATOM    302  C   CYS A  39       8.015   7.538   7.573  1.00  2.69           C
ATOM    303  O   CYS A  39       6.957   8.180   7.696  1.00  2.90           O
ATOM    304  CB  CYS A  39       8.518   8.006   5.167  1.00  2.63           C
ATOM    305  SG  CYS A  39       9.909   9.101   5.614  1.00  3.02           S
ATOM    306  N   PRO A  40       8.872   7.421   8.563  1.00  3.42           N
ATOM    307  CA  PRO A  40       8.561   7.954   9.916  1.00  3.79           C
ATOM    308  C   PRO A  40       8.613   9.480   9.954  1.00  3.62           C
ATOM    309  O   PRO A  40       8.148  10.060  10.929  1.00  5.67           O
ATOM    310  CB  PRO A  40       9.655   7.382  10.800  1.00  4.40           C
ATOM    311  CG  PRO A  40      10.750   7.051   9.885  1.00  5.10           C
ATOM    312  CD  PRO A  40      10.139   6.613   8.599  1.00  3.49           C
ATOM    313  N   ILE A  41       9.158  10.104   8.922  1.00  3.58           N
ATOM    314  CA  ILE A  41       9.257  11.549   8.867  1.00  3.55           C
ATOM    315  C   ILE A  41       8.095  12.177   8.150  1.00  3.59           C
ATOM    316  O   ILE A  41       7.399  13.021   8.678  1.00  4.62           O
ATOM    317  CB  ILE A  41      10.614  11.965   8.241  1.00  3.97           C
ATOM    318  CG1 ILE A  41      11.762  11.284   8.954  1.00  4.60           C
ATOM    319  CG2 ILE A  41      10.737  13.460   8.152  1.00  5.47           C
ATOM    320  CD1 ILE A  41      11.735  11.381  10.447  1.00  7.33           C
ATOM    321  N   CYS A  42       7.817  11.737   6.928  1.00  3.39           N
ATOM    322  CA  CYS A  42       6.813  12.362   6.092  1.00  2.32           C
ATOM    323  C   CYS A  42       5.546  11.599   5.836  1.00  2.80           C
ATOM    324  O   CYS A  42       4.605  12.109   5.233  1.00  3.24           O
ATOM    325  CB  CYS A  42       7.416  12.822   4.772  1.00  3.36           C
ATOM    326  SG  CYS A  42       7.709  11.447   3.585  1.00  3.47           S
ATOM    327  N   GLY A  43       5.529  10.292   6.224  1.00  2.45           N
ATOM    328  CA  GLY A  43       4.397   9.433   6.041  1.00  3.04           C
ATOM    329  C   GLY A  43       4.241   8.817   4.697  1.00  3.00           C
ATOM    330  O   GLY A  43       3.225   8.146   4.401  1.00  4.16           O
ATOM    331  N   ALA A  44       5.237   8.944   3.824  1.00  3.45           N
ATOM    332  CA  ALA A  44       5.165   8.344   2.481  1.00  3.40           C
ATOM    333  C   ALA A  44       5.111   6.836   2.552  1.00  3.56           C
ATOM    334  O   ALA A  44       5.790   6.203   3.365  1.00  3.89           O
ATOM    335  CB  ALA A  44       6.378   8.747   1.654  1.00  3.81           C
ATOM    336  N   PRO A  45       4.346   6.215   1.643  1.00  3.69           N
ATOM    337  CA  PRO A  45       4.298   4.754   1.585  1.00  4.02           C
ATOM    338  C   PRO A  45       5.576   4.171   1.101  1.00  3.00           C
ATOM    339  O   PRO A  45       6.468   4.880   0.549  1.00  3.12           O
ATOM    340  CB  PRO A  45       3.135   4.457   0.624  1.00  5.15           C
ATOM    341  CG  PRO A  45       2.917   5.686  -0.103  1.00  6.41           C
ATOM    342  CD  PRO A  45       3.406   6.842   0.655  1.00  4.15           C
ATOM    358  N   GLU A  48       7.481   4.933  -2.272  1.00  2.79           N
ATOM    359  CA  GLU A  48       8.221   6.151  -2.488  1.00  3.20           C
ATOM    360  C   GLU A  48       9.676   6.047  -2.095  1.00  2.89           C
ATOM    361  O   GLU A  48      10.378   7.082  -2.151  1.00  3.92           O
ATOM    362  CB  GLU A  48       7.533   7.364  -1.863  1.00  3.74           C
ATOM    363  CG  GLU A  48       6.066   7.505  -2.237  1.00  3.93           C
ATOM    364  CD  GLU A  48       5.870   7.509  -3.737  1.00  3.85           C
ATOM    365  OE1 GLU A  48       6.135   8.569  -4.357  1.00  5.88           O
ATOM    366  OE2 GLU A  48       5.533   6.437  -4.316  1.00  4.39           O
ATOM    367  N   PHE A  49      10.156   4.875  -1.748  1.00  2.54           N
ATOM    368  CA  PHE A  49      11.543   4.634  -1.421  1.00  3.16           C
ATOM    369  C   PHE A  49      12.287   4.035  -2.598  1.00  3.23           C
ATOM    370  O   PHE A  49      11.739   3.260  -3.399  1.00  4.57           O
ATOM    371  CB  PHE A  49      11.623   3.646  -0.221  1.00  3.07           C
ATOM    372  CG  PHE A  49      11.333   4.340   1.105  1.00  2.65           C
ATOM    373  CD1 PHE A  49      10.072   4.713   1.442  1.00  2.91           C
ATOM    374  CD2 PHE A  49      12.383   4.683   1.940  1.00  2.86           C
ATOM    375  CE1 PHE A  49       9.830   5.384   2.660  1.00  3.94           C
ATOM    376  CE2 PHE A  49      12.152   5.386   3.132  1.00  3.21           C
ATOM    377  CZ  PHE A  49      10.853   5.721   3.494  1.00  3.61           C
TER
HETATM  422 FE    FE A  55       9.794  10.608   3.873  1.00  2.90          FE
HETATM  429  O   HOH A 107       5.079   8.927   9.538  1.00 12.57           O
HETATM  445  O   HOH A 123       9.120  16.310   5.902  1.00 10.75           O
HETATM  446  O   HOH A 124      12.913  16.679   4.634  1.00  5.39           O
HETATM  470  O   HOH A 148      13.863  16.664   0.689  1.00 19.89           O
HETATM  480  O   HOH A 158      13.763  14.830  -1.239  1.00 16.61           O
HETATM  486  O   HOH A 164       8.536  15.979   3.044  1.00  8.33           O
HETATM  505  O   HOH A 183      11.838  17.652   2.302  1.00 20.40           O
HETATM  559  O   HOH A 307       9.914  17.434   0.789  0.66 11.76           O
HETATM  604  O   HOH A 422       9.389  18.060   1.780  0.34  9.73           O
"""


def _run_endo_exo_on_string(pdb_str, radius=None, include=None, selection=None):
  """Drive ``mmtbx.programs.endo_exo.Program`` in-memory on a PDB string
  with default settings (metal scan, radius=5.0, depth=3).  Parses the
  string with ``iotbx.pdb`` directly -- no disk roundtrip.  Returns
  the single result dict produced for the seed.  *radius* overrides
  ``params.buffer.radius`` when given; *include*, when given, is a
  ``(selection, scope, proximity)`` tuple configuring
  ``residues_to_include``; *selection*, when given, seeds from that CCTBX
  selection string instead of scanning for metals."""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  dm = DataManager(["model"])
  dm.add_model("model", model)
  dm.set_default_model("model")

  master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
  params = master.extract()
  params.write_files = False
  if radius is not None:
    params.buffer.radius = radius
  if selection is not None:
    params.selection = [selection]
  if include is not None:
    (params.residues_to_include.selection,
     params.residues_to_include.scope,
     params.residues_to_include.proximity) = include

  prog = EndoexoProgram(dm, params, master_phil=master, logger=io.StringIO())
  prog.validate()
  prog.run()
  results = prog.get_results()
  assert len(results) == 1, (
    f"expected 1 submodel (1 seed atom in input); got {len(results)}")
  return results[0]


def _residue_atom_names(result, resseq):
  """Return the sorted set of atom names for residue *resseq* in the
  materialized submodel of *result* (empty set if the residue is absent)."""
  names = set()
  for rg in result["model"].get_hierarchy().residue_groups():
    if rg.resseq.strip() == resseq:
      for ag in rg.atom_groups():
        for a in ag.atoms():
          names.add(a.name.strip())
  return names


def exercise_residues_to_include():
  """residues_to_include pulls a residue into the region whole, exempt from
  the sidechain cut rules, gated by the per_seed proximity sphere.

  Target: Lys 7 of 1BQ8.  Its closest atom sits 5.88 A from the Fe, so it
  is absent from the default region and straddles the proximity threshold
  -- ideal for exercising both the gate and the include path."""
  full_lys7 = {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"}

  # Baseline: Lys 7 is not in the default region.
  base = _run_endo_exo_on_string(_1BQ8_FE_SPHERE_PDB)
  assert base["model"].get_number_of_atoms() == 72
  assert _residue_atom_names(base, "7") == set()

  # per_seed, proximity below the 5.88 A closest approach -> still excluded.
  near_excl = _run_endo_exo_on_string(
    _1BQ8_FE_SPHERE_PDB, include=("resseq 7", "per_seed", 5.0))
  assert near_excl["model"].get_number_of_atoms() == 72
  assert _residue_atom_names(near_excl, "7") == set()

  # per_seed, proximity above 5.88 A -> included whole (all 9 heavy atoms),
  # regardless of the sidechain cut rules.
  near_incl = _run_endo_exo_on_string(
    _1BQ8_FE_SPHERE_PDB, include=("resseq 7", "per_seed", 7.0))
  assert near_incl["model"].get_number_of_atoms() == 77
  assert _residue_atom_names(near_incl, "7") == full_lys7

  # global ignores proximity: included even with a tiny sphere.
  glob = _run_endo_exo_on_string(
    _1BQ8_FE_SPHERE_PDB, include=("resseq 7", "global", 1.0))
  assert glob["model"].get_number_of_atoms() == 77
  assert _residue_atom_names(glob, "7") == full_lys7

  # Whole-residue expansion: a single-atom selection still pulls the
  # complete residue.
  expand = _run_endo_exo_on_string(
    _1BQ8_FE_SPHERE_PDB, include=("resseq 7 and name NZ", "global", 1.0))
  assert _residue_atom_names(expand, "7") == full_lys7


def exercise_submodel_shape():
  """Lock in the structural shape of the submodel: total atom count,
  element distribution, seed and cap iseq counts."""
  result = _run_endo_exo_on_string(_1BQ8_FE_SPHERE_PDB)

  atoms = list(result["model"].get_hierarchy().atoms())
  assert len(atoms) == 72, (
    f"submodel atom count drifted: expected 72, got {len(atoms)}")

  elements = {}
  for a in atoms:
    el = a.element.strip().upper()
    elements[el] = elements.get(el, 0) + 1
  expected_elements = {"C": 40, "FE": 1, "H": 10, "N": 8, "O": 9, "S": 4}
  assert elements == expected_elements, (
    f"submodel element distribution drifted:\n"
    f"  expected: {expected_elements}\n"
    f"  got     : {elements}")

  assert len(result["seed_iseqs"]) == 1
  assert len(result["cap_iseqs"]) == 10
  seed_iseq = result["seed_iseqs"][0]
  assert atoms[seed_iseq].element.strip().upper() == "FE"

  # Each cap's heavy-atom anchor is recorded (caps sit 1.1 A from their
  # anchor): anchors are heavy atoms, and every cap has one within bonding
  # distance.
  anchors = result["cap_anchor_iseqs"]
  assert anchors, "no cap anchors recorded"
  for ai in anchors:
    assert atoms[ai].element.strip().upper() != "H", (
      f"cap anchor {ai} is not a heavy atom")
  for ci in result["cap_iseqs"]:
    cap_xyz = matrix.col(atoms[ci].xyz)
    assert any((matrix.col(atoms[ai].xyz) - cap_xyz).length() < 1.2
               for ai in anchors), (
      f"cap {ci} has no recorded anchor within bonding distance")


def exercise_cys_coordination():
  """Verify the chemistry: the four Cys side chains coordinate the Fe
  with Sg atoms within 3 A of the seed."""
  result = _run_endo_exo_on_string(_1BQ8_FE_SPHERE_PDB)
  hier = result["model"].get_hierarchy()
  atoms = list(hier.atoms())

  cys_residues = [
    ag for ag in hier.atom_groups()
    if ag.resname.strip().upper() == "CYS"]
  assert len(cys_residues) == 4, (
    f"expected 4 Cys residues in submodel; got {len(cys_residues)}")

  seed_iseq = result["seed_iseqs"][0]
  fe_xyz = matrix.col(atoms[seed_iseq].xyz)
  sg_distances = []
  for ag in cys_residues:
    sg = next((a for a in ag.atoms()
               if a.name.strip().upper() == "SG"), None)
    assert sg is not None, (
      f"Cys residue {ag.parent().resseq.strip()} has no SG atom")
    d = (matrix.col(sg.xyz) - fe_xyz).length()
    sg_distances.append(d)
    assert d < 3.0, (
      f"Cys-SG -> Fe distance {d:.2f} A is outside coordination range")

  mean_d = sum(sg_distances) / len(sg_distances)
  assert approx_equal(mean_d, 2.27, eps=0.1)


# Asp 93 / Arg 89 / Arg 92 / Phe 90 + Fe 1210 + three Fe-shell waters
# from 2C2U (HsFt-like ferritin core, P 2 3, a=b=c=90.369).  Fe sits on
# a 3-fold special position: three symmetry operators (1_555, 6_566 and
# 12_665) map ASU Asp 93 / HOH 2154 onto the metal at 2.5-2.6 A, so the
# materialized QM region must contain three Asp 93 copies and three HOH
# 2154 copies plus a single deduplicated Fe.  Fe occupancy is 0.5 in
# 2C2U; bumped to 1.00 here so the BFS metal-detection threshold (which
# uses crystallographic occupancy) is not the thing being exercised.
_2C2U_FE_SPHERE_PDB = """\
CRYST1   90.369   90.369   90.369  90.00  90.00  90.00 P 2 3
SCALE1      0.011066  0.000000  0.000000        0.00000
SCALE2      0.000000  0.011066  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011066        0.00000
ATOM    481  N   ARG A  89      29.832  67.632  19.953  1.00  6.64           N
ATOM    482  CA  ARG A  89      29.359  67.682  21.326  1.00  7.99           C
ATOM    483  C   ARG A  89      29.123  66.317  21.969  1.00  8.16           C
ATOM    484  O   ARG A  89      28.978  66.252  23.181  1.00  8.85           O
ATOM    485  CB  ARG A  89      28.123  68.590  21.438  1.00  8.08           C
ATOM    486  CG  ARG A  89      26.860  67.986  20.809  1.00  9.06           C
ATOM    487  CD  ARG A  89      25.899  69.081  20.404  1.00  9.62           C
ATOM    488  NE  ARG A  89      24.608  68.586  19.936  1.00  9.83           N
ATOM    489  CZ  ARG A  89      24.381  68.082  18.723  1.00  9.19           C
ATOM    490  NH1 ARG A  89      25.372  68.013  17.839  1.00 10.52           N
ATOM    491  NH2 ARG A  89      23.190  67.671  18.374  1.00 11.31           N
ATOM    492  N   PHE A  90      29.143  65.230  21.179  1.00  7.41           N
ATOM    493  CA  PHE A  90      29.067  63.858  21.670  1.00  8.07           C
ATOM    494  C   PHE A  90      30.440  63.196  21.800  1.00  6.75           C
ATOM    495  O   PHE A  90      30.528  62.012  22.041  1.00  6.99           O
ATOM    496  CB  PHE A  90      28.107  63.005  20.854  1.00  9.51           C
ATOM    497  CG  PHE A  90      26.747  63.660  20.727  1.00  9.88           C
ATOM    498  CD1 PHE A  90      25.917  63.757  21.804  1.00 12.35           C
ATOM    499  CD2 PHE A  90      26.356  64.247  19.523  1.00 11.28           C
ATOM    500  CE1 PHE A  90      24.652  64.408  21.670  1.00 13.72           C
ATOM    501  CE2 PHE A  90      25.134  64.895  19.415  1.00 13.59           C
ATOM    502  CZ  PHE A  90      24.329  64.988  20.488  1.00 12.95           C
ATOM    514  N   ARG A  92      32.664  63.114  24.285  1.00  6.21           N
ATOM    515  CA  ARG A  92      32.906  62.301  25.482  1.00  8.24           C
ATOM    516  C   ARG A  92      32.082  61.069  25.480  1.00  7.94           C
ATOM    517  O   ARG A  92      32.550  59.985  25.866  1.00 10.38           O
ATOM    518  CB AARG A  92      32.838  63.133  26.771  0.50  8.24           C
ATOM    519  CB BARG A  92      32.507  63.088  26.731  0.50  8.15           C
ATOM    520  CG AARG A  92      33.908  64.231  26.908  0.50  8.00           C
ATOM    521  CG BARG A  92      33.428  64.191  27.114  0.50  7.60           C
ATOM    522  CD AARG A  92      35.345  63.768  27.104  0.50  8.06           C
ATOM    523  CD BARG A  92      34.674  63.668  27.784  0.50  8.27           C
ATOM    524  NE AARG A  92      35.471  63.103  28.390  0.50  7.84           N
ATOM    525  NE BARG A  92      34.344  63.138  29.101  0.50  7.94           N
ATOM    526  CZ AARG A  92      36.586  62.621  28.904  0.50  9.02           C
ATOM    527  CZ BARG A  92      35.190  62.540  29.954  0.50  8.28           C
ATOM    528  NH1AARG A  92      37.706  62.705  28.233  0.50 14.66           N
ATOM    529  NH1BARG A  92      36.484  62.428  29.633  0.50  7.37           N
ATOM    530  NH2AARG A  92      36.555  62.048  30.101  0.50 11.64           N
ATOM    531  NH2BARG A  92      34.767  62.156  31.175  0.50  6.81           N
ATOM    532  N   ASP A  93      30.815  61.167  25.066  1.00  7.08           N
ATOM    533  CA  ASP A  93      29.930  59.993  25.026  1.00  7.52           C
ATOM    534  C   ASP A  93      30.579  58.852  24.237  1.00  6.59           C
ATOM    535  O   ASP A  93      30.664  57.719  24.694  1.00  8.25           O
ATOM    536  CB  ASP A  93      28.578  60.370  24.419  1.00  9.20           C
ATOM    537  CG  ASP A  93      27.912  61.536  25.111  1.00 11.51           C
ATOM    538  OD1 ASP A  93      27.283  61.299  26.128  1.00 18.61           O
ATOM    539  OD2 ASP A  93      28.068  62.711  24.724  1.00 16.53           O
TER
HETATM 1512 FE   FE  A1210      26.685  63.748  26.686  1.00  8.65          FE
HETATM 1665  O   HOH A2153      30.815  65.185  24.984  1.00 16.78           O
HETATM 1666  O   HOH A2154      26.856  65.362  24.718  1.00 13.80           O
HETATM 1674  O   HOH A2162      28.848  59.920  28.597  1.00 33.42           O
END
"""


def exercise_2c2u_symmetry_materialization():
  """Lock in the symmetry-expanded shape of the 2C2U Fe region.

  Asserts:

  * exactly one Fe in the materialized model (the three sym_ops put the
    metal at the same point; only one survives deduplication);
  * three Asp 93 atom groups, one per coordinating sym_op;
  * three HOH 2154 oxygens, one per coordinating sym_op;
  * chain IDs are all single character (PDB-compatible) and include
    the parent "A" plus at least two additional chains for the
    symmetry images.
  """
  result = _run_endo_exo_on_string(_2C2U_FE_SPHERE_PDB)
  hier = result["model"].get_hierarchy()

  fe_atoms = 0
  asp93_groups = 0
  hoh2154_groups = 0
  chain_ids = set()
  for ch in hier.chains():
    chain_ids.add(ch.id.strip())
    assert len(ch.id.strip()) == 1, (
      f"chain id {ch.id!r} is not a single character; this would "
      f"overflow the PDB chain field and break downstream restraint "
      f"reconstruction")
    for rg in ch.residue_groups():
      for ag in rg.atom_groups():
        rn = ag.resname.strip().upper()
        if rn == "FE":
          fe_atoms += len(list(ag.atoms()))
        elif rn == "ASP" and rg.resseq.strip() == "93":
          asp93_groups += 1
        elif rn == "HOH" and rg.resseq.strip() == "2154":
          hoh2154_groups += 1

  assert fe_atoms == 1, (
    f"expected 1 Fe atom after special-position dedup; got {fe_atoms}")
  assert asp93_groups == 3, (
    f"expected 3 Asp 93 atom groups (one per sym_op); got {asp93_groups}")
  assert hoh2154_groups == 3, (
    f"expected 3 HOH 2154 atom groups (one per sym_op); got "
    f"{hoh2154_groups}")

  assert "A" in chain_ids, (
    f"parent chain 'A' missing from materialized region: {chain_ids}")
  assert len(chain_ids) >= 3, (
    f"expected >=3 chain ids (identity + two sym images); got "
    f"{sorted(chain_ids)}")


def exercise_2c2u_fe_coordination_distances():
  """The three Asp 93 OD1/OD2 and three HOH 2154 oxygens that survive
  materialization should all be within Fe coordination range
  (<= 3.0 A from the deduplicated Fe)."""
  result = _run_endo_exo_on_string(_2C2U_FE_SPHERE_PDB)
  hier = result["model"].get_hierarchy()
  atoms = list(hier.atoms())

  seed_iseq = result["seed_iseqs"][0]
  fe_xyz = matrix.col(atoms[seed_iseq].xyz)
  assert atoms[seed_iseq].element.strip().upper() == "FE"

  asp_o_distances = []
  hoh_o_distances = []
  for a in atoms:
    rn = a.parent().resname.strip().upper()
    rg = a.parent().parent()
    name = a.name.strip().upper()
    if rn == "ASP" and rg.resseq.strip() == "93" and name in ("OD1", "OD2"):
      asp_o_distances.append((matrix.col(a.xyz) - fe_xyz).length())
    elif rn == "HOH" and rg.resseq.strip() == "2154":
      hoh_o_distances.append((matrix.col(a.xyz) - fe_xyz).length())

  # 3 Asp93 * (OD1 + OD2) = 6 carboxylate oxygens; at least 3 of them
  # (one per Asp93 image) coordinate the Fe.
  close_asp_o = [d for d in asp_o_distances if d < 3.0]
  assert len(close_asp_o) >= 3, (
    f"expected >=3 Asp 93 oxygens within 3 A of Fe; got "
    f"{sorted(asp_o_distances)}")

  assert len(hoh_o_distances) == 3
  for d in hoh_o_distances:
    assert d < 3.0, (
      f"HOH 2154 O-Fe distance {d:.2f} A outside coordination range")


def exercise_2c2u_symmetry_truncation_consistency():
  """At radius=6 the three symmetry copies of Asp 93 must be truncated
  identically.

  Asp 93's CA (5.2 A) and CB (4.5 A) both sit inside a 6 A sphere around
  the Fe.  The radius search is symmetry-aware, so every copy's CA/CB are
  seeded and protected -- none is cut at its preferred CA-CB site.  Before
  the symmetry-aware seeding the identity copy kept CA/CB (its atoms were
  ASU seeds) while the symmetry images, reached only by BFS, were cut at
  CA-CB; this pins that asymmetry shut."""
  result = _run_endo_exo_on_string(_2C2U_FE_SPHERE_PDB, radius=6.0)
  hier = result["model"].get_hierarchy()

  asp_copies = []
  for ch in hier.chains():
    for rg in ch.residue_groups():
      for ag in rg.atom_groups():
        if ag.resname.strip().upper() == "ASP" and rg.resseq.strip() == "93":
          # real (non-cap) heavy-atom names; a capped atom carries element H
          real_heavy = {a.name.strip() for a in ag.atoms()
                        if a.element.strip().upper() != "H"}
          asp_copies.append((ch.id.strip(), real_heavy))

  assert len(asp_copies) == 3, (
    f"expected 3 Asp 93 copies (one per 3-fold sym_op); got "
    f"{len(asp_copies)}")
  for chain_id, real_heavy in asp_copies:
    assert {"CA", "CB"} <= real_heavy, (
      f"Asp 93 copy in chain {chain_id} is missing real CA/CB (cut at "
      f"CA-CB instead of keeping the in-radius backbone): "
      f"{sorted(real_heavy)}")
  # All three copies must keep the *same* heavy-atom set.
  atom_sets = {frozenset(real_heavy) for _id, real_heavy in asp_copies}
  assert len(atom_sets) == 1, (
    f"Asp 93 symmetry copies truncated differently: "
    f"{[(cid, sorted(s)) for cid, s in asp_copies]}")


def exercise_residue_composition():
  """The default buffer (radius=5, depth=3) pulls in a stable scaffold
  around the four coordinating Cys.  Pin the residue type counts."""
  result = _run_endo_exo_on_string(_1BQ8_FE_SPHERE_PDB)
  hier = result["model"].get_hierarchy()

  resname_counts = {}
  for ag in hier.atom_groups():
    rn = ag.resname.strip().upper()
    resname_counts[rn] = resname_counts.get(rn, 0) + 1

  expected = {
    "ALA": 1, "CYS": 4, "FE": 1, "GLY": 2,
    "ILE": 2, "TYR": 1, "VAL": 2,
  }
  assert resname_counts == expected, (
    f"residue composition drifted:\n"
    f"  expected: {expected}\n"
    f"  got     : {resname_counts}")


def exercise_selection_seed_terminates():
  """A non-metal ``selection`` seed on a symmetric crystal grows a bounded,
  covalent-only region containing the seeded residue with its side chain
  kept whole."""
  result = _run_endo_exo_on_string(model_1yjp, selection="resid 3")

  # Finite region, not an unbounded lattice walk.
  n = result["model"].get_number_of_atoms()
  assert 0 < n < 200, f"resid 3 region atom count looks unbounded: {n}"

  # The seeded ASN 3 is present with its side chain kept whole.
  asn3 = _residue_atom_names(result, "3")
  assert {"CB", "CG", "OD1", "ND2"} <= asn3, (
    f"seeded ASN 3 side chain incomplete: {sorted(asn3)}")


# ===========================================================================
# Engine unit tests
#
# These drive the pure engine pieces directly (no metal scan, no full
# Program pipeline) and cover the subtle logic -- symmetry-op
# canonicalisation, cap placement, bond-cut heuristics, and adjacency
# construction -- that the count-based integration tests above do not.
# ===========================================================================

def _atoms_by_name(pdb_str):
  """Parse *pdb_str* and return ``{atom_name: atom}`` with i_seqs reset to
  positional order."""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  hier = pdb_in.construct_hierarchy()
  atoms = hier.atoms()
  atoms.reset_i_seq()
  return {a.name.strip(): a for a in atoms}


def exercise_canon_op():
  """``_canon_op`` must give a canonical, hash-stable representative: two
  rt_mx values that compare equal (same xyz) must hash equal and satisfy
  set membership after canonicalisation, regardless of how they were
  built.  This is the contract the symmetry-aware BFS relies on when it
  keys ``visited`` / ``cap_candidates`` on ``(i_seq, rt_mx)`` nodes."""
  op = sgtbx.rt_mx("x,y,z+1")
  # Same operation reached via composition with the identity -- this can
  # leave a different internal denominator / representation.
  op_composed = op.multiply(sgtbx.rt_mx())
  assert op.as_xyz() == op_composed.as_xyz()

  c1 = _canon_op(op)
  c2 = _canon_op(op_composed)

  # canonical form preserves the operation ...
  assert c1.as_xyz() == op.as_xyz()
  # ... is idempotent ...
  assert _canon_op(c1).as_xyz() == c1.as_xyz()
  # ... and is hash-stable: equal ops hash equal and live as one set member.
  assert c1 == c2
  assert hash(c1) == hash(c2)
  assert c2 in {c1}

  # The inverse of a pure translation negates it; canonicalisation must
  # survive that too (the reverse-edge case in build_adjacency).
  inv = _canon_op(op.inverse())
  assert inv.as_xyz() == sgtbx.rt_mx("x,y,z-1").as_xyz()


_CAP_PDB = """\
ATOM      1  C1  XXX A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C2  XXX A   1       1.500   0.000   0.000  1.00  0.00           C
"""


def exercise_hydrogen_capper():
  """``cap_atom`` retypes the cap to hydrogen and moves it to 1.1 A along
  the anchor->cap direction; ``None`` arguments are a no-op."""
  capper = HydrogenCapper(log=io.StringIO())

  by_name = _atoms_by_name(_CAP_PDB)
  anchor, cap = by_name["C1"], by_name["C2"]
  capper.cap_atom(anchor, cap)

  assert cap.element.strip().upper() == "H"
  d = (matrix.col(cap.xyz) - matrix.col(anchor.xyz)).length()
  assert approx_equal(d, 1.1, eps=1e-6), d
  # cap lies on the original +x direction from the anchor
  assert approx_equal(cap.xyz, (1.1, 0.0, 0.0), eps=1e-6), cap.xyz

  # None on either side is a no-op (must not raise / must not mutate).
  fresh = _atoms_by_name(_CAP_PDB)
  before = fresh["C2"].xyz
  capper.cap_atom(None, fresh["C2"])
  capper.cap_atom(fresh["C1"], None)
  assert fresh["C2"].xyz == before
  assert fresh["C2"].element.strip().upper() == "C"


# Lysine fragment: backbone N-CA-C plus the CD-CE sidechain bond that is a
# PREFERRED_CUTS site for LYS.  Geometry is only nominal; the preferred /
# backbone checks are name + adjacency based.
_LYS_PDB = """\
ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  LYS A   1       1.450   0.000   0.000  1.00  0.00           C
ATOM      3  C   LYS A   1       2.000   1.400   0.000  1.00  0.00           C
ATOM      4  CD  LYS A   1       3.000   0.000   0.000  1.00  0.00           C
ATOM      5  CE  LYS A   1       4.500   0.000   0.000  1.00  0.00           C
"""

# Two sp3 carbons 1.54 A apart for the geometric C-C heuristic branch.
_CC_PDB = """\
ATOM      1  CA  XXX A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CB  XXX A   1       1.540   0.000   0.000  1.00  0.00           C
"""


def _edge(adj, i, j):
  """Add an undirected bare edge (the sym_op is irrelevant to the cut
  checks, which drop it via _neighbour_iseqs)."""
  adj.setdefault(i, set()).add((j, None))
  adj.setdefault(j, set()).add((i, None))


def exercise_bond_cut_preferred():
  """With ``use_preferred_cuts=True`` the PREFERRED_CUTS table decides:
  LYS CD-CE is a cut site, LYS CA-CD is not."""
  by_name = _atoms_by_name(_LYS_PDB)
  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())

  adj = {}
  _edge(adj, by_name["CD"].i_seq, by_name["CE"].i_seq)
  _edge(adj, by_name["CA"].i_seq, by_name["CD"].i_seq)

  assert det.is_cc_single_sp3_bond(
    "LYS", by_name["CD"], by_name["CE"], adj) is True
  # CA is not in LYS's preferred {CD, CE} set.
  assert det.is_cc_single_sp3_bond(
    "LYS", by_name["CA"], by_name["CD"], adj) is False


def exercise_bond_cut_heuristic():
  """For a residue absent from PREFERRED_CUTS the geometric heuristic
  applies: two degree-4 sp3 carbons 1.42-1.68 A apart are a cut site;
  shifting the distance out of range disqualifies the bond."""
  by_name = _atoms_by_name(_CC_PDB)
  ca, cb = by_name["CA"], by_name["CB"]
  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())

  # Degree-4 carbons: the real C-C bond plus three dummy neighbours each.
  adj = {}
  _edge(adj, ca.i_seq, cb.i_seq)
  for d in (101, 102, 103):
    _edge(adj, ca.i_seq, d)
  for d in (201, 202, 203):
    _edge(adj, cb.i_seq, d)

  # 'XXX' is not in PREFERRED_CUTS -> heuristic branch; 1.54 A is in range.
  assert det.is_cc_single_sp3_bond("XXX", ca, cb, adj) is True

  # Move CB out to 2.0 A: now outside the 1.42-1.68 A window.
  cb.set_xyz((2.0, 0.0, 0.0))
  assert det.is_cc_single_sp3_bond("XXX", ca, cb, adj) is False


def exercise_bond_cut_backbone():
  """``is_ca_c_bond`` / ``is_ca_n_bond`` fire only on a genuine, bonded
  CA->C and CA->N pair (direction and adjacency both matter)."""
  by_name = _atoms_by_name(_LYS_PDB)
  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())

  adj = {}
  _edge(adj, by_name["CA"].i_seq, by_name["C"].i_seq)
  _edge(adj, by_name["CA"].i_seq, by_name["N"].i_seq)

  assert det.is_ca_c_bond(by_name["CA"], by_name["C"], adj) is True
  assert det.is_ca_n_bond(by_name["CA"], by_name["N"], adj) is True
  # Wrong direction / wrong names.
  assert det.is_ca_c_bond(by_name["C"], by_name["CA"], adj) is False
  assert det.is_ca_n_bond(by_name["CA"], by_name["C"], adj) is False
  # Right names but not adjacent -> not a bond.
  assert det.is_ca_c_bond(by_name["CA"], by_name["C"], {}) is False


# Real 1RYO chain-A LYS 296 (heavy atoms + the CD/CG hydrogens).  CD and
# NZ are 4.8 / 4.3 A from the Fe in the full structure; CE is 5.0 A out.
# The CD/CG hydrogens give those carbons covalent degree 4 so the
# geometric C-C heuristic (used by the fallback) accepts the CD-CG bond;
# CD-CG is 1.52 A apart, inside the 1.42-1.68 A sp3 window.
_LYS_296_PDB = """\
CRYST1   44.092   57.252  135.988  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   LYS A 296      46.335  42.452  48.863  1.00  4.54           N
ATOM      2  CA  LYS A 296      47.390  42.927  47.971  1.00  4.20           C
ATOM      3  C   LYS A 296      48.088  41.772  47.253  1.00  4.42           C
ATOM      4  O   LYS A 296      48.337  40.715  47.847  1.00  4.42           O
ATOM      5  CB  LYS A 296      48.426  43.737  48.755  1.00  5.02           C
ATOM      6  CG  LYS A 296      47.999  45.167  49.009  1.00  4.18           C
ATOM      7  CD  LYS A 296      48.823  45.848  50.096  1.00  4.88           C
ATOM      8  CE  LYS A 296      50.305  45.946  49.760  1.00  5.79           C
ATOM      9  NZ  LYS A 296      50.995  46.848  50.737  1.00  7.00           N
ATOM     10  HD2 LYS A 296      48.728  45.278  51.020  1.00  4.88           H
ATOM     11  HD3 LYS A 296      48.447  46.860  50.242  1.00  4.88           H
ATOM     12  HG2 LYS A 296      48.116  45.741  48.090  1.00  4.18           H
ATOM     13  HG3 LYS A 296      46.955  45.175  49.324  1.00  4.18           H
"""


def exercise_preferred_cut_fallback():
  """Reproduce the 1RYO Fe / LYS A 296 case: the radius search seeds CD
  and NZ (4.8 / 4.3 A from the Fe) while CE sits just outside (5.0 A), so
  the LYS preferred cut CD-CE ends up with both endpoints interior and can
  no longer be made.

  Baseline ``grow_by_depth`` then walks inward past CD into the backbone
  (CG, CB, CA, N pulled in; C capped).  ``grow_region`` with
  ``preferred_cut_fallback=True`` detects the consumed cut and re-cuts at
  CD-CG with the geometric heuristic, trimming the region to just the
  {CD, CE, NZ} tip with CG as the cap."""
  pdb_in = iotbx.pdb.input(
    source_info=None, lines=_LYS_296_PDB.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  atoms = model.get_hierarchy().atoms()
  atoms.reset_i_seq()
  iseq = {a.name.strip(): a.i_seq for a in atoms}
  name_of = {a.i_seq: a.name.strip() for a in atoms}
  elem_of = {a.i_seq: a.element.strip().upper() for a in atoms}

  # Tagged adjacency must carry a real rt_mx on each edge (the BFS
  # composes ops); the identity stands in for an intra-ASU bond.
  identity = _canon_op(sgtbx.rt_mx())
  adjacency = defaultdict(set)
  def _bond(a, b):
    adjacency[iseq[a]].add((iseq[b], identity))
    adjacency[iseq[b]].add((iseq[a], identity))
  for a, b in [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"),
               ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"),
               ("CD", "HD2"), ("CD", "HD3"), ("CG", "HG2"), ("CG", "HG3")]:
    _bond(a, b)

  # CD and NZ are the atoms the radius search would seed.
  seeds = {iseq["CD"], iseq["NZ"]}

  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())
  grower = QMRegionGrower(det, log=io.StringIO())

  def _split(visited, caps):
    cap_iseqs = {c[0] for c in caps}
    interior = {name_of[i] for (i, _op) in visited
                if i not in cap_iseqs and elem_of[i] not in ("H", "D")}
    return interior, {name_of[i] for i in cap_iseqs}

  # Baseline: preferred cut consumed -> overgrowth into the backbone.
  v0, c0 = grower.grow_by_depth(seeds, adjacency, model)
  interior0, _caps0 = _split(v0, c0)
  assert {"N", "CA", "CB", "CG"} <= interior0, (
    f"baseline should pull the backbone in; interior={sorted(interior0)}")

  # Fallback: re-cut at CD-CG, trimming everything inward of it.
  v1, c1 = grower.grow_region(
    seeds, adjacency, model, preferred_cut_fallback=True)
  interior1, caps1 = _split(v1, c1)
  assert interior1 == {"CD", "CE", "NZ"}, (
    f"fallback interior drifted: expected {{CD, CE, NZ}}, got "
    f"{sorted(interior1)}")
  assert caps1 == {"CG"}, (
    f"fallback cap drifted: expected {{CG}}, got {sorted(caps1)}")


# Real 1RYO chain-A LEU 62 - ASP 63 dipeptide (peptide bond LEU62 C - ASP63
# N).  In the full structure ASP 63 coordinates the Fe (it is radius-seeded)
# while LEU 62 is >6 A away -- it only enters the region as backbone
# overgrowth off ASP 63.
_LEU62_ASP63_PDB = """\
CRYST1   44.092   57.252  135.988  90.00  90.00  90.00 P 21 21 21
ATOM    867  N   LEU A  62      42.072  46.636  54.794  1.00  4.17           N
ATOM    868  CA  LEU A  62      42.362  48.055  54.653  1.00  4.15           C
ATOM    869  C   LEU A  62      43.824  48.433  54.764  1.00  4.47           C
ATOM    870  O   LEU A  62      44.579  47.852  55.543  1.00  4.55           O
ATOM    871  CB  LEU A  62      41.613  48.852  55.726  1.00  5.28           C
ATOM    872  CG  LEU A  62      40.088  48.813  55.774  1.00  5.85           C
ATOM    873  CD1 LEU A  62      39.590  49.620  56.961  1.00  6.23           C
ATOM    874  CD2 LEU A  62      39.523  49.367  54.480  1.00  5.86           C
ATOM    875  H   LEU A  62      42.395  46.287  55.502  1.00  4.17           H
ATOM    876  HA  LEU A  62      42.047  48.357  53.787  1.00  4.15           H
ATOM    877  HB2 LEU A  62      41.924  48.540  56.590  1.00  5.28           H
ATOM    878  HB3 LEU A  62      41.862  49.783  55.624  1.00  5.28           H
ATOM    886  N   ASP A  63      44.202  49.436  53.978  1.00  4.55           N
ATOM    887  CA  ASP A  63      45.538  50.008  54.033  1.00  4.24           C
ATOM    888  C   ASP A  63      45.620  50.561  55.466  1.00  4.40           C
ATOM    889  O   ASP A  63      44.602  50.957  56.042  1.00  4.34           O
ATOM    890  CB  ASP A  63      45.644  51.143  53.015  1.00  4.42           C
ATOM    891  CG  ASP A  63      46.806  52.071  53.293  1.00  4.13           C
ATOM    892  OD1 ASP A  63      47.954  51.770  52.885  1.00  4.01           O
ATOM    893  OD2 ASP A  63      46.556  53.112  53.936  1.00  4.67           O
ATOM    895  HA  ASP A  63      46.214  49.346  53.894  1.00  4.24           H
ATOM    896  HB2 ASP A  63      45.767  50.764  52.131  1.00  4.42           H
ATOM    897  HB3 ASP A  63      44.829  51.668  53.041  1.00  4.42           H
"""


def exercise_overgrowth_geometric_cut():
  """A residue with no seed (radius) atom is pure backbone overgrowth and,
  under preferred_cut_fallback, defers to the geometric C-C heuristic.

  ASP 63 is the radius-seeded ligand; LEU 62 enters only via the peptide
  bond.  Baseline (preferred cuts) trims LEU at its preferred CB-CG site,
  keeping CB.  With the overgrowth rule LEU is cut at the first sp3 C-C
  bond instead -- CA-CB -- so CB becomes the cap and the rest of the
  sidechain (CG, CD1, CD2) is dropped."""
  pdb_in = iotbx.pdb.input(
    source_info=None, lines=_LEU62_ASP63_PDB.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  atoms = model.get_hierarchy().atoms()
  atoms.reset_i_seq()

  by_res = {}
  for a in atoms:
    resseq = a.parent().parent().resseq.strip()
    by_res.setdefault(resseq, {})[a.name.strip()] = a.i_seq
  leu, asp = by_res["62"], by_res["63"]
  name_of = {a.i_seq: a.name.strip() for a in atoms}
  resseq_of = {a.i_seq: a.parent().parent().resseq.strip() for a in atoms}
  elem_of = {a.i_seq: a.element.strip().upper() for a in atoms}

  identity = _canon_op(sgtbx.rt_mx())
  adjacency = defaultdict(set)
  def _bond(i, j):
    adjacency[i].add((j, identity))
    adjacency[j].add((i, identity))
  for a, b in [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"),
               ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("N", "H"),
               ("CA", "HA"), ("CB", "HB2"), ("CB", "HB3")]:
    _bond(leu[a], leu[b])
  for a, b in [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"),
               ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2"),
               ("CA", "HA"), ("CB", "HB2"), ("CB", "HB3")]:
    _bond(asp[a], asp[b])
  _bond(leu["C"], asp["N"])  # peptide bond

  # ASP 63 is the radius-seeded ligand; LEU 62 carries no seed atom.
  seeds = set(asp.values())

  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())
  grower = QMRegionGrower(det, log=io.StringIO())

  def _leu_atoms(visited, caps):
    cap_iseqs = {c[0] for c in caps}
    present = {name_of[i] for (i, _o) in visited if resseq_of[i] == "62"}
    caps_leu = {name_of[i] for i in cap_iseqs if resseq_of[i] == "62"}
    interior = {n for n in present
                if n not in caps_leu and elem_of[leu[n]] not in ("H", "D")}
    return interior, caps_leu, present

  # Baseline: LEU keeps CB and caps CG (its preferred CB-CG cut).
  v0, c0 = grower.grow_by_depth(seeds, adjacency, model)
  interior0, caps0, _present0 = _leu_atoms(v0, c0)
  assert "CB" in interior0 and "CG" in caps0, (
    f"baseline LEU should cut at CB-CG; interior={sorted(interior0)} "
    f"caps={sorted(caps0)}")

  # Overgrowth rule: LEU defers to the geometric heuristic and cuts CA-CB.
  v1, c1 = grower.grow_region(
    seeds, adjacency, model, preferred_cut_fallback=True)
  interior1, caps1, present1 = _leu_atoms(v1, c1)
  assert "CB" in caps1, (
    f"overgrowth LEU should cap CB (CA-CB cut); caps={sorted(caps1)}")
  assert not ({"CG", "CD1", "CD2"} & present1), (
    f"overgrowth LEU sidechain beyond CB should be dropped; "
    f"present={sorted(present1)}")


def _simple_bond_proxies(pdb_str):
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  model.process(
    pdb_interpretation_params=model.get_current_pdb_interpretation_params(),
    make_restraints=True)
  grm = model.get_restraints_manager().geometry
  simple, _asu = grm.get_all_bond_proxies(sites_cart=model.get_sites_cart())
  return simple


def _asu_bond_proxies(pdb_str, distance_cutoff=3.2):
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  pair_asu_table = model.get_xray_structure().pair_asu_table(
    distance_cutoff=distance_cutoff)
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    pair_asu_table=pair_asu_table)
  return list(sorted_asu_proxies.asu), pair_asu_table.asu_mappings()


def exercise_build_adjacency():
  """Intra-ASU bonds carry the identity op on both directed edges;
  symmetry-crossing bonds carry the rt_mx forward and its inverse on the
  reverse edge."""
  builder = AtomGraphBuilder()
  identity = _canon_op(sgtbx.rt_mx())

  simple = _simple_bond_proxies(_1BQ8_FE_SPHERE_PDB)
  assert len(simple) > 0
  adj = builder.build_adjacency(simple, [], None)
  for proxy in simple:
    i_seq, j_seq = proxy.i_seqs
    assert (j_seq, identity) in adj[i_seq]
    assert (i_seq, identity) in adj[j_seq]

  asu, asu_mappings = _asu_bond_proxies(_2C2U_FE_SPHERE_PDB)
  assert len(asu) > 0
  adj = builder.build_adjacency([], asu, asu_mappings)
  n_non_identity = 0
  for proxy in asu:
    op = _canon_op(asu_mappings.get_rt_mx_ji(proxy))
    if op != identity:
      n_non_identity += 1
    assert (proxy.j_seq, op) in adj[proxy.i_seq]
    assert (proxy.i_seq, _canon_op(op.inverse())) in adj[proxy.j_seq]
  assert n_non_identity > 0


def run():
  exercise_submodel_shape()
  exercise_cys_coordination()
  exercise_residue_composition()
  exercise_residues_to_include()
  exercise_2c2u_symmetry_materialization()
  exercise_2c2u_fe_coordination_distances()
  exercise_2c2u_symmetry_truncation_consistency()
  exercise_selection_seed_terminates()
  # engine unit tests
  exercise_canon_op()
  exercise_hydrogen_capper()
  exercise_bond_cut_preferred()
  exercise_bond_cut_heuristic()
  exercise_bond_cut_backbone()
  exercise_preferred_cut_fallback()
  exercise_overgrowth_geometric_cut()
  exercise_build_adjacency()
  print(format_cpu_times())
  print("OK")


if __name__ == "__main__":
  run()
