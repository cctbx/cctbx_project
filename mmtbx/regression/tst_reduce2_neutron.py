"""reduce2 with neutron distances: every step uses the same setting.

use_neutron_distances=True (and qmready=True, which switches it on) must give
H at the neutron restraint targets, water O-H from the HOH restraints, and an
Optimizer run with neutron radii (main run and flipkins). Bond lengths are
checked against each restraint's own target; bonds whose X-ray and neutron
targets differ by more than 0.04 A must also be far from the wrong target.
"""
from __future__ import absolute_import, division, print_function
import os
from six.moves import cStringIO as StringIO
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from scitbx import matrix
from cctbx import geometry_restraints
from cctbx.geometry_restraints.auto_linking_types import iterate_covalent
from cctbx.geometry_restraints.linking_class import linking_class
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.reduce import Optimizers
from mmtbx.regression.model.tst_model_idealize_h_riding_water import \
  pdb_water, layouts, hoh_cif_096, make_model, check_oh

# 1yjp: Asn/Gln for flipkins, Tyr OH, waters, P21
pdb_1yjp = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
ATOM     21  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     22  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     23  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     24  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.944   1.458   1.213  1.00 10.25           C
ATOM     27  CD  GLN A   4       2.504   2.044  -0.089  1.00 12.43           C
ATOM     28  OE1 GLN A   4       2.744   3.268  -0.190  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.750   1.161  -1.091  1.00  9.05           N
ATOM     30  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     31  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     32  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     33  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     34  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     35  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     36  CD  GLN A   5       1.344   2.414   8.954  1.00 13.10           C
ATOM     37  OE1 GLN A   5       0.774   1.325   9.002  1.00 10.65           O
ATOM     38  NE2 GLN A   5       1.549   3.187  10.039  1.00 12.30           N
ATOM     39  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
HETATM   61  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   62  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   63  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   64  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM   65  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   66  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   67  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
END
"""

# Lys 3 stops at CD: its HD2/HD3 cannot be parameterized, so placement drops
# them with model.select()
pdb_truncated_lys = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      11.204   2.315  10.312  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.555   2.835  10.528  1.00 20.00           C
ATOM      3  C   ALA A   1      13.588   1.726  10.381  1.00 20.00           C
ATOM      4  O   ALA A   1      13.325   0.567  10.721  1.00 20.00           O
ATOM      5  CB  ALA A   1      12.866   3.965   9.556  1.00 20.00           C
ATOM      6  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      7  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      8  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      9  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM     10  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM     11  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM     12  N   LYS A   3      17.560   0.282  11.199  1.00 20.00           N
ATOM     13  CA  LYS A   3      18.643   0.370  12.172  1.00 20.00           C
ATOM     14  C   LYS A   3      19.945  -0.148  11.570  1.00 20.00           C
ATOM     15  O   LYS A   3      19.964  -1.176  10.892  1.00 20.00           O
ATOM     16  CB  LYS A   3      18.290  -0.424  13.432  1.00 20.00           C
ATOM     17  CG  LYS A   3      19.308  -0.301  14.556  1.00 20.00           C
ATOM     18  CD  LYS A   3      18.913  -1.137  15.762  1.00 20.00           C
ATOM     19  OXT LYS A   3      20.995   0.451  11.763  1.00 20.00           O
HETATM   20  O   HOH A   4      24.000   6.000  20.000  1.00 20.00           O
END
"""

INFO_TRUE = "useNeutronDistances = True"
INFO_FALSE = "useNeutronDistances = False"

def run_reduce2(prefix, pdb_str, args=()):
  fn_in, fn_out = prefix + ".pdb", prefix + "_out.pdb"
  fn_txt = prefix + ".txt"
  with open(fn_in, "w") as f:
    f.write(pdb_str)
  for fn in (fn_out, fn_txt):
    if os.path.exists(fn): os.remove(fn)
  log = StringIO()
  run_program(program_class=reduce2.Program, logger=log,
    args=[fn_in, "output.filename=" + fn_out, "overwrite=True",
          "output.description_file_name=" + fn_txt] + list(args))
  with open(fn_txt) as f:
    description = f.read()
  return fn_out, log.getvalue(), description

def interpret(file_name, use_neutron_distances):
  inp = iotbx.pdb.input(file_name)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  p = reduce_hydrogen.get_reduce_pdb_interpretation_params(
    use_neutron_distances)
  model.process(make_restraints=True, pdb_interpretation_params=p)
  return model

def h_bond_lengths(model):
  """(i_seq, j_seq, j_sym) -> (ideal, actual) for covalent bonds with an H/D
  end; asu bonds with their symmetry operator"""
  metal = linking_class().get_origin_id('metal coordination')
  covalent = set(iterate_covalent()) - set([metal])
  sites_cart = model.get_sites_cart()
  hd = model.get_hd_selection()
  bond_proxies = model.get_restraints_manager().geometry.pair_proxies(
    sites_cart=sites_cart).bond_proxies
  result = {}
  for p in bond_proxies.simple:
    i, j = p.i_seqs
    if p.origin_id not in covalent or not (hd[i] or hd[j]): continue
    b = geometry_restraints.bond(sites_cart=sites_cart, proxy=p)
    result[(i, j, 0)] = (p.distance_ideal, b.distance_model)
  asu_mappings = bond_proxies.asu_mappings()
  for p in bond_proxies.asu:
    if p.origin_id not in covalent or not (hd[p.i_seq] or hd[p.j_seq]):
      continue
    b = geometry_restraints.bond(sites_cart=sites_cart,
      asu_mappings=asu_mappings, proxy=p)
    result[(p.i_seq, p.j_seq, p.j_sym)] = (p.distance_ideal, b.distance_model)
  return result

def check_bond_lengths(file_name, use_neutron_distances):
  right = h_bond_lengths(interpret(file_name, use_neutron_distances))
  wrong = h_bond_lengths(interpret(file_name, not use_neutron_distances))
  assert len(right) > 0
  n_discriminating = 0
  for key, (ideal, actual) in right.items():
    assert abs(actual - ideal) < 0.02, (file_name, key, ideal, actual)
    if key not in wrong: continue
    ideal_wrong = wrong[key][0]
    if abs(ideal - ideal_wrong) > 0.04:
      assert abs(actual - ideal_wrong) > 0.02, (file_name, key, ideal_wrong)
      n_discriminating += 1
  assert n_discriminating > 0, file_name
  return n_discriminating

def water_geometry(file_name):
  """per water O and conformer: ([O-H/D lengths], H-O-H angle or None);
  conformer residues combine blank-altloc O with altloc H/D"""
  h = iotbx.pdb.input(file_name).construct_hierarchy()
  seen = set()
  result = []
  for chain in h.only_model().chains():
    for conformer in chain.conformers():
      for residue in conformer.residues():
        if residue.resname not in ("HOH", "DOD"): continue
        ra = residue.atoms()
        o = [a for a in ra if a.element.strip() == "O"]
        hs = [a for a in ra if a.element.strip() in ("H", "D")]
        if len(o) != 1 or not hs: continue
        key = tuple(sorted(a.i_seq for a in [o[0]] + hs))
        if key in seen: continue
        seen.add(key)
        vs = [matrix.col(a.xyz) - matrix.col(o[0].xyz) for a in hs]
        angle = None
        if len(vs) == 2: angle = vs[0].angle(vs[1], deg=True)
        result.append(([v.length() for v in vs], angle))
  return result

def check_water(file_name, d_ideal, n_expected):
  waters = water_geometry(file_name)
  assert len(waters) == n_expected, (file_name, len(waters))
  for lengths, angle in waters:
    for d in lengths:
      assert approx_equal(d, d_ideal, eps=0.005, out=None), (file_name, d)
    if angle is not None:
      assert approx_equal(angle, 103.91, eps=0.1, out=None), (file_name, angle)

def test_xray_default():
  # control: the default run stays at X-ray distances
  fn, log, description = run_reduce2("tst_reduce2_neutron_xray", pdb_1yjp,
    ["exclude_water=False"])
  check_bond_lengths(fn, use_neutron_distances=False)
  check_water(fn, 0.85, 7)
  assert INFO_FALSE in description and INFO_TRUE not in description

def test_neutron_add():
  fn, log, description = run_reduce2("tst_reduce2_neutron_add", pdb_1yjp,
    ["use_neutron_distances=True", "exclude_water=False"])
  check_bond_lengths(fn, use_neutron_distances=True)
  check_water(fn, 0.98, 7)
  assert INFO_TRUE in description and INFO_FALSE not in description

def test_neutron_place():
  fn, log, description = run_reduce2("tst_reduce2_neutron_place", pdb_1yjp,
    ["use_neutron_distances=True", "exclude_water=False", "approach=place"])
  check_bond_lengths(fn, use_neutron_distances=True)
  check_water(fn, 0.98, 7)

def test_neutron_optimize():
  # input water H at 0.85 A (X-ray run), optimized with neutron distances
  fn_in, _, _ = run_reduce2("tst_reduce2_neutron_opt_in", pdb_1yjp,
    ["exclude_water=False"])
  check_water(fn_in, 0.85, 7)
  with open(fn_in) as f:
    pdb_str = f.read()
  fn, log, description = run_reduce2("tst_reduce2_neutron_opt", pdb_str,
    ["use_neutron_distances=True", "exclude_water=False",
     "approach=optimize"])
  check_water(fn, 0.98, 7)
  assert INFO_TRUE in description and INFO_FALSE not in description

def test_neutron_flipkins():
  if not os.path.isdir("tst_reduce2_neutron_flipkins"):
    os.mkdir("tst_reduce2_neutron_flipkins")
  fn, log, description = run_reduce2("tst_reduce2_neutron_fk", pdb_1yjp,
    ["use_neutron_distances=True", "add_flip_movers=True",
     "output.flipkin_directory=tst_reduce2_neutron_flipkins"])
  assert INFO_TRUE in description and INFO_FALSE not in description
  # each flipkin configuration logs its own Optimizer run
  assert log.count(INFO_TRUE) >= 2, log.count(INFO_TRUE)
  assert INFO_FALSE not in log

def test_qmready():
  fn, log, description = run_reduce2("tst_reduce2_neutron_qmready", pdb_1yjp,
    ["qmready=True"])
  assert "qmready=True: use_neutron_distances set to True" in log
  check_bond_lengths(fn, use_neutron_distances=True)
  check_water(fn, 0.98, 7)
  assert INFO_TRUE in description and INFO_FALSE not in description

def test_neutron_after_select():
  # placement drops the unparameterizable Lys HD2/HD3 with model.select();
  # the water must still get the neutron target afterwards
  fn, log, description = run_reduce2("tst_reduce2_neutron_select",
    pdb_truncated_lys, ["use_neutron_distances=True", "exclude_water=False"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  names = [a.name.strip() for a in h.atoms() if a.parent().resname == "LYS"]
  assert "HD2" not in names and "HD3" not in names, names
  check_water(fn, 0.98, 1)

def test_water_layouts():
  # 9l: waters the atom_group rule skips, from 0.90 A inputs; per O and
  # conformer: (number of H/D, expected number of groups)
  cryst1 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1\n"
  for key, n_groups in (("a", 1), ("b", 1), ("c", 2), ("d", 2), ("e", 1)):
    pdb_str = cryst1 + layouts[key] + "END\n"
    prefix = "tst_reduce2_neutron_layout_" + key
    with open(prefix + "_in.pdb", "w") as f:
      f.write(pdb_str)
    for lengths, angle in water_geometry(prefix + "_in.pdb"):
      for d in lengths: assert approx_equal(d, 0.90, eps=0.001, out=None)
    fn, log, description = run_reduce2(prefix, pdb_str,
      ["approach=optimize", "use_neutron_distances=True"])
    check_water(fn, 0.98, n_groups)
    if key == "c":
      # A and B keep distinct H positions
      h = iotbx.pdb.input(fn).construct_hierarchy()
      xyz = dict(((a.name.strip(), a.parent().altloc), matrix.col(a.xyz))
                 for a in h.atoms())
      for name in ("H1", "H2"):
        assert (xyz[(name, "A")] - xyz[(name, "B")]).length() > 0.5, name

def test_optimizer_passes_flag():
  # the Optimizer's own setup_riding_h_manager() call idealizes the water
  for flag, d in ((True, 0.96), (False, 0.85)):
    model = make_model(pdb_water, True, "n_gaussian", cif=hoh_cif_096)
    Optimizers.Optimizer(Optimizers._philLike(), False, model,
      useNeutronDistances=flag, verbosity=0)
    check_oh(model, d)
  # default useNeutronDistances: scattering-table rule as before
  model = make_model(pdb_water, False, "electron")
  Optimizers.Optimizer(Optimizers._philLike(), False, model, verbosity=0)
  check_oh(model, 0.98)

if __name__ == "__main__":
  test_xray_default()
  test_neutron_add()
  test_neutron_place()
  test_neutron_optimize()
  test_neutron_flipkins()
  test_qmready()
  test_neutron_after_select()
  test_water_layouts()
  test_optimizer_passes_flag()
  print("OK")
