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
from libtbx.utils import null_out, Sorry
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

FLOATING = "H/D without a covalent partner, removed before optimization:"

def h_input_without_tyr_oh(prefix, args=()):
  """1yjp with H from reduce2, then Tyr OH deleted: its HH has no partner"""
  fn, _, _ = run_reduce2(prefix + "_h", pdb_1yjp, list(args))
  with open(fn) as f:
    lines = [l for l in f.read().splitlines() if " OH  TYR" not in l]
  return "\n".join(lines) + "\n"

def atom_labels(file_name):
  h = iotbx.pdb.input(file_name).construct_hierarchy()
  return sorted((a.parent().parent().resseq, a.parent().resname, a.name,
                 a.parent().altloc) for a in h.atoms())

def test_floating_h_removed():
  # 5a (none): approach=optimize on input with an H whose parent heavy atom
  # is missing: the run completes, the H is gone and listed, all other atoms
  # unchanged in count and names
  pdb_str = h_input_without_tyr_oh("tst_reduce2_neutron_float")
  with open("tst_reduce2_neutron_float_in.pdb", "w") as f:
    f.write(pdb_str)
  fn, log, description = run_reduce2("tst_reduce2_neutron_float", pdb_str,
    ["approach=optimize"])
  assert FLOATING in description and FLOATING in log
  listed = description.split(FLOATING)[1].strip().splitlines()[0]
  assert listed.split() == ["HH", "TYR", "A", "7"], listed
  before = [l for l in atom_labels("tst_reduce2_neutron_float_in.pdb")
            if l[1:3] != ("TYR", " HH ")]
  assert atom_labels(fn) == before

def test_floating_h_water_098():
  # 9k(b): water at 0.85 A in the input, a floating H removed with select();
  # with neutron distances the written water O-H is 0.98 A
  pdb_str = h_input_without_tyr_oh("tst_reduce2_neutron_float_w",
    ["exclude_water=False"])
  fn, log, description = run_reduce2("tst_reduce2_neutron_float_w", pdb_str,
    ["approach=optimize", "use_neutron_distances=True",
     "exclude_water=False"])
  assert FLOATING in description
  assert not [l for l in atom_labels(fn) if l[1:3] == ("TYR", " HH ")]
  check_water(fn, 0.98, 7)

def test_no_floating_h_unchanged():
  # input without floating H/D: no removal, output as without the step
  fn_in, _, _ = run_reduce2("tst_reduce2_neutron_nofloat_h", pdb_1yjp, [])
  with open(fn_in) as f:
    pdb_str = f.read()
  fn, log, description = run_reduce2("tst_reduce2_neutron_nofloat", pdb_str,
    ["approach=optimize"])
  assert FLOATING not in description
  original = reduce2.Program._RemoveFloatingHD
  reduce2.Program._RemoveFloatingHD = lambda self: ("", [])
  try:
    fn_ref, _, _ = run_reduce2("tst_reduce2_neutron_nofloat_ref", pdb_str,
      ["approach=optimize"])
  finally:
    reduce2.Program._RemoveFloatingHD = original
  def atom_lines(file_name):
    with open(file_name) as f:
      return [l for l in f if l.startswith(("ATOM", "HETATM", "ANISOU"))]
  assert atom_lines(fn) == atom_lines(fn_ref)

def test_floating_h_flipkins():
  # orphan water H (water 8 without its O): removed in the main run and not
  # brought back by the flipkin passes, which restart from the initial model
  fn, _, _ = run_reduce2("tst_reduce2_neutron_orphan_h", pdb_1yjp,
    ["exclude_water=False"])
  with open(fn) as f:
    pdb_str = "\n".join(l for l in f.read().splitlines()
                        if " O   HOH A   8" not in l) + "\n"
  assert pdb_str.count("HOH A   8") == 2       # H1, H2
  def orphans(hierarchy):
    return [a for a in hierarchy.atoms() if a.parent().resname == "HOH" and
            a.parent().parent().resseq.strip() == "8"]
  seen = []
  original = Optimizers.Optimizer.__init__
  def init(self, probePhil, addFlipMovers, model, *args, **kwargs):
    seen.append(len(orphans(model.get_hierarchy())))
    return original(self, probePhil, addFlipMovers, model, *args, **kwargs)
  if not os.path.isdir("tst_reduce2_neutron_orphan_fk"):
    os.mkdir("tst_reduce2_neutron_orphan_fk")
  for keep in ("True", "False"):
    del seen[:]
    Optimizers.Optimizer.__init__ = init
    try:
      fn, log, description = run_reduce2(
        "tst_reduce2_neutron_orphan_" + keep, pdb_str,
        ["approach=optimize", "keep_existing_H=" + keep, "exclude_water=False",
         "add_flip_movers=True",
         "output.flipkin_directory=tst_reduce2_neutron_orphan_fk"])
    finally:
      Optimizers.Optimizer.__init__ = original
    # main run and both flipkin configurations, none with the orphan H
    assert seen == [0, 0, 0], (keep, seen)
    assert description.count(FLOATING) == 1
    listed = description.split(FLOATING)[1].strip().splitlines()[:2]
    assert sorted(l.split()[0] for l in listed) == ["H1", "H2"], listed
    assert all(l.split()[1:] == ["HOH", "A", "8"] for l in listed), listed
    assert orphans(iotbx.pdb.input(fn).construct_hierarchy()) == []

def test_water_without_o_add():
  # approach=add keep_existing_H=True on 1yjp with water H and water 8's O
  # deleted: no crash in workaround_002(); the two H are removed and reported
  # as missing their heavy atom, not as not parameterized
  fn, _, _ = run_reduce2("tst_reduce2_neutron_no_o_h", pdb_1yjp,
    ["exclude_water=False"])
  with open(fn) as f:
    pdb_str = "\n".join(l for l in f.read().splitlines()
                        if " O   HOH A   8" not in l) + "\n"
  fn, log, description = run_reduce2("tst_reduce2_neutron_no_o", pdb_str,
    ["approach=add", "keep_existing_H=True", "exclude_water=False"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  assert not [a for a in h.atoms() if a.parent().resname == "HOH" and
              a.parent().parent().resseq.strip() == "8"]
  report = log.split("were not placed because a neighbouring heavy atom is")
  assert len(report) == 2 and "HOH A 8" in report[1].split("\n\n")[0], log
  for block in log.split("could not be parameterized")[1:]:
    assert "HOH A   8" not in block.split("\n\n")[0]

# ------------------------------------------------------------------------------
# step 3.1: phil, validate(), file names, _RebuildModel()

OPTIONS = ("none", "all_h", "all_d", "hd_and_h", "hd_and_d", "all_hd")

def make_program(prefix, pdb_str, cif_str=None, **params):
  """reduce2.Program as run_program builds it, validate() not yet called;
  params: name=value, with __ for the scope separator"""
  from iotbx.cli_parser import CCTBXParser
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb_str)
  args = [prefix + ".pdb"]
  if cif_str is not None:
    with open(prefix + ".cif", "w") as f:
      f.write(cif_str)
    args.append(prefix + ".cif")
  args += ["%s=%s" % (name.replace("__", "."), value)
           for name, value in params.items()]
  parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
  parser.parse_args(args)
  log = StringIO()
  prog = reduce2.Program(parser.data_manager, parser.working_phil.extract(),
    logger=log)
  return prog, log

def test_output_names():
  # 13: default output and description file names
  for option in OPTIONS:
    for flips, pad in ((False, "H"), (True, "FH")):
      prog, log = make_program("tst_reduce2_neutron_names", pdb_1yjp,
        neutron_option=option, add_flip_movers=flips)
      prog.validate()
      suffix = "" if option == "none" else "_" + option
      expected = "tst_reduce2_neutron_names" + pad + suffix
      assert prog.params.output.filename == expected + ".pdb", (
        option, flips, prog.params.output.filename)
      assert prog.params.output.description_file_name == expected + ".txt"

def check_sorry(prefix, pdb_str, text, **params):
  prog, log = make_program(prefix, pdb_str, **params)
  try:
    prog.validate()
  except Sorry as e:
    assert text in str(e), str(e)
  else:
    raise AssertionError("no Sorry for %s" % params)

pdb_ser_dg = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      2  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      3  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      4  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM      5  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM      6  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM      7  DG  SER A   2      15.856   0.874   6.441  1.00 20.00           D
END
"""

def test_validate_sorrys():
  prefix = "tst_reduce2_neutron_validate"
  for option in OPTIONS[1:]:
    check_sorry(prefix, pdb_1yjp, "cannot be used with approach=remove",
      neutron_option=option, approach="remove")
    check_sorry(prefix, pdb_1yjp, "cannot be used with qmready=True",
      neutron_option=option, qmready=True)
    # interim guards until step 4
    check_sorry(prefix, pdb_1yjp,
      "neutron_option with approach=optimize is not supported yet",
      neutron_option=option, approach="optimize")
    check_sorry(prefix, pdb_ser_dg, "keep_existing_H=True and D in the input "
      "is not supported yet", neutron_option=option, keep_existing_H=True)
    # H-only input with keep_existing_H=True is fine
    prog, log = make_program(prefix, pdb_1yjp, neutron_option=option,
      keep_existing_H=True)
    prog.validate()
  # none: the same runs pass validate()
  for params in (dict(approach="remove"), dict(qmready=True),
                 dict(approach="optimize"), dict(keep_existing_H=True)):
    prog, log = make_program(prefix, pdb_ser_dg, neutron_option="none",
      **params)
    prog.validate()

def test_optimize_guard_reached():
  # the optimize guard stops a real run before any work is done
  try:
    run_reduce2("tst_reduce2_neutron_guard", pdb_1yjp,
      ["approach=optimize", "neutron_option=hd_and_h"])
  except Sorry as e:
    assert "not supported yet" in str(e)
  else:
    raise AssertionError("approach=optimize neutron_option=hd_and_h ran")

def test_none_distances():
  # 9b: none keeps today's distances; no D anywhere
  for args, neutron in (([], False), (["use_neutron_distances=True"], True),
                        (["qmready=True"], True)):
    fn, log, description = run_reduce2("tst_reduce2_neutron_9b", pdb_1yjp,
      ["neutron_option=none"] + args)
    check_bond_lengths(fn, use_neutron_distances=neutron)
    h = iotbx.pdb.input(fn).construct_hierarchy()
    assert "D" not in [a.element.strip() for a in h.atoms()]
    assert "DOD" not in [ag.resname for ag in h.atom_groups()]

def test_resolved_flag_all_h():
  # neutron_option switches neutron distances on; output still H-only
  fn, log, description = run_reduce2("tst_reduce2_neutron_all_h", pdb_1yjp,
    ["neutron_option=all_h", "approach=place", "exclude_water=False"])
  assert "neutron_option=all_h: use_neutron_distances set to True" in log
  check_water(fn, 0.98, 7)
  check_bond_lengths(fn, use_neutron_distances=True)
  h = iotbx.pdb.input(fn).construct_hierarchy()
  assert "D" not in [a.element.strip() for a in h.atoms()]

# ligand unknown to the CCD and the monomer library; C1-O1 1.4321 A
pdb_lig = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  C1  ZQ1 A   1      10.000  10.000  10.000  1.00 20.00           C
HETATM    2  O1  ZQ1 A   1      11.432  10.000  10.000  1.00 20.00           O
HETATM    3  H11 ZQ1 A   1       9.640  11.020  10.000  1.00 20.00           H
HETATM    4  H12 ZQ1 A   1       9.640   9.490  10.880  1.00 20.00           H
HETATM    5  H13 ZQ1 A   1       9.640   9.490   9.120  1.00 20.00           H
HETATM    6  HO1 ZQ1 A   1      11.750  10.900  10.000  1.00 20.00           H
END
"""
cif_lig = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZQ1  ZQ1  'test methanol'  ligand  6  2  .

data_comp_ZQ1
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZQ1  C1   C  CH3  0   0.000   0.000   0.000
 ZQ1  O1   O  OH1  0   1.432   0.000   0.000
 ZQ1  H11  H  HCH3 0  -0.360   1.020   0.000
 ZQ1  H12  H  HCH3 0  -0.360  -0.510   0.880
 ZQ1  H13  H  HCH3 0  -0.360  -0.510  -0.880
 ZQ1  HO1  H  HOH1 0   1.750   0.900   0.000

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 ZQ1  C1  O1   single  1.4321  0.020  1.4321
 ZQ1  C1  H11  single  0.970   0.020  1.090
 ZQ1  C1  H12  single  0.970   0.020  1.090
 ZQ1  C1  H13  single  0.970   0.020  1.090
 ZQ1  O1  HO1  single  0.850   0.020  0.980

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 ZQ1  H11  C1  H12  109.5  3.0
 ZQ1  H11  C1  H13  109.5  3.0
 ZQ1  H12  C1  H13  109.5  3.0
 ZQ1  O1   C1  H11  109.5  3.0
 ZQ1  O1   C1  H12  109.5  3.0
 ZQ1  O1   C1  H13  109.5  3.0
 ZQ1  C1   O1  HO1  109.5  3.0
"""

def test_rebuild_model():
  # _RebuildModel(): user CIF restraints, crystal symmetry, atom count
  prog, log = make_program("tst_reduce2_neutron_rebuild", pdb_lig, cif_lig,
    stop_on_any_missing_hydrogen=True)
  prog.validate()
  prog.model = prog.data_manager.get_model()
  hierarchy = prog.model.get_hierarchy().deep_copy()
  for make_restraints in (True, False):
    rebuilt = prog._RebuildModel(hierarchy.deep_copy(), make_restraints)
    assert rebuilt.get_number_of_atoms() == prog.model.get_number_of_atoms()
    assert rebuilt.crystal_symmetry().is_similar_symmetry(
      prog.model.crystal_symmetry())
    assert (rebuilt.get_restraints_manager() is not None) == make_restraints
  # the CIF's C1-O1 target is used: the restraints come from the user CIF
  simple, asu = prog._RebuildModel(hierarchy.deep_copy(),
    True).get_restraints_manager().geometry.get_all_bond_proxies()
  ideal = sorted(round(p.distance_ideal, 4) for p in simple)
  assert 1.4321 in ideal, ideal
  # without the CIF the ligand cannot be interpreted
  prog, log = make_program("tst_reduce2_neutron_rebuild_nocif", pdb_lig,
    stop_on_any_missing_hydrogen=True)
  prog.validate()
  prog.model = prog.data_manager.get_model()
  try:
    prog._RebuildModel(prog.model.get_hierarchy().deep_copy(), True)
  except Sorry:
    pass
  else:
    raise AssertionError("ZQ1 interpreted without its CIF")

# ------------------------------------------------------------------------------
# step 3.2: output step for add and place

NEUTRON = ("all_h", "all_d", "hd_and_h", "hd_and_d", "all_hd")

def run_prog(prefix, pdb_str, cif_str=None, **params):
  """validate() and run() a reduce2.Program; returns the program"""
  params.setdefault("output__filename", prefix + "_out.pdb")
  params.setdefault("output__description_file_name", prefix + ".txt")
  params.setdefault("output__overwrite", True)
  prog, log = make_program(prefix, pdb_str, cif_str, **params)
  prog.validate()
  prog.run()
  prog.log_text = log.getvalue()
  return prog

def hd_pairs(hierarchy):
  """[(H, D)] from exchangeable_hd_selections()"""
  atoms = hierarchy.atoms()
  result = []
  for x in hierarchy.exchangeable_hd_selections():
    a, b = atoms[x[0][0]], atoms[x[1][0]]
    if a.element.strip() == "D": a, b = b, a
    result.append((a, b))
  return result

def heavy_atoms(hierarchy):
  from mmtbx.hydrogens import hd_exchange
  return dict((hd_exchange.atom_key(a), (a.name, a.parent().altloc, a.xyz,
    a.occ, a.b, tuple(a.uij) if a.uij_is_defined() else None))
    for a in hierarchy.atoms() if a.element.strip() not in ("H", "D"))

def test_distances_all_modes():
  # neutron targets for every mode, add and place; water 0.98 A; H and D of
  # a pair at the same position
  for option in NEUTRON:
    for approach in ("add", "place"):
      prefix = "tst_reduce2_neutron_32_%s_%s" % (option, approach)
      fn, log, description = run_reduce2(prefix, pdb_1yjp,
        ["neutron_option=" + option, "approach=" + approach,
         "exclude_water=False"])
      check_bond_lengths(fn, use_neutron_distances=True)
      check_water(fn, 0.98, 7)
      h = iotbx.pdb.input(fn).construct_hierarchy()
      waters = set(ag.resname for ag in h.atom_groups()
                   if ag.resname in ("HOH", "DOD"))
      assert waters == set(["HOH" if option == "all_h" else "DOD"]), waters
      pairs = hd_pairs(h)
      assert (len(pairs) > 0) == (option in ("hd_and_h", "hd_and_d",
                                             "all_hd")), (option, len(pairs))
      for hh, dd in pairs:
        assert approx_equal(hh.xyz, dd.xyz, eps=1e-3), (hh.id_str())

def test_n_terminal_nh3():
  # 1: Gly 1 H1, H2, H3 all split in hd_and_h
  fn, log, description = run_reduce2("tst_reduce2_neutron_nh3", pdb_1yjp,
    ["neutron_option=hd_and_h"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  gly = [rg for rg in h.residue_groups() if rg.resseq.strip() == "1"][0]
  by_alt = dict((ag.altloc, sorted(a.name.strip() for a in ag.atoms()))
                for ag in gly.atom_groups())
  assert by_alt["A"] == ["H1", "H2", "H3"], by_alt
  assert by_alt["B"] == ["D1", "D2", "D3"], by_alt

def test_conformers():
  # 8: Ser with OG in A and B: backbone H on the blank N split into C/D,
  # HG in A and B written as D and listed as conformer sites
  from mmtbx.hydrogens.tst_hd_exchange import pdb_pep_og_ab
  for approach in ("add", "place"):
    prog = run_prog("tst_reduce2_neutron_conf_" + approach, pdb_pep_og_ab,
      neutron_option="hd_and_h", approach=approach)
    h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
    ser = [rg for rg in h.residue_groups() if rg.resseq.strip() == "2"][0]
    hd = dict((ag.altloc, sorted(a.name.strip() for a in ag.atoms()
               if a.element.strip() in ("H", "D"))) for ag in ser.atom_groups())
    assert hd["C"] == ["H"] and hd["D"] == ["D"], hd
    assert "DG" in hd["A"] and "DG" in hd["B"], hd
    sites = prog._neutron_result.conformer_sites
    assert sorted(x.site[5:9].strip() for x in sites
                  if x.reason == "conformer") == ["HG", "HG"], sites
    # 8 split sites and the two conformer D (water excluded)
    assert prog._neutron_result.n_exchangeable == 10
    with open(prog.params.output.description_file_name) as f:
      assert "Exchangeable H sites (not water): 10" in f.read()

def test_output_symmetry():
  # 9d: none in -> none out; cryo-EM placeholder -> 1 1 1 P 1; cell -> cell
  no_cs = "\n".join(l for l in pdb_1yjp.splitlines()
                    if not l.startswith("CRYST1"))
  placeholder = pdb_1yjp.replace(
    "CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1",
    "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1")
  for label, pdb_str, cryst1 in (("nocs", no_cs, None),
      ("placeholder", placeholder, "CRYST1    1.000    1.000    1.000"),
      ("cell", pdb_1yjp, "CRYST1   21.937    4.866   23.477")):
    fn, log, description = run_reduce2("tst_reduce2_neutron_cs_" + label,
      pdb_str, ["neutron_option=hd_and_h"])
    with open(fn) as f:
      lines = [l for l in f.read().splitlines() if l.startswith("CRYST1")]
    if cryst1 is None:
      assert lines == [], lines
    else:
      assert len(lines) == 1 and lines[0].startswith(cryst1), lines

def test_output_model():
  # 9e, 9f, 9h: results and heavy atoms
  prog = run_prog("tst_reduce2_neutron_out", pdb_1yjp,
    neutron_option="hd_and_h", exclude_water=False)
  results = prog.get_results()
  assert results.model is prog._output_model
  assert hasattr(results, "total_charge")
  out_h = results.model.get_hierarchy()
  assert "D" in [a.element.strip() for a in out_h.atoms()]
  # self.model stays H-only
  assert "D" not in [a.element.strip()
                     for a in prog.model.get_hierarchy().atoms()]
  # 9h: heavy atoms identical (xyz to round-off from rebuilding the model)
  out_heavy = heavy_atoms(out_h)
  in_heavy = heavy_atoms(prog.model.get_hierarchy())
  assert sorted(out_heavy) == sorted(in_heavy)
  for key, (name, altloc, xyz, occ, b, uij) in in_heavy.items():
    o = out_heavy[key]
    assert (o[0], o[1], o[3], o[4], o[5]) == (name, altloc, occ, b, uij), key
    assert approx_equal(o[2], xyz, eps=1e-9, out=None), key
  # 9e: TYR ring names as in self.model
  def ring(h):
    return dict((a.name.strip(), a.xyz) for a in h.atoms()
      if a.parent().resname == "TYR" and a.name.strip() in
      ("CD1", "CD2", "CE1", "CE2"))
  assert ring(out_h) == ring(prog.model.get_hierarchy())

def test_write_files_false():
  # 9g: no files, and the results model is the H/D model
  prefix = "tst_reduce2_neutron_nowrite"
  for ext in ("_out.pdb", ".txt"):
    if os.path.exists(prefix + ext): os.remove(prefix + ext)
  prog = run_prog(prefix, pdb_1yjp, neutron_option="hd_and_h",
    output__write_files=False)
  assert not os.path.exists(prefix + "_out.pdb")
  assert not os.path.exists(prefix + ".txt")
  model = prog.get_results().model
  assert "D" in [a.element.strip() for a in model.get_hierarchy().atoms()]

def test_water_after_select_hd():
  # 9k(a): placement drops the unparameterizable Lys HD2/HD3 (select());
  # the water O-D stays 0.98 A
  fn, log, description = run_reduce2("tst_reduce2_neutron_9ka",
    pdb_truncated_lys, ["neutron_option=hd_and_h", "exclude_water=False"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  names = [a.name.strip() for a in h.atoms() if a.parent().resname == "LYS"]
  assert "HD2" not in names and "HD3" not in names, names
  assert [ag.resname for ag in h.atom_groups()
          if ag.resname in ("HOH", "DOD")] == ["DOD"]
  check_water(fn, 0.98, 1)

def test_occupancy_groups():
  # 10: one exchangeable pair per split site; phenix.refine's
  # occupancy_selections() puts each split site in one two-atom group
  from mmtbx.refinement.occupancies import occupancy_selections
  prog = run_prog("tst_reduce2_neutron_occ", pdb_1yjp,
    neutron_option="hd_and_h")
  fn = prog.params.output.filename
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(fn), log=null_out())
  model.process(make_restraints=True)
  h = model.get_hierarchy()
  pairs = h.exchangeable_hd_selections()
  assert len(pairs) == prog._neutron_result.n_split
  groups = occupancy_selections(model=model, as_flex_arrays=False)
  for x in pairs:
    i_h, i_d = x[0][0], x[1][0]
    found = [g for g in groups if any(i_h in s for s in g)]
    assert len(found) == 1, (h.atoms()[i_h].id_str(), found)
    members = sorted(i for s in found[0] for i in s)
    assert members == sorted([i_h, i_d]), (h.atoms()[i_h].id_str(), members)

def test_partial_parent():
  # 12: Tyr OH at 0.8 -> HH 0.4 / DH 0.4, listed
  pdb_str = pdb_1yjp.replace(
    "OH  TYR A   7       3.766   0.589  10.291  1.00",
    "OH  TYR A   7       3.766   0.589  10.291  0.80")
  prog = run_prog("tst_reduce2_neutron_partial", pdb_str,
    neutron_option="hd_and_h")
  h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
  occ = dict((a.name.strip(), a.occ) for a in h.atoms()
             if a.parent().resname == "TYR" and a.name.strip() in ("HH", "DH"))
  assert approx_equal(occ["HH"], 0.4) and approx_equal(occ["DH"], 0.4), occ
  sites = prog._neutron_result.partial_parent_sites
  assert [(x.site[5:9].strip(), x.parent_occ) for x in sites] == [
    ("HH", 0.8)], sites

def test_ambiguous_parent():
  # 5b: an H/D with two covalent partners stops the run before the Optimizer.
  # reduce2 takes no geometry_restraints.edits, and placement drops an H
  # that a CIF bonds to two heavy atoms, so the check is called directly on
  # a model whose edits make Ser HG ambiguous
  from mmtbx.hydrogens import tst_hd_exchange as t
  model = t.make_model(t.pdb_ser_hd.replace(
    "ATOM     11  DG BSER A   2      15.856   0.874   6.441  0.50 20.00"
    "           D\n", ""), [t.edits_bond("name HG", "name CB", 1.1)])
  prog, log = make_program("tst_reduce2_neutron_amb", pdb_1yjp,
    neutron_option="hd_and_h")
  prog.validate()
  prog.model = model
  try:
    prog._CheckAmbiguousParents()
  except Sorry as e:
    text = str(e)
    assert "more than one covalent bond partner" in text
    assert 'HG ASER' in text and 'OG  SER' in text and 'CB  SER' in text
    assert "covalent geometry" in text and "edits" in text
  else:
    raise AssertionError("ambiguous HG not reported")
  # a normal model passes the check
  prog.model = make_model_from_str(pdb_1yjp_h())
  prog._CheckAmbiguousParents()
  # end to end, a CIF bond from HO1 to a second heavy atom: placement drops
  # HO1 (not parameterizable), so the run completes with none and hd_and_h
  cif = cif_lig.replace(" ZQ1  O1  HO1  single  0.850   0.020  0.980\n",
    " ZQ1  O1  HO1  single  0.850   0.020  0.980\n"
    " ZQ1  C1  HO1  single  1.050   0.020  1.050\n")
  assert cif != cif_lig
  for option in ("none", "hd_and_h"):
    prog = run_prog("tst_reduce2_neutron_amb_" + option, pdb_lig, cif,
      neutron_option=option, keep_existing_H=True)
    h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
    assert "HO1" not in [a.name.strip() for a in h.atoms()]

_pdb_1yjp_h = []
def pdb_1yjp_h():
  if not _pdb_1yjp_h:
    fn, _, _ = run_reduce2("tst_reduce2_neutron_1yjp_h", pdb_1yjp, [])
    with open(fn) as f: _pdb_1yjp_h.append(f.read())
  return _pdb_1yjp_h[0]

def make_model_from_str(pdb_str):
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(
    lines=pdb_str.split("\n"), source_info=None), log=null_out())
  model.process(make_restraints=True)
  return model

def test_flipkins_hd():
  # hd_and_h with flipkins: completes, the written model is the H/D model,
  # and the flipkins are those of the H-only run with the same (neutron)
  # distances: they are made from self.model
  kins = {}
  for option, args in (("none", dict(use_neutron_distances=True)),
                       ("hd_and_h", {})):
    d = "tst_reduce2_neutron_fk_" + option
    if not os.path.isdir(d): os.mkdir(d)
    prefix = "tst_reduce2_neutron_fkrun_" + option
    prog = run_prog(prefix, pdb_1yjp, neutron_option=option,
      add_flip_movers=True, output__flipkin_directory=d, **args)
    h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
    has_d = "D" in [a.element.strip() for a in h.atoms()]
    assert has_d == (option != "none")
    kins[option] = {}
    for f in sorted(os.listdir(d)):
      with open(os.path.join(d, f)) as fh:
        text = fh.read().replace(prefix, "BASE")
      kins[option][f.replace(prefix, "BASE")] = [l for l in text.splitlines()
        if not l.startswith("@caption")]
  assert len(kins["none"]) > 0
  assert kins["none"] == kins["hd_and_h"]

def test_output_check_fires():
  # a provenance map missing an entry is an internal error
  from mmtbx.hydrogens import hd_exchange
  original = hd_exchange.apply_neutron_option
  def broken(*args, **kwargs):
    r = original(*args, **kwargs)
    r.provenance.pop(sorted(r.provenance)[0])
    return r
  hd_exchange.apply_neutron_option = broken
  try:
    run_reduce2("tst_reduce2_neutron_broken", pdb_1yjp,
      ["neutron_option=hd_and_h"])
  except Sorry as e:
    assert "Internal error: H/D output check failed" in str(e), str(e)
    assert "no provenance" in str(e)
  else:
    raise AssertionError("output check did not fire")
  finally:
    hd_exchange.apply_neutron_option = original

pdb_1yjp_dod = pdb_1yjp.replace("HOH A", "DOD A")

def test_dod_input_all_h():
  # input waters named DOD (O only), all_h: placement adds H, the written
  # waters are HOH with H and there is no DOD anywhere
  assert pdb_1yjp_dod.count("DOD A") == 7
  for approach in ("add", "place"):
    fn, log, description = run_reduce2("tst_reduce2_neutron_dod_" + approach,
      pdb_1yjp_dod, ["neutron_option=all_h", "approach=" + approach,
      "exclude_water=False", "keep_existing_H=False"])
    h = iotbx.pdb.input(fn).construct_hierarchy()
    waters = [ag for ag in h.atom_groups() if ag.resname in ("HOH", "DOD")]
    assert len(waters) == 7
    assert set(ag.resname for ag in waters) == set(["HOH"])
    assert set(a.element.strip() for ag in waters for a in ag.atoms()) == \
      set(["O", "H"])
    assert "DOD" not in open(fn).read()
    check_water(fn, 0.98, 7)

def test_output_check_water_names():
  # the water-name check: renaming left out, a water with H named DOD
  from mmtbx.hydrogens import hd_exchange
  original = hd_exchange.apply_neutron_option
  def no_rename(hierarchy, *args, **kwargs):
    r = original(hierarchy, *args, **kwargs)
    for ag in hierarchy.atom_groups():
      if ag.resname == "HOH": ag.resname = "DOD"
    return r
  hd_exchange.apply_neutron_option = no_rename
  try:
    run_reduce2("tst_reduce2_neutron_wname", pdb_1yjp_dod,
      ["neutron_option=all_h", "exclude_water=False"])
  except Sorry as e:
    assert "Internal error: H/D output check failed" in str(e), str(e)
    assert "water name DOD does not match its H/D" in str(e), str(e)
  else:
    raise AssertionError("water-name check did not fire")
  finally:
    hd_exchange.apply_neutron_option = original

RENAMED = "water(s) renamed DOD -> HOH (their D were replaced by H)"

def water_groups(fn):
  h = iotbx.pdb.input(fn).construct_hierarchy()
  return [(ag.resname, sorted(set(a.element.strip() for a in ag.atoms())))
          for ag in h.atom_groups() if ag.resname in ("HOH", "DOD")]

def test_dod_input_none():
  # none, add and place: placement replaced the D, so the waters are HOH
  for approach in ("add", "place"):
    fn, log, description = run_reduce2("tst_reduce2_neutron_dodn_" + approach,
      pdb_1yjp_dod, ["approach=" + approach, "exclude_water=False"])
    assert water_groups(fn) == [("HOH", ["H", "O"])] * 7, water_groups(fn)
    with open(fn) as f:
      assert "DOD" not in f.read()
    assert "7 " + RENAMED in log and "7 " + RENAMED in description

def test_dod_input_none_exclude_water():
  # no H added to water: the waters stay DOD with O only, nothing logged
  fn, log, description = run_reduce2("tst_reduce2_neutron_dodx", pdb_1yjp_dod,
    [])
  assert water_groups(fn) == [("DOD", ["O"])] * 7, water_groups(fn)
  assert RENAMED not in log and RENAMED not in description

def test_dod_input_with_d_kept():
  # keep_existing_H=True and a DOD water with D1/D2: stays DOD with its D
  # (today's behaviour until step 4)
  pdb_str = pdb_1yjp_dod.replace(
    "HETATM   61  O   DOD A   8      -6.471   5.227   7.124  1.00 22.62"
    "           O\n",
    "HETATM   61  O   DOD A   8      -6.471   5.227   7.124  1.00 22.62"
    "           O\n"
    "HETATM   61  D1  DOD A   8      -5.621   5.227   7.124  1.00 22.62"
    "           D\n"
    "HETATM   61  D2  DOD A   8      -6.675   6.052   7.124  1.00 22.62"
    "           D\n")
  assert pdb_str.count("D1  DOD") == 1
  fn, log, description = run_reduce2("tst_reduce2_neutron_dodd", pdb_str,
    ["keep_existing_H=True", "exclude_water=False"])
  # existing behaviour: placement keeps D1/D2 and also adds H1/H2 (it does
  # not match D1/D2 to the dictionary's H1/H2); the water has D, so the
  # rename leaves it as DOD; the six other waters are renamed
  groups = water_groups(fn)
  assert groups.count(("DOD", ["D", "H", "O"])) == 1, groups
  assert groups.count(("HOH", ["H", "O"])) == 6, groups
  h = iotbx.pdb.input(fn).construct_hierarchy()
  w8 = sorted(a.name.strip() for a in h.atoms()
              if a.parent().parent().resseq.strip() == "8")
  assert w8 == ["D1", "D2", "H1", "H2", "O"], w8
  assert "6 " + RENAMED in log

def test_dod_input_hd_and_h():
  # the earlier rename to HOH does not break the neutron path or the
  # resname invariant: written as DOD with D, output check passes
  fn, log, description = run_reduce2("tst_reduce2_neutron_dodhd",
    pdb_1yjp_dod, ["neutron_option=hd_and_h", "exclude_water=False"])
  assert water_groups(fn) == [("DOD", ["D", "O"])] * 7, water_groups(fn)
  check_water(fn, 0.98, 7)

def test_dod_input_flipkins():
  d = "tst_reduce2_neutron_dodfk"
  if not os.path.isdir(d): os.mkdir(d)
  fn, log, description = run_reduce2("tst_reduce2_neutron_dodfk_run",
    pdb_1yjp_dod, ["exclude_water=False", "add_flip_movers=True",
    "output.flipkin_directory=" + d])
  assert water_groups(fn) == [("HOH", ["H", "O"])] * 7, water_groups(fn)
  assert len(os.listdir(d)) > 0

# ------------------------------------------------------------------------------
# step 3.3: per-mode labels, His, ligands, glycan, mmCIF, description

def atom_label_set(src, is_file=True):
  """{(resseq, name, altloc, element, occupancy)}"""
  inp = iotbx.pdb.input(file_name=src) if is_file else iotbx.pdb.input(
    source_info=None, lines=src)
  h = inp.construct_hierarchy()
  return set((a.parent().parent().resseq.strip(), a.name.strip(),
              a.parent().altloc, a.element.strip(), round(a.occ, 2))
             for a in h.atoms())

def hd_table(text):
  """{(resseq, name, altloc, element, occupancy)} from lines
  'resseq name altloc element occupancy', altloc '-' for blank"""
  result = set()
  for line in text.strip().splitlines():
    resseq, name, altloc, element, occ = line.split()
    result.add((resseq, name, "" if altloc == "-" else altloc, element,
                float(occ)))
  return result

def hd_labels(fn):
  return set(x for x in atom_label_set(fn) if x[3] in ("H", "D"))

# Ser-Leu dipeptide; Leu stops after CA
dipeptide = """\
CRYST1   14.228   15.246   13.577  90.00  90.00  90.00 P 1
ATOM    577  N   SER A  37       8.872   5.000   5.995  1.00 25.33           N
ATOM    578  CA  SER A  37       7.827   5.915   6.440  1.00 31.94           C
ATOM    579  CB  SER A  37       7.709   5.892   7.965  1.00 42.16           C
ATOM    580  OG  SER A  37       8.926   6.283   8.577  1.00 54.74           O
ATOM    581  C   SER A  37       8.095   7.336   5.956  1.00 32.18           C
ATOM    582  O   SER A  37       9.228   7.814   6.002  1.00 31.99           O
ATOM    588  N   LEU A  38       7.045   8.005   5.492  1.00 34.04           N
ATOM    589  CA  LEU A  38       7.165   9.372   5.000  1.00 36.58           C
ATOM    594  C   LEU A  38       6.026  10.246   5.515  1.00 37.03           C
ATOM    595  O   LEU A  38       5.000   9.741   5.969  1.00 35.00           O
END
"""
# expected H/D per mode: the chain-start N-H (Ser 37 H, not resseq 1, so no
# NH3) and HA on the truncated Leu are placed like any other H
dipeptide_hd = {
"all_h": """
37 H - H 1.0
37 HA - H 1.0
37 HB2 - H 1.0
37 HB3 - H 1.0
37 HG - H 1.0
38 H - H 1.0
38 HA - H 1.0""",
"all_d": """
37 D - D 1.0
37 DA - D 1.0
37 DB2 - D 1.0
37 DB3 - D 1.0
37 DG - D 1.0
38 D - D 1.0
38 DA - D 1.0""",
"hd_and_h": """
37 H A H 0.5
37 D B D 0.5
37 HG A H 0.5
37 DG B D 0.5
37 HA - H 1.0
37 HB2 - H 1.0
37 HB3 - H 1.0
38 H A H 0.5
38 D B D 0.5
38 HA - H 1.0""",
"hd_and_d": """
37 H A H 0.5
37 D B D 0.5
37 HG A H 0.5
37 DG B D 0.5
37 DA - D 1.0
37 DB2 - D 1.0
37 DB3 - D 1.0
38 H A H 0.5
38 D B D 0.5
38 DA - D 1.0""",
"all_hd": """
37 H A H 0.5
37 D B D 0.5
37 HA A H 0.5
37 DA B D 0.5
37 HB2 A H 0.5
37 DB2 B D 0.5
37 HB3 A H 0.5
37 DB3 B D 0.5
37 HG A H 0.5
37 DG B D 0.5
38 H A H 0.5
38 D B D 0.5
38 HA A H 0.5
38 DA B D 0.5""",
}

def test_dipeptide_modes():
  # checked against phenix.ready_set on 2026-10-17; differences explained in
  # the neutron_option plan
  for option in NEUTRON:
    fn, log, description = run_reduce2("tst_reduce2_neutron_dipep_" + option,
      dipeptide, ["neutron_option=" + option])
    ours, expected = hd_labels(fn), hd_table(dipeptide_hd[option])
    assert ours == expected, (option, sorted(ours ^ expected))

his_no_h = """\
CRYST1   12.670   14.862   15.186  90.00  90.00  90.00 P 1
ATOM    780  N   HIS A  45       5.209   7.219   5.800  1.00  9.38           N
ATOM    781  CA  HIS A  45       6.450   6.586   6.229  1.00 10.30           C
ATOM    782  C   HIS A  45       6.629   5.221   5.572  1.00 11.06           C
ATOM    783  O   HIS A  45       7.669   4.940   4.976  1.00 12.76           O
ATOM    784  CB  HIS A  45       6.482   6.444   7.752  1.00 11.71           C
ATOM    785  CG  HIS A  45       6.466   7.752   8.482  1.00 17.20           C
ATOM    786  ND1 HIS A  45       6.374   7.837   9.854  1.00 28.53           N
ATOM    787  CD2 HIS A  45       6.532   9.027   8.030  1.00 22.78           C
ATOM    788  CE1 HIS A  45       6.383   9.107  10.217  1.00 26.38           C
ATOM    789  NE2 HIS A  45       6.478   9.850   9.129  1.00 20.23           N
END
"""
his_hd = """
45 H A H 0.5
45 D B D 0.5
45 HD1 A H 0.5
45 DD1 B D 0.5
45 HE2 A H 0.5
45 DE2 B D 0.5
45 HA - H 1.0
45 HB2 - H 1.0
45 HB3 - H 1.0
45 HD2 - H 1.0
45 HE1 - H 1.0"""

def test_his_split():
  # His in hd_and_h: HD1 and HE2 split into A/B with DD1/DE2. approach=place
  # keeps both ring N-H (no Optimizer, so no His flip or protonation choice).
  # checked against phenix.ready_set on 2026-10-17; differences explained in
  # the neutron_option plan
  fn, log, description = run_reduce2("tst_reduce2_neutron_his", his_no_h,
    ["neutron_option=hd_and_h", "approach=place"])
  ours, expected = hd_labels(fn), hd_table(his_hd)
  assert ours == expected, sorted(ours ^ expected)

hc4_no_h = """\
CRYST1   67.145   67.145   41.119  90.00  90.00 120.00 P 63          6
HETATM 1987  C1  HC4 A 169      12.270  -2.851   7.690  1.00  2.51           C
HETATM 1988  O1  HC4 A 169      12.412  -3.051   8.915  1.00  5.58           O
HETATM 1989  C2  HC4 A 169      13.073  -1.969   6.902  1.00  3.66           C
HETATM 1990  C3  HC4 A 169      13.948  -1.115   7.480  1.00  4.10           C
HETATM 1991  C1' HC4 A 169      14.960  -0.357   6.862  1.00  3.50           C
HETATM 1992  C2' HC4 A 169      15.734   0.515   7.642  1.00  7.24           C
HETATM 1993  C3' HC4 A 169      16.796   1.251   7.105  1.00  7.06           C
HETATM 1994  C4' HC4 A 169      17.093   1.090   5.743  1.00  6.80           C
HETATM 1995  C5' HC4 A 169      16.317   0.232   4.947  1.00  7.54           C
HETATM 1996  C6' HC4 A 169      15.255  -0.474   5.508  1.00  4.60           C
HETATM 1997  O4' HC4 A 169      18.097   1.748   5.166  1.00  7.92           O
END
"""
hc4_all_d = """
169 D2 - D 1.0
169 D3 - D 1.0
169 D2' - D 1.0
169 D3' - D 1.0
169 D5' - D 1.0
169 D6' - D 1.0
169 DO4' - D 1.0"""

def test_hc4_all_d():
  # HC4 ligand, all_d: 7 D and no H
  # checked against phenix.ready_set on 2026-10-17; differences explained in
  # the neutron_option plan
  fn, log, description = run_reduce2("tst_reduce2_neutron_hc4_d", hc4_no_h,
    ["neutron_option=all_d"])
  ours, expected = hd_labels(fn), hd_table(hc4_all_d)
  assert ours == expected, sorted(ours ^ expected)

def test_hc4_hd_and_h():
  # 6: ligand OH: HO4' split, the six C-H stay H
  fn, log, description = run_reduce2("tst_reduce2_neutron_hc4_hd",
    hc4_no_h, ["neutron_option=hd_and_h"])
  labels = atom_label_set(fn)
  hd = sorted((x[1], x[2]) for x in labels if x[3] in ("H", "D"))
  assert ("HO4'", "A") in hd and ("DO4'", "B") in hd, hd
  assert len([x for x in hd if x[1] == ""]) == 6, hd

def test_glycan_link():
  # 7: NAG-ASN (fixture of tst_add_hydrogen_16, 1a7s ASN 114 - NAG 401):
  # ASN ND2 H split, ring OH split, no H on the O1 the link replaces
  from mmtbx.hydrogens.tst_add_hydrogen_16 import pdb_str_019_control
  fn, log, description = run_reduce2("tst_reduce2_neutron_nag",
    pdb_str_019_control, ["neutron_option=hd_and_h"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  def hd(resname):
    return sorted((a.name.strip(), a.parent().altloc) for a in h.atoms()
      if a.parent().resname == resname and a.element.strip() in ("H", "D"))
  asn, nag = hd("ASN"), hd("NAG")
  assert ("HD21", "A") in asn and ("DD21", "B") in asn, asn
  assert not [x for x in asn if x[0][1:] == "D22"], asn
  for name in ("O3", "O4", "O6"):
    assert ("H" + name, "A") in nag and ("D" + name, "B") in nag, nag
  assert not [x for x in nag if x[0][1:] == "O1"], nag

def make_cif(code, atoms, bonds, angles):
  """monomer CIF: atoms (name, element, energy type, x, y, z), bonds
  (a, b, x-ray, neutron), angles (a, b, c, value)"""
  lines = ["data_comp_list", "loop_", "_chem_comp.id",
    "_chem_comp.three_letter_code", "_chem_comp.name", "_chem_comp.group",
    "_chem_comp.number_atoms_all", "_chem_comp.number_atoms_nh",
    "_chem_comp.desc_level",
    " %s  %s  'test ligand'  ligand  %d  %d  ." % (code, code, len(atoms),
      len([x for x in atoms if x[1] != "H"])), "", "data_comp_" + code,
    "loop_", "_chem_comp_atom.comp_id", "_chem_comp_atom.atom_id",
    "_chem_comp_atom.type_symbol", "_chem_comp_atom.type_energy",
    "_chem_comp_atom.charge", "_chem_comp_atom.x", "_chem_comp_atom.y",
    "_chem_comp_atom.z"]
  lines += [" %s  %-4s  %s  %-4s  0  %.3f  %.3f  %.3f" % ((code,) + x)
            for x in atoms]
  lines += ["", "loop_", "_chem_comp_bond.comp_id",
    "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2",
    "_chem_comp_bond.type", "_chem_comp_bond.value_dist",
    "_chem_comp_bond.value_dist_esd", "_chem_comp_bond.value_dist_neutron"]
  lines += [" %s  %-4s  %-4s  single  %.3f  0.020  %.3f" % ((code,) + x)
            for x in bonds]
  lines += ["", "loop_", "_chem_comp_angle.comp_id",
    "_chem_comp_angle.atom_id_1", "_chem_comp_angle.atom_id_2",
    "_chem_comp_angle.atom_id_3", "_chem_comp_angle.value_angle",
    "_chem_comp_angle.value_angle_esd"]
  lines += [" %s  %-4s  %-4s  %-4s  %.1f  3.0" % ((code,) + x)
            for x in angles]
  return "\n".join(lines) + "\n"

def heavy_pdb(code, atoms, offset=(10.0, 10.0, 10.0), cryst1=True):
  lines = []
  if cryst1:
    lines.append("CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1")
  for i, (name, element, energy, x, y, z) in enumerate(atoms):
    if element == "H": continue
    lines.append("HETATM%5d %-4s %s A   1    %8.3f%8.3f%8.3f  1.00 20.00"
      "          %2s" % (i + 1, name if len(name) == 4 else " " + name, code,
      x + offset[0], y + offset[1], z + offset[2], element))
  return "\n".join(lines) + "\n"

def _tetra(center, a, b):
  """angles for the H on a CH2 between a and b"""
  return [(a, center, "H%s1" % center[-1]), (a, center, "H%s2" % center[-1]),
          (b, center, "H%s1" % center[-1]), (b, center, "H%s2" % center[-1]),
          ("H%s1" % center[-1], center, "H%s2" % center[-1])]

# ZQ3: ethylene glycol; one hydroxyl H named 1HO2 (digit first)
zq3_atoms = [
  ("O2", "O", "OH1", -0.700, 1.210, 0.0), ("C1", "C", "CH2", 0.0, 0.0, 0.0),
  ("C2", "C", "CH2", 1.520, 0.0, 0.0), ("O3", "O", "OH1", 2.220, 1.210, 0.0),
  ("1HO2", "H", "HOH1", -1.600, 1.000, 0.0),
  ("H11", "H", "HCH2", -0.330, -0.550, 0.890),
  ("H12", "H", "HCH2", -0.330, -0.550, -0.890),
  ("H21", "H", "HCH2", 1.850, -0.550, 0.890),
  ("H22", "H", "HCH2", 1.850, -0.550, -0.890),
  ("HO3", "H", "HOH1", 3.120, 1.000, 0.0)]
zq3_cif = make_cif("ZQ3", zq3_atoms,
  [("C1", "O2", 1.43, 1.43), ("C1", "C2", 1.52, 1.52),
   ("C2", "O3", 1.43, 1.43), ("O2", "1HO2", 0.85, 0.98),
   ("O3", "HO3", 0.85, 0.98), ("C1", "H11", 0.97, 1.09),
   ("C1", "H12", 0.97, 1.09), ("C2", "H21", 0.97, 1.09),
   ("C2", "H22", 0.97, 1.09)],
  [(x[0], x[1], x[2], 109.5) for x in
    [("O2", "C1", "C2"), ("C1", "C2", "O3"), ("C1", "O2", "1HO2"),
     ("C2", "O3", "HO3")] + _tetra("C1", "O2", "C2") +
    _tetra("C2", "C1", "O3")])

# ZQ4: N-methyl ethanolamine, NH and OH (a terminal NH2 has no dihedral
# anchor for its riding H and is not placed: one heavy neighbour)
zq4_atoms = [
  ("C3", "C", "CH3", -1.430, 1.800, 0.0), ("N1", "N", "NH1", -0.700, 1.210, 0.0),
  ("C1", "C", "CH2", 0.0, 0.0, 0.0), ("C2", "C", "CH2", 1.520, 0.0, 0.0),
  ("O1", "O", "OH1", 2.220, 1.210, 0.0),
  ("HN1", "H", "HNH1", -0.200, 1.900, 0.500),
  ("H31", "H", "HCH3", -2.400, 1.400, 0.0),
  ("H32", "H", "HCH3", -1.300, 2.400, 0.880),
  ("H33", "H", "HCH3", -1.300, 2.400, -0.880),
  ("H11", "H", "HCH2", -0.330, -0.550, 0.890),
  ("H12", "H", "HCH2", -0.330, -0.550, -0.890),
  ("H21", "H", "HCH2", 1.850, -0.550, 0.890),
  ("H22", "H", "HCH2", 1.850, -0.550, -0.890),
  ("HO1", "H", "HOH1", 3.120, 1.000, 0.0)]
zq4_cif = make_cif("ZQ4", zq4_atoms,
  [("C3", "N1", 1.47, 1.47), ("C1", "N1", 1.47, 1.47),
   ("C1", "C2", 1.52, 1.52), ("C2", "O1", 1.43, 1.43),
   ("N1", "HN1", 0.89, 1.01), ("O1", "HO1", 0.85, 0.98),
   ("C3", "H31", 0.96, 1.09), ("C3", "H32", 0.96, 1.09),
   ("C3", "H33", 0.96, 1.09), ("C1", "H11", 0.97, 1.09),
   ("C1", "H12", 0.97, 1.09), ("C2", "H21", 0.97, 1.09),
   ("C2", "H22", 0.97, 1.09)],
  [(x[0], x[1], x[2], 109.5) for x in
    [("C3", "N1", "C1"), ("N1", "C1", "C2"), ("C1", "C2", "O1"),
     ("C3", "N1", "HN1"), ("C1", "N1", "HN1"), ("C2", "O1", "HO1"),
     ("N1", "C3", "H31"), ("N1", "C3", "H32"), ("N1", "C3", "H33"),
     ("H31", "C3", "H32"), ("H31", "C3", "H33"), ("H32", "C3", "H33")] +
    _tetra("C1", "N1", "C2") + _tetra("C2", "C1", "O1")])

def test_user_cif_ligand():
  # 9c (add): invented ligand with OH and NH2 from a CIF written here: the
  # polar H are split, no missing restraints; without the CIF the run fails
  prog = run_prog("tst_reduce2_neutron_zq4", heavy_pdb("ZQ4", zq4_atoms),
    zq4_cif, neutron_option="hd_and_h")
  h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
  hd = sorted((a.name.strip(), a.parent().altloc) for a in h.atoms()
              if a.element.strip() in ("H", "D"))
  for name in ("HN1", "HO1"):
    assert (name, "A") in hd and ("D" + name[1:], "B") in hd, hd
  assert len([x for x in hd if x[1] == ""]) == 7, hd
  with open(prog.params.output.description_file_name) as f:
    description = f.read()
  for text in ("Restraints were not found", "Insufficient restraints",
               "not placed", "Could not find atom info"):
    assert text not in description and text not in prog.log_text, text
  try:
    run_prog("tst_reduce2_neutron_zq4_nocif", heavy_pdb("ZQ4", zq4_atoms),
      neutron_option="hd_and_h")
  except Sorry:
    pass
  else:
    raise AssertionError("ZQ4 ran without its CIF")

def test_name_not_pairable_end_to_end():
  # 9i: 1HO2 stays H (never D) and is listed; the other OH is split
  # (hd_and_h) or D (all_d); the written model reinterprets with restraints
  for option, ho3 in (("hd_and_h", "split"), ("all_d", "D")):
    prog = run_prog("tst_reduce2_neutron_zq3_" + option,
      heavy_pdb("ZQ3", zq3_atoms), zq3_cif, neutron_option=option)
    fn = prog.params.output.filename
    h = iotbx.pdb.input(fn).construct_hierarchy()
    hd = sorted((a.name.strip(), a.parent().altloc, a.element.strip())
                for a in h.atoms() if a.element.strip() in ("H", "D"))
    assert ("1HO2", "", "H") in hd, (option, hd)
    assert not [x for x in hd if x[0] == "1DO2"], hd
    if ho3 == "split":
      assert ("HO3", "A", "H") in hd and ("DO3", "B", "D") in hd, hd
    else:
      assert ("DO3", "", "D") in hd, hd
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    line = [l for l in description.splitlines() if "1HO2 ZQ3" in l]
    assert len(line) == 1 and "name not pairable by phenix.refine" in line[0]
    split = sorted(tuple(sorted(h.atoms()[i[0]].name.strip() for i in x))
                   for x in h.exchangeable_hd_selections())
    assert split == ([("DO3", "HO3")] if ho3 == "split" else []), split
    # reads back and reinterprets with restraints (and the CIF)
    from iotbx.cif import reader as cif_reader
    model = mmtbx.model.manager(model_input=iotbx.pdb.input(fn),
      restraint_objects=[("zq3.cif", cif_reader(
        input_string=zq3_cif).model())], log=null_out())
    model.process(make_restraints=True)
    from mmtbx.hydrogens import hd_exchange
    table = hd_exchange.h_parent_table(model)
    assert set(e.status for e in table.values()) == set(["ok"])

def test_section_in_log():
  # the H/D section also goes to the log, so it is there without files
  prefix = "tst_reduce2_neutron_zq3_log"
  for ext in ("_out.pdb", ".txt"):
    if os.path.exists(prefix + ext): os.remove(prefix + ext)
  prog = run_prog(prefix, heavy_pdb("ZQ3", zq3_atoms), zq3_cif,
    neutron_option="all_d", output__write_files=False)
  assert not os.path.exists(prefix + "_out.pdb")
  log = prog.log_text
  assert log.count("H/D model (neutron_option=all_d):") == 1
  line = [l for l in log.splitlines() if "1HO2 ZQ3" in l]
  assert len(line) == 1 and "kept as H" in line[0], line
  # the counts are printed once
  assert log.count("Split into H/D pairs") == 1

def test_mmcif_output():
  # 11: .cif output reads back with the PDB output's H/D labels and
  # occupancies and the same exchangeable pairs
  results = {}
  for ext in ("pdb", "cif"):
    fn, log, description = run_reduce2("tst_reduce2_neutron_mmcif_" + ext,
      pdb_1yjp, ["neutron_option=hd_and_h",
      "output.filename=tst_reduce2_neutron_mmcif_out." + ext])
    fn = "tst_reduce2_neutron_mmcif_out." + ext
    h = iotbx.pdb.input(fn).construct_hierarchy()
    labels = set((a.parent().parent().resseq.strip(), a.parent().resname,
      a.name.strip(), a.parent().altloc, a.element.strip(), round(a.occ, 2))
      for a in h.atoms())
    atoms = h.atoms()
    pairs = set(frozenset(atoms[i[0]].parent().parent().resseq.strip() + ":" +
      atoms[i[0]].name.strip() for i in x)
      for x in h.exchangeable_hd_selections())
    results[ext] = (labels, pairs)
  assert results["cif"][0] == results["pdb"][0], sorted(
    results["cif"][0] ^ results["pdb"][0])
  assert results["cif"][1] == results["pdb"][1]
  assert len(results["pdb"][1]) == 20

def test_description_file():
  # every category, by written name and with its reason: conformer site,
  # merged copies, partial parent, name not pairable, water renamed
  from mmtbx.hydrogens.tst_hd_exchange import pdb_pep_og_ab
  pdb_str = pdb_pep_og_ab.replace(
    "NZ  LYS A   3      19.531  -1.812  18.059  1.00",
    "NZ  LYS A   3      19.531  -1.812  18.059  0.80").replace("END\n", "")
  assert "0.80" in pdb_str
  pdb_str += heavy_pdb("ZQ3", zq3_atoms, offset=(5.0, 15.0, 5.0),
    cryst1=False) + "END\n"
  prog = run_prog("tst_reduce2_neutron_desc", pdb_str, zq3_cif,
    neutron_option="hd_and_h", exclude_water=False)
  with open(prog.params.output.description_file_name) as f:
    description = f.read()
  section = description[description.index("H/D model (neutron_option="
    "hd_and_h)"):]
  def entry(label):
    found = [l for l in section.splitlines() if l.strip().startswith(label)]
    assert len(found) == 1, (label, found)
    return found[0]
  assert "Neutron distances used (neutron_option=hd_and_h)" in section
  assert "Waters renamed DOD: 1" in section
  assert entry("DG ASER A   2").endswith("conformer")
  assert entry("DG BSER A   2").endswith("conformer")
  assert "merged conformer copies" in entry("H  CSER A   2")
  assert "merged conformer copies" in entry("HA  SER A   2")
  for name in ("HZ1", "HZ2", "HZ3"):
    assert "parent occupancy 0.80" in entry("%sALYS A   3" % name)
  assert "phenix.refine currently refines H + D = 1" in section
  assert "name not pairable by phenix.refine" in entry("1HO2 ZQ3 A   1")
  # none: no H/D section
  fn, log, description = run_reduce2("tst_reduce2_neutron_desc_none",
    pdb_1yjp, [])
  assert "H/D model" not in description

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
  test_floating_h_removed()
  test_floating_h_water_098()
  test_no_floating_h_unchanged()
  test_floating_h_flipkins()
  test_water_without_o_add()
  test_output_names()
  test_validate_sorrys()
  test_optimize_guard_reached()
  test_none_distances()
  test_resolved_flag_all_h()
  test_rebuild_model()
  test_distances_all_modes()
  test_n_terminal_nh3()
  test_conformers()
  test_output_symmetry()
  test_output_model()
  test_write_files_false()
  test_water_after_select_hd()
  test_occupancy_groups()
  test_partial_parent()
  test_ambiguous_parent()
  test_flipkins_hd()
  test_output_check_fires()
  test_dod_input_all_h()
  test_output_check_water_names()
  test_dod_input_none()
  test_dod_input_none_exclude_water()
  test_dod_input_with_d_kept()
  test_dod_input_hd_and_h()
  test_dod_input_flipkins()
  test_dipeptide_modes()
  test_his_split()
  test_hc4_all_d()
  test_hc4_hd_and_h()
  test_glycan_link()
  test_user_cif_ligand()
  test_name_not_pairable_end_to_end()
  test_section_in_log()
  test_mmcif_output()
  test_description_file()
  print("OK")
