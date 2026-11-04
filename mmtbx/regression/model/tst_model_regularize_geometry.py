from __future__ import absolute_import, division, print_function
import random
import iotbx.pdb
import mmtbx.model
from mmtbx.refinement.geometry_minimization import geometry_violations
from libtbx.utils import null_out
from scitbx.array_family import flex

# 1yjp with riding H (reduce2), geometry-minimized in its cell: no bond
# outliers and no clashes, in the asymmetric unit or with symmetry mates.
pdb_str = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.184   4.057   6.379  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.321   3.478   5.012  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.127   2.630   4.621  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.823   1.632   5.273  1.00 16.78           O
ATOM      5  H1  GLY A   1      -9.903   4.546   6.568  1.00 16.77           H
ATOM      6  H2  GLY A   1      -9.108   3.399   6.974  1.00 16.77           H
ATOM      7  H3  GLY A   1      -8.458   4.572   6.412  1.00 16.77           H
ATOM      8  HA2 GLY A   1     -10.115   2.923   4.977  1.00 16.57           H
ATOM      9  HA3 GLY A   1      -9.409   4.197   4.367  1.00 16.57           H
ATOM     10  N   ASN A   2      -7.449   3.028   3.547  1.00 15.02           N
ATOM     11  CA  ASN A   2      -6.278   2.293   3.093  1.00 14.10           C
ATOM     12  C   ASN A   2      -5.059   2.671   3.925  1.00 13.13           C
ATOM     13  O   ASN A   2      -4.965   3.782   4.454  1.00 11.91           O
ATOM     14  CB  ASN A   2      -6.016   2.576   1.614  1.00 15.38           C
ATOM     15  CG  ASN A   2      -7.194   2.208   0.733  1.00 14.08           C
ATOM     16  OD1 ASN A   2      -7.711   1.094   0.803  1.00 17.46           O
ATOM     17  ND2 ASN A   2      -7.625   3.146  -0.103  1.00 11.72           N
ATOM     18  H   ASN A   2      -7.646   3.716   3.069  1.00 15.02           H
ATOM     19  HA  ASN A   2      -6.434   1.341   3.195  1.00 14.10           H
ATOM     20  HB2 ASN A   2      -5.839   3.523   1.499  1.00 15.38           H
ATOM     21  HB3 ASN A   2      -5.249   2.058   1.323  1.00 15.38           H
ATOM     22 HD21 ASN A   2      -8.290   2.986  -0.623  1.00 11.72           H
ATOM     23 HD22 ASN A   2      -7.238   3.914  -0.123  1.00 11.72           H
ATOM     24  N   ASN A   3      -4.120   1.735   4.040  1.00 12.26           N
ATOM     25  CA  ASN A   3      -2.903   1.924   4.816  1.00 11.74           C
ATOM     26  C   ASN A   3      -1.741   1.254   4.100  1.00 11.10           C
ATOM     27  O   ASN A   3      -1.909   0.191   3.494  1.00 10.42           O
ATOM     28  CB  ASN A   3      -3.056   1.350   6.228  1.00 12.15           C
ATOM     29  CG  ASN A   3      -1.929   1.767   7.152  1.00 12.82           C
ATOM     30  OD1 ASN A   3      -1.912   2.886   7.664  1.00 15.05           O
ATOM     31  ND2 ASN A   3      -0.978   0.865   7.369  1.00 13.48           N
ATOM     32  H   ASN A   3      -4.169   0.962   3.666  1.00 12.26           H
ATOM     33  HA  ASN A   3      -2.710   2.872   4.891  1.00 11.74           H
ATOM     34  HB2 ASN A   3      -3.890   1.666   6.610  1.00 12.15           H
ATOM     35  HB3 ASN A   3      -3.060   0.382   6.177  1.00 12.15           H
ATOM     36 HD21 ASN A   3      -0.317   1.054   7.886  1.00 13.48           H
ATOM     37 HD22 ASN A   3      -1.023   0.093   6.993  1.00 13.48           H
ATOM     38  N   GLN A   4      -0.566   1.877   4.171  1.00 10.29           N
ATOM     39  CA  GLN A   4       0.648   1.356   3.551  1.00 10.53           C
ATOM     40  C   GLN A   4       1.840   2.098   4.138  1.00 10.24           C
ATOM     41  O   GLN A   4       2.636   2.685   3.398  1.00  8.86           O
ATOM     42  CB  GLN A   4       0.603   1.511   2.030  1.00  9.80           C
ATOM     43  CG  GLN A   4       1.704   0.759   1.298  1.00 10.25           C
ATOM     44  CD  GLN A   4       2.114   1.437   0.006  1.00 12.43           C
ATOM     45  OE1 GLN A   4       2.408   2.633  -0.015  1.00 14.62           O
ATOM     46  NE2 GLN A   4       2.138   0.675  -1.081  1.00  9.05           N
ATOM     47  H   GLN A   4      -0.445   2.622   4.584  1.00 10.29           H
ATOM     48  HA  GLN A   4       0.739   0.413   3.760  1.00 10.53           H
ATOM     49  HB2 GLN A   4      -0.248   1.175   1.709  1.00  9.80           H
ATOM     50  HB3 GLN A   4       0.693   2.452   1.811  1.00  9.80           H
ATOM     51  HG2 GLN A   4       2.486   0.705   1.870  1.00 10.25           H
ATOM     52  HG3 GLN A   4       1.389  -0.133   1.083  1.00 10.25           H
ATOM     53 HE21 GLN A   4       2.364   1.014  -1.839  1.00  9.05           H
ATOM     54 HE22 GLN A   4       1.927  -0.157  -1.028  1.00  9.05           H
ATOM     55  N   GLN A   5       1.961   2.074   5.462  1.00 10.38           N
ATOM     56  CA  GLN A   5       3.058   2.722   6.169  1.00 11.39           C
ATOM     57  C   GLN A   5       4.373   2.083   5.746  1.00 11.52           C
ATOM     58  O   GLN A   5       4.739   1.015   6.247  1.00 12.05           O
ATOM     59  CB  GLN A   5       2.876   2.609   7.684  1.00 11.96           C
ATOM     60  CG  GLN A   5       1.680   3.380   8.221  1.00 10.81           C
ATOM     61  CD  GLN A   5       1.282   2.942   9.616  1.00 13.10           C
ATOM     62  OE1 GLN A   5       0.766   1.842   9.811  1.00 10.65           O
ATOM     63  NE2 GLN A   5       1.521   3.804  10.598  1.00 12.30           N
ATOM     64  H   GLN A   5       1.405   1.678   5.985  1.00 10.38           H
ATOM     65  HA  GLN A   5       3.085   3.663   5.936  1.00 11.39           H
ATOM     66  HB2 GLN A   5       2.752   1.675   7.914  1.00 11.96           H
ATOM     67  HB3 GLN A   5       3.670   2.955   8.119  1.00 11.96           H
ATOM     68  HG2 GLN A   5       1.901   4.324   8.253  1.00 10.81           H
ATOM     69  HG3 GLN A   5       0.922   3.237   7.633  1.00 10.81           H
ATOM     70 HE21 GLN A   5       1.314   3.602  11.407  1.00 12.30           H
ATOM     71 HE22 GLN A   5       1.884   4.564  10.423  1.00 12.30           H
ATOM     72  N   ASN A   6       5.076   2.738   4.826  1.00 11.99           N
ATOM     73  CA  ASN A   6       6.341   2.244   4.285  1.00 12.30           C
ATOM     74  C   ASN A   6       7.438   2.541   5.303  1.00 13.40           C
ATOM     75  O   ASN A   6       8.034   3.623   5.287  1.00 13.92           O
ATOM     76  CB  ASN A   6       6.645   2.887   2.935  1.00 12.13           C
ATOM     77  CG  ASN A   6       5.482   2.788   1.967  1.00 12.77           C
ATOM     78  OD1 ASN A   6       5.274   1.755   1.332  1.00 14.27           O
ATOM     79  ND2 ASN A   6       4.715   3.866   1.852  1.00 10.07           N
ATOM     80  H   ASN A   6       4.834   3.491   4.489  1.00 11.99           H
ATOM     81  HA  ASN A   6       6.287   1.284   4.165  1.00 12.30           H
ATOM     82  HB2 ASN A   6       6.845   3.827   3.069  1.00 12.13           H
ATOM     83  HB3 ASN A   6       7.407   2.439   2.535  1.00 12.13           H
ATOM     84 HD21 ASN A   6       4.043   3.861   1.316  1.00 10.07           H
ATOM     85 HD22 ASN A   6       4.892   4.571   2.313  1.00 10.07           H
ATOM     86  N   TYR A   7       7.703   1.582   6.184  1.00 14.70           N
ATOM     87  CA  TYR A   7       8.747   1.737   7.190  1.00 15.18           C
ATOM     88  C   TYR A   7      10.119   1.442   6.594  1.00 15.91           C
ATOM     89  O   TYR A   7      10.230   0.840   5.526  1.00 15.76           O
ATOM     90  CB  TYR A   7       8.485   0.817   8.383  1.00 15.35           C
ATOM     91  CG  TYR A   7       7.187   1.099   9.104  1.00 14.45           C
ATOM     92  CD1 TYR A   7       7.095   2.127  10.032  1.00 14.80           C
ATOM     93  CD2 TYR A   7       6.054   0.335   8.859  1.00 15.68           C
ATOM     94  CE1 TYR A   7       5.910   2.387  10.695  1.00 14.33           C
ATOM     95  CE2 TYR A   7       4.865   0.587   9.516  1.00 13.46           C
ATOM     96  CZ  TYR A   7       4.799   1.615  10.433  1.00 15.09           C
ATOM     97  OH  TYR A   7       3.618   1.870  11.091  1.00 14.39           O
ATOM     98  OXT TYR A   7      11.149   1.800   7.166  1.00 17.49           O
ATOM     99  H   TYR A   7       7.290   0.829   6.221  1.00 14.70           H
ATOM    100  HA  TYR A   7       8.745   2.652   7.512  1.00 15.18           H
ATOM    101  HB2 TYR A   7       8.455  -0.100   8.068  1.00 15.35           H
ATOM    102  HB3 TYR A   7       9.207   0.924   9.023  1.00 15.35           H
ATOM    103  HD1 TYR A   7       7.843   2.650  10.211  1.00 14.80           H
ATOM    104  HD2 TYR A   7       6.096  -0.359   8.241  1.00 15.68           H
ATOM    105  HE1 TYR A   7       5.863   3.080  11.314  1.00 14.33           H
ATOM    106  HE2 TYR A   7       4.114   0.067   9.341  1.00 13.46           H
ATOM    107  HH  TYR A   7       3.716   2.518  11.617  1.00 14.39           H
TER
END
"""

def get_model(pdb_str, make_restraints=True):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  if(make_restraints):
    model.process(make_restraints=True)
  return model

def violations(model):
  model_ = model
  if(model.get_restraints_manager() is None):
    model_ = model.deep_copy()
    model_.process(make_restraints=True)
  return geometry_violations(
    geometry_restraints_manager = model_.get_restraints_manager().geometry,
    sites_cart                  = model_.get_sites_cart())

def shifts(model, sites_cart_start, selection):
  d = (model.get_sites_cart() - sites_cart_start).norms().select(selection)
  return flex.mean(d), flex.max(d)

def move_atom(model, atom_name, resseq, shift):
  sel = model.selection("name %s and resseq %d" % (atom_name, resseq))
  assert sel.count(True) == 1
  sites_cart = model.get_sites_cart()
  i = sel.iselection()[0]
  sites_cart[i] = tuple(x + s for x, s in zip(sites_cart[i], shift))
  model.set_sites_cart(sites_cart)

def residue_selection(model, resseqs):
  return model.selection("resseq %s" % " or resseq ".join(
    [str(r) for r in resseqs]))

def exercise_clean():
  """No violations: the model is not touched at all."""
  model = get_model(pdb_str)
  grm = model.get_restraints_manager()
  sites_cart_start = model.get_sites_cart().deep_copy()
  v = violations(model)
  assert len(v.bonds) == 0 and len(v.clashes) == 0, (v.bonds, v.clashes)
  r = model.regularize_geometry()
  assert not r.needed
  assert r.n_violations_start == 0
  assert r.n_violations_end == 0
  assert r.n_iterations == 0
  assert r.mean_shift == 0 and r.max_shift == 0
  assert (model.get_sites_cart() - sites_cart_start).norms().all_eq(0)
  assert model.get_restraints_manager() is grm

def perturb_sparse(model):
  # A stretched bond at the N-terminus and, at the C-terminus, OH of Tyr 7
  # pulled to 1.8 A from NE2 of Gln 5 (a clash and a stretched bond): the
  # peptide is a straight strand along a, the two ends are > 15 A apart.
  from scitbx import matrix
  move_atom(model, "C",  1, (0.30, 0.0, 0.0))
  sites_cart = model.get_sites_cart()
  oh = matrix.col(sites_cart[
    model.selection("resseq 7 and name OH").iselection()[0]])
  ne2 = matrix.col(sites_cart[
    model.selection("resseq 5 and name NE2").iselection()[0]])
  move_atom(model, "OH", 7, tuple(ne2 + (oh - ne2).normalize() * 1.8 - oh))

def exercise_sparse():
  """Few violations far apart: only their surroundings move."""
  model = get_model(pdb_str)
  grm = model.get_restraints_manager()
  perturb_sparse(model)
  sites_cart_start = model.get_sites_cart().deep_copy()
  v = violations(model)
  assert len(v.bonds) > 0, v.bonds
  assert len(v.clashes) > 0, v.clashes
  not_hd = ~model.get_hd_selection()
  r = model.regularize_geometry()
  assert r.needed
  assert r.n_violations_start == len(v.bonds) + len(v.clashes)
  assert r.n_violations_end == 0
  v = violations(model)
  assert len(v.bonds) == 0 and len(v.clashes) == 0, (v.bonds, v.clashes)
  assert r.max_shift <= 1.0 and r.mean_shift <= 0.3, (r.mean_shift, r.max_shift)
  # Residue 4 is more than 5 A from both violations: unchanged.
  far = residue_selection(model, [4]) & not_hd
  mean_far, max_far = shifts(model, sites_cart_start, far)
  assert max_far == 0, max_far
  # The ends moved.
  near = residue_selection(model, [1, 7]) & not_hd
  mean_near, max_near = shifts(model, sites_cart_start, near)
  assert max_near > 0.01, max_near
  # The restraints of the model are the ones it had.
  assert model.get_restraints_manager() is grm

def exercise_everywhere():
  """Violations throughout: the whole model moves."""
  random.seed(0)
  flex.set_random_seed(0)
  model = get_model(pdb_str)
  not_hd = ~model.get_hd_selection()
  sites_cart = model.get_sites_cart()
  n = sites_cart.size()
  noise = flex.vec3_double(flex.random_double(n*3)*0.3-0.15)
  noise.set_selected(~not_hd, (0.0, 0.0, 0.0))
  model.set_sites_cart(sites_cart + noise)
  sites_cart_start = model.get_sites_cart().deep_copy()
  v = violations(model)
  assert len(v.bonds) > 20, len(v.bonds)
  r = model.regularize_geometry()
  assert r.needed
  assert r.n_violations_end == 0
  v = violations(model)
  assert len(v.bonds) == 0 and len(v.clashes) == 0, (v.bonds, v.clashes)
  d = (model.get_sites_cart() - sites_cart_start).norms().select(not_hd)
  assert (d > 0.001).count(True) > 0.9 * d.size(), (d > 0.001).count(True)

def exercise_symmetry_clash():
  """Clashes with symmetry mates only (the whole model translated)."""
  model = get_model(pdb_str)
  model.set_sites_cart(model.get_sites_cart() + (-1.0, 0.0, 0.0))
  v = violations(model)
  assert len(v.bonds) == 0, v.bonds
  assert len(v.clashes) > 0
  for i, j, overlap, rt_mx in v.clashes:
    assert rt_mx is not None and not rt_mx.is_unit_mx(), rt_mx
  r = model.regularize_geometry()
  assert r.needed
  assert r.n_violations_start == len(v.clashes)
  assert r.n_violations_end == 0
  v = violations(model)
  assert len(v.bonds) == 0 and len(v.clashes) == 0, (v.bonds, v.clashes)

def exercise_mild():
  """Mild outliers (a 5 sigma bond, a 0.6 A overlap) are below the default
  thresholds: the model is not touched."""
  model = water_model(distance=2.2)
  move_atom(model, "C", 1, (0.12, 0.0, 0.0))
  sites_cart_start = model.get_sites_cart().deep_copy()
  grm = model.get_restraints_manager().geometry
  v = geometry_violations(grm, model.get_sites_cart(),
    bond_z_max=4.0, clash_overlap_max=0.4)
  assert len(v.bonds) == 2 and len(v.clashes) == 1, (v.bonds, v.clashes)
  assert max([z for i, j, z in v.bonds]) > 5, v.bonds
  assert v.clashes[0][2] > 0.5, v.clashes
  v = geometry_violations(grm, model.get_sites_cart(),
    bond_z_max=6.0, clash_overlap_max=1.0)
  assert len(v.bonds) == 0 and len(v.clashes) == 0, (v.bonds, v.clashes)
  r = model.regularize_geometry()
  assert not r.needed
  assert r.n_violations_start == 0
  assert (model.get_sites_cart() - sites_cart_start).norms().all_eq(0)

def exercise_force():
  """force=True minimizes a clean model within the shift budget."""
  model = get_model(pdb_str)
  sites_cart_start = model.get_sites_cart().deep_copy()
  r = model.regularize_geometry(force=True)
  assert not r.needed
  assert r.n_violations_start == 0
  assert r.n_violations_end == 0
  assert r.n_iterations > 0
  assert r.max_shift <= 1.0 and r.mean_shift <= 0.3, (r.mean_shift, r.max_shift)
  not_hd = ~model.get_hd_selection()
  mean_shift, max_shift = shifts(model, sites_cart_start, not_hd)
  assert max_shift > 0

def exercise_no_h():
  """Works without H atoms (no riding H)."""
  lines = [l for l in pdb_str.splitlines() if not l.endswith(" H")]
  model = get_model("\n".join(lines))
  assert not model.has_hd()
  perturb_sparse(model)
  v = violations(model)
  assert len(v.bonds) > 0
  r = model.regularize_geometry()
  assert r.needed
  assert r.n_violations_end == 0
  v = violations(model)
  assert len(v.bonds) == 0 and len(v.clashes) == 0, (v.bonds, v.clashes)

def water_model(distance=1.7):
  """The clean model in a P 1 box (no symmetry contacts) plus a water whose
  O is distance (A) from OE1 of Gln 5, a clash with the water O."""
  from scitbx import matrix
  pdb_str_box = pdb_str.replace(pdb_str.splitlines()[0],
    "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1")
  model = get_model(pdb_str_box, make_restraints=False)
  sites_cart = model.get_sites_cart()
  oe1 = matrix.col(sites_cart[
    model.selection("resseq 5 and name OE1").iselection()[0]])
  cd = matrix.col(sites_cart[
    model.selection("resseq 5 and name CD").iselection()[0]])
  d = (oe1 - cd).normalize()
  o = oe1 + distance * d
  h1 = o + 0.85 * d
  h2 = o + 0.85 * d.ortho().normalize()
  lines = []
  for name, site in [("O ", o), ("H1", h1), ("H2", h2)]:
    lines.append("HETATM  %3d  %-3s HOH A 201    %8.3f%8.3f%8.3f  1.00 20.00"
      "           %s" % (108 + len(lines), name, site[0], site[1], site[2],
      name[0]))
  return get_model(pdb_str_box.replace("END", "\n".join(lines) + "\nEND"))

def exercise_water():
  """Water O never moves (it would fly away), water H do: the clash with
  the water is resolved by the protein."""
  model = water_model()
  sites_cart_start = model.get_sites_cart().deep_copy()
  water_o = model.selection("water and not element H")
  assert water_o.count(True) == 1
  v = violations(model)
  assert len(v.clashes) > 0, v.clashes
  r = model.regularize_geometry()
  assert r.needed
  assert r.n_violations_end == 0
  v = violations(model)
  assert len(v.bonds) == 0 and len(v.clashes) == 0, (v.bonds, v.clashes)
  d = (model.get_sites_cart() - sites_cart_start).norms()
  assert d.select(water_o).all_eq(0)
  protein = ~model.selection("water") & ~model.get_hd_selection()
  assert flex.max(d.select(protein)) > 0.01

def exercise_no_restraints():
  """A model without restraints is regularized, and left without them."""
  model = get_model(pdb_str, make_restraints=False)
  assert model.get_restraints_manager() is None
  perturb_sparse(model)
  r = model.regularize_geometry()
  assert r.needed
  assert r.n_violations_end == 0
  assert model.get_restraints_manager() is None
  v = violations(model)
  assert len(v.bonds) == 0 and len(v.clashes) == 0, (v.bonds, v.clashes)

if (__name__ == "__main__"):
  exercise_clean()
  exercise_sparse()
  exercise_everywhere()
  exercise_symmetry_clash()
  exercise_mild()
  exercise_force()
  exercise_no_h()
  exercise_water()
  exercise_no_restraints()
  print("OK")
