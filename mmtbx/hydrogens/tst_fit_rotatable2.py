"""
fit_rotatable2: one OMIT map for the whole model (all rotatable H omitted at
once), every residue fitted in it.
"""
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
import mmtbx.hydrogens
from mmtbx.hydrogens import map_manager
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from scitbx.matrix import rotate_point_around_axis, col

def dist(a, b):
  return (col(a) - col(b)).length()

# The first six residues of 2lvr (NMR, with H): MET LYS PRO TYR VAL CYS.
pdb_str = """\
CRYST1   51.084   33.664   29.435  90.00  90.00  90.00 P 1
ATOM      1  N   MET A   1      39.628  19.616  12.383  1.00  4.50           N
ATOM      2  CA  MET A   1      40.082  18.496  13.243  1.00  3.96           C
ATOM      3  C   MET A   1      38.919  17.903  14.026  1.00  2.83           C
ATOM      4  O   MET A   1      38.503  16.772  13.771  1.00  3.15           O
ATOM      5  CB  MET A   1      41.174  18.951  14.218  1.00  4.62           C
ATOM      6  CG  MET A   1      42.527  19.192  13.574  1.00  5.46           C
ATOM      7  SD  MET A   1      43.782  19.679  14.775  1.00  6.29           S
ATOM      8  CE  MET A   1      45.227  19.787  13.724  1.00  7.14           C
ATOM      9  H1  MET A   1      40.440  20.057  11.898  1.00  4.66           H
ATOM     10  H2  MET A   1      39.145  20.338  12.953  1.00  4.62           H
ATOM     11  H3  MET A   1      38.963  19.260  11.660  1.00  5.04           H
ATOM     12  HA  MET A   1      40.488  17.727  12.601  1.00  4.40           H
ATOM     13  HB2 MET A   1      40.858  19.871  14.687  1.00  4.77           H
ATOM     14  HB3 MET A   1      41.294  18.195  14.981  1.00  4.83           H
ATOM     15  HG2 MET A   1      42.845  18.281  13.089  1.00  5.77           H
ATOM     16  HG3 MET A   1      42.427  19.975  12.839  1.00  5.67           H
ATOM     17  HE1 MET A   1      45.405  18.828  13.260  1.00  7.30           H
ATOM     18  HE2 MET A   1      46.084  20.064  14.318  1.00  7.44           H
ATOM     19  HE3 MET A   1      45.061  20.532  12.961  1.00  7.52           H
ATOM     20  N   LYS A   2      38.392  18.662  14.982  1.00  2.07           N
ATOM     21  CA  LYS A   2      37.311  18.171  15.822  1.00  1.44           C
ATOM     22  C   LYS A   2      35.954  18.619  15.284  1.00  1.07           C
ATOM     23  O   LYS A   2      35.679  19.814  15.157  1.00  1.15           O
ATOM     24  CB  LYS A   2      37.494  18.631  17.276  1.00  2.11           C
ATOM     25  CG  LYS A   2      37.436  20.140  17.464  1.00  2.75           C
ATOM     26  CD  LYS A   2      37.625  20.532  18.918  1.00  3.31           C
ATOM     27  CE  LYS A   2      37.438  22.028  19.121  1.00  4.33           C
ATOM     28  NZ  LYS A   2      36.034  22.458  18.874  1.00  4.97           N
ATOM     29  H   LYS A   2      38.733  19.575  15.122  1.00  2.52           H
ATOM     30  HA  LYS A   2      37.345  17.092  15.794  1.00  1.78           H
ATOM     31  HB2 LYS A   2      36.716  18.189  17.881  1.00  2.38           H
ATOM     32  HB3 LYS A   2      38.453  18.285  17.630  1.00  2.79           H
ATOM     33  HG2 LYS A   2      38.216  20.596  16.876  1.00  3.31           H
ATOM     34  HG3 LYS A   2      36.473  20.497  17.127  1.00  2.97           H
ATOM     35  HD2 LYS A   2      36.900  20.004  19.520  1.00  3.33           H
ATOM     36  HD3 LYS A   2      38.622  20.257  19.228  1.00  3.54           H
ATOM     37  HE2 LYS A   2      37.706  22.275  20.138  1.00  4.70           H
ATOM     38  HE3 LYS A   2      38.091  22.552  18.440  1.00  4.68           H
ATOM     39  HZ1 LYS A   2      35.387  21.966  19.526  1.00  5.32           H
ATOM     40  HZ2 LYS A   2      35.752  22.240  17.892  1.00  5.02           H
ATOM     41  HZ3 LYS A   2      35.944  23.485  19.024  1.00  5.43           H
ATOM     42  N   PRO A   3      35.107  17.662  14.900  1.00  0.85           N
ATOM     43  CA  PRO A   3      33.736  17.933  14.509  1.00  0.57           C
ATOM     44  C   PRO A   3      32.772  17.826  15.685  1.00  0.42           C
ATOM     45  O   PRO A   3      33.168  17.505  16.809  1.00  0.58           O
ATOM     46  CB  PRO A   3      33.470  16.821  13.505  1.00  0.65           C
ATOM     47  CG  PRO A   3      34.236  15.659  14.039  1.00  1.14           C
ATOM     48  CD  PRO A   3      35.419  16.226  14.793  1.00  1.08           C
ATOM     49  HA  PRO A   3      33.631  18.895  14.029  1.00  0.64           H
ATOM     50  HB2 PRO A   3      32.408  16.613  13.463  1.00  0.67           H
ATOM     51  HB3 PRO A   3      33.830  17.114  12.528  1.00  0.88           H
ATOM     52  HG2 PRO A   3      33.609  15.085  14.707  1.00  1.43           H
ATOM     53  HG3 PRO A   3      34.571  15.044  13.222  1.00  1.54           H
ATOM     54  HD2 PRO A   3      35.492  15.773  15.771  1.00  1.26           H
ATOM     55  HD3 PRO A   3      36.332  16.070  14.235  1.00  1.20           H
ATOM     56  N   TYR A   4      31.513  18.093  15.411  1.00  0.26           N
ATOM     57  CA  TYR A   4      30.448  17.858  16.364  1.00  0.19           C
ATOM     58  C   TYR A   4      29.562  16.750  15.828  1.00  0.15           C
ATOM     59  O   TYR A   4      28.773  16.964  14.911  1.00  0.21           O
ATOM     60  CB  TYR A   4      29.640  19.135  16.624  1.00  0.27           C
ATOM     61  CG  TYR A   4      30.368  20.120  17.509  1.00  0.46           C
ATOM     62  CD1 TYR A   4      31.246  21.052  16.976  1.00  0.56           C
ATOM     63  CD2 TYR A   4      30.184  20.105  18.885  1.00  0.68           C
ATOM     64  CE1 TYR A   4      31.921  21.942  17.787  1.00  0.76           C
ATOM     65  CE2 TYR A   4      30.854  20.994  19.703  1.00  0.88           C
ATOM     66  CZ  TYR A   4      31.720  21.911  19.148  1.00  0.91           C
ATOM     67  OH  TYR A   4      32.395  22.797  19.956  1.00  1.11           O
ATOM     68  H   TYR A   4      31.284  18.446  14.524  1.00  0.36           H
ATOM     69  HA  TYR A   4      30.898  17.528  17.289  1.00  0.25           H
ATOM     70  HB2 TYR A   4      29.411  19.626  15.682  1.00  0.27           H
ATOM     71  HB3 TYR A   4      28.718  18.868  17.117  1.00  0.37           H
ATOM     72  HD1 TYR A   4      31.401  21.076  15.908  1.00  0.57           H
ATOM     73  HD2 TYR A   4      29.505  19.385  19.316  1.00  0.76           H
ATOM     74  HE1 TYR A   4      32.600  22.659  17.352  1.00  0.89           H
ATOM     75  HE2 TYR A   4      30.696  20.970  20.770  1.00  1.07           H
ATOM     76  HH  TYR A   4      32.824  22.312  20.675  1.00  1.47           H
ATOM     77  N   VAL A   5      29.732  15.559  16.372  1.00  0.22           N
ATOM     78  CA  VAL A   5      29.068  14.387  15.841  1.00  0.18           C
ATOM     79  C   VAL A   5      27.718  14.159  16.501  1.00  0.14           C
ATOM     80  O   VAL A   5      27.535  14.420  17.691  1.00  0.19           O
ATOM     81  CB  VAL A   5      29.939  13.124  16.000  1.00  0.25           C
ATOM     82  CG1 VAL A   5      31.305  13.341  15.370  1.00  0.35           C
ATOM     83  CG2 VAL A   5      30.078  12.725  17.464  1.00  0.30           C
ATOM     84  H   VAL A   5      30.318  15.464  17.156  1.00  0.34           H
ATOM     85  HA  VAL A   5      28.911  14.550  14.783  1.00  0.17           H
ATOM     86  HB  VAL A   5      29.453  12.321  15.472  1.00  0.25           H
ATOM     87 HG11 VAL A   5      31.820  14.132  15.893  1.00  1.09           H
ATOM     88 HG12 VAL A   5      31.180  13.619  14.335  1.00  1.02           H
ATOM     89 HG13 VAL A   5      31.881  12.430  15.433  1.00  1.04           H
ATOM     90 HG21 VAL A   5      30.692  11.839  17.542  1.00  1.09           H
ATOM     91 HG22 VAL A   5      29.099  12.523  17.876  1.00  1.04           H
ATOM     92 HG23 VAL A   5      30.541  13.531  18.014  1.00  1.00           H
ATOM     93  N   CYS A   6      26.770  13.704  15.702  1.00  0.12           N
ATOM     94  CA  CYS A   6      25.453  13.343  16.196  1.00  0.13           C
ATOM     95  C   CYS A   6      25.544  12.108  17.072  1.00  0.16           C
ATOM     96  O   CYS A   6      25.970  11.051  16.615  1.00  0.20           O
ATOM     97  CB  CYS A   6      24.516  13.057  15.032  1.00  0.13           C
ATOM     98  SG  CYS A   6      22.808  12.741  15.514  1.00  0.21           S
ATOM     99  H   CYS A   6      26.960  13.626  14.742  1.00  0.16           H
ATOM    100  HA  CYS A   6      25.067  14.168  16.777  1.00  0.17           H
ATOM    101  HB2 CYS A   6      24.512  13.898  14.365  1.00  0.19           H
ATOM    102  HB3 CYS A   6      24.871  12.183  14.502  1.00  0.14           H
END
"""

def get_model(lines):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=lines.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  model.setup_scattering_dictionaries(scattering_table="n_gaussian")
  return model

def get_fmodel(model, sites_cart):
  """fmodel with f_obs = |Fcalc| of model and the structure at sites_cart."""
  xrs = model.get_xray_structure().deep_copy_scatterers()
  f_obs = abs(xrs.structure_factors(d_min=1.0, algorithm="direct").f_calc())
  flags = f_obs.generate_r_free_flags(fraction=0.1)
  xrs.set_sites_cart(sites_cart)
  return mmtbx.f_model.manager(f_obs=f_obs, r_free_flags=flags,
    xray_structure=xrs)

def i_seq(model, resname, name):
  found = [a.i_seq for a in model.get_hierarchy().atoms()
    if a.parent().resname == resname and a.name.strip() == name]
  assert len(found) == 1, (resname, name, found)
  return found[0]

def rotated_tyr_hh(model, angle):
  """Sites of model with the TYR HH rotated by angle about CZ-OH."""
  sites = model.get_sites_cart().deep_copy()
  cz, oh, hh = [i_seq(model, "TYR", n) for n in ("CZ", "OH", "HH")]
  sites[hh] = rotate_point_around_axis(axis_point_1=sites[cz],
    axis_point_2=sites[oh], point=sites[hh], angle=angle, deg=True)
  return sites

class counting(object):
  """Counts map_manager.update_omit_map calls and records the selections."""
  def __init__(self):
    self.selections = []
    self.orig = map_manager.update_omit_map
  def __enter__(self):
    def wrapped(mm, omit_selection):
      self.selections.append(list(omit_selection))
      return self.orig(mm, omit_selection)
    map_manager.update_omit_map = wrapped
    return self
  def __exit__(self, *args):
    map_manager.update_omit_map = self.orig

def exercise_one_map():
  """One OMIT map, with all 16 rotatable H of the fixture omitted."""
  model = get_model(pdb_str)
  fmodel = get_fmodel(model, rotated_tyr_hh(model, 120))
  with counting() as c:
    mmtbx.hydrogens.fit_rotatable2(model=model, fmodel=fmodel, reset_h=False)
  assert len(c.selections) == 1, len(c.selections)
  assert len(c.selections[0]) == 16, len(c.selections[0])
  assert len(set(c.selections[0])) == 16

def exercise_recovers_rotation():
  """The TYR HH rotated by 120 degrees is fitted back to the true site."""
  model = get_model(pdb_str)
  true = model.get_sites_cart().deep_copy()
  rotated = rotated_tyr_hh(model, 120)
  hh = i_seq(model, "TYR", "HH")
  assert dist(true[hh], rotated[hh]) > 1.0
  fmodel = get_fmodel(model, rotated) # f_obs from the true sites
  model.set_sites_cart(rotated)
  mmtbx.hydrogens.fit_rotatable2(model=model, fmodel=fmodel, reset_h=False)
  fitted = model.get_sites_cart()
  assert dist(fitted[hh], true[hh]) < 0.15, dist(fitted[hh], true[hh])
  # Only rotatable H move; the others by the 1-degree scan noise at most.
  rotatable = set()
  for rn, names in [("MET", ["H1","H2","H3","HE1","HE2","HE3"]),
                    ("LYS", ["HZ1","HZ2","HZ3"]), ("TYR", ["HH"]),
                    ("VAL", ["HG11","HG12","HG13","HG21","HG22","HG23"])]:
    rotatable.update(i_seq(model, rn, n) for n in names)
  assert len(rotatable) == 16
  moved = [i for i in range(true.size()) if dist(fitted[i], rotated[i]) > 1e-6]
  assert set(moved).issubset(rotatable), sorted(set(moved) - rotatable)
  big = [i for i in moved if dist(fitted[i], rotated[i]) > 0.1]
  assert big == [hh], big
  assert approx_equal(list(fmodel.xray_structure.sites_cart().as_double()),
    list(fitted.as_double()), 1.e-9)

def exercise_no_rotatable():
  """A model without rotatable H: no map, sites unchanged."""
  lines = [l for l in pdb_str.split("\n")
    if l.startswith("CRYST1") or (l.startswith("ATOM") and l[17:20] == "PRO")]
  model = get_model("\n".join(lines + ["END"]))
  assert model.size() > 10
  sites = model.get_sites_cart().deep_copy()
  fmodel = get_fmodel(model, sites)
  with counting() as c:
    mmtbx.hydrogens.fit_rotatable2(model=model, fmodel=fmodel, reset_h=False)
  assert len(c.selections) == 0, len(c.selections)
  assert approx_equal(list(model.get_sites_cart().as_double()),
    list(sites.as_double()), 1.e-9)

def run():
  exercise_one_map()
  exercise_recovers_rotation()
  exercise_no_rotatable()

if(__name__ == "__main__"):
  run()
  print("OK")
