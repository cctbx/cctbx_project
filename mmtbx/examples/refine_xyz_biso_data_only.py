from __future__ import absolute_import, division, print_function
import iotbx.pdb
import random
import mmtbx.f_model
from scitbx.array_family import flex
from libtbx import group_args
from cctbx import xray
import scitbx.lbfgs
from cctbx import adptbx
from six.moves import range

pdb_str = """
CRYST1    5.827    9.541    6.239  90.00  90.00  90.00 P 1
ATOM      1 OH   TYR     1       3.566   1.332   4.100  1.00 10.00           O
ATOM      2 CZ   TYR     1       3.227   2.567   3.610  1.00 10.00           C
ATOM      3 CE1  TYR     1       3.773   3.048   2.423  1.00 10.00           C
ATOM      4 CE2  TYR     1       2.305   3.332   4.316  1.00 10.00           C
ATOM      5 CD1  TYR     1       3.395   4.305   1.954  1.00 10.00           C
ATOM      6 CD2  TYR     1       1.920   4.578   3.822  1.00 10.00           C
ATOM      7 CG   TYR     1       2.457   5.086   2.634  1.00 10.00           C
ATOM      8 CB   TYR     1       2.033   6.443   2.116  1.00 10.00           C
ATOM      9 CA   TYR     1       2.859   7.627   2.655  1.00 10.00           C
ATOM     10 N    TYR     1       4.203   7.676   2.028  1.00 10.00           N
ATOM     11 HN1  TYR     1       4.827   8.410   2.401  1.00 10.00           H
ATOM     12 HN2  TYR     1       4.763   6.829   2.214  1.00 10.00           H
ATOM     13 HN3  TYR     1       4.148   7.735   1.000  1.00 10.00           H
ATOM     14 HH   TYR     1       4.395   1.000   3.753  1.00 10.00           H
ATOM     15 HD1  TYR     1       3.841   4.670   1.051  1.00 10.00           H
ATOM     16 HE1  TYR     1       4.506   2.464   1.906  1.00 10.00           H
ATOM     17 HD2  TYR     1       1.232   5.163   4.402  1.00 10.00           H
ATOM     18 HE2  TYR     1       1.940   2.925   5.239  1.00 10.00           H
ATOM     19 HB1  TYR     1       2.121   6.457   1.034  1.00 10.00           H
ATOM     20 HB2  TYR     1       1.000   6.608   2.385  1.00 10.00           H
ATOM     21 HA   TYR     1       2.327   8.541   2.425  1.00 10.00           H
ATOM     22 O    TYR     1       4.138   7.163   4.634  1.00 10.00           O
ATOM     23 C    TYR     1       3.064   7.506   4.165  1.00 10.00           C
TER
END
"""

class minimizer(object):
  def __init__(self,
        fmodel,
        max_iterations=100,
        sites = False,
        u_iso = False):
    self.fmodel = fmodel
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.x_target_functor = self.fmodel.target_functor()
    self.sites = sites
    self.u_iso = u_iso
    if(self.sites):
      self.x = self.fmodel.xray_structure.sites_cart().as_double()
    if(self.u_iso):
      assert self.fmodel.xray_structure.scatterers().size() == \
        self.fmodel.xray_structure.use_u_iso().count(True)
      self.x = self.fmodel.xray_structure.extract_u_iso_or_u_equiv()
    if(self.sites):
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        site       = True)
    if(self.u_iso):
      sel = flex.bool(
        self.fmodel.xray_structure.scatterers().size(), True).iselection()
      self.fmodel.xray_structure.scatterers().flags_set_grad_u_iso(
        iselection = sel)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))
    self.fmodel.xray_structure.tidy_us()
    self.fmodel.xray_structure.apply_symmetry_sites()
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)

  def compute_functional_and_gradients(self):
    if(self.sites):
      self.fmodel.xray_structure.set_sites_cart(
        sites_cart = flex.vec3_double(self.x))
    if(self.u_iso):
      self.fmodel.xray_structure.set_u_iso(values = self.x)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)
    tgx = self.x_target_functor(compute_gradients=True)
    if(self.sites):
      tx = tgx.target_work()
      gx = flex.vec3_double(tgx.\
        gradients_wrt_atomic_parameters(site=True).packed())
      f = tx
      g = gx
    if(self.u_iso):
      tx = tgx.target_work()
      gx = tgx.gradients_wrt_atomic_parameters(u_iso=True)
      f = tx
      g = gx
    return f, g.as_double()

def get_inputs(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  ph = pdb_inp.construct_hierarchy()
  xrs = ph.extract_xray_structure(crystal_symmetry =
    pdb_inp.crystal_symmetry())
  return group_args(pdb_hierarchy = ph, xray_structure = xrs)

def run(refine_xyz=True, refine_adp=False):
  # get xray_structure from PDB file
  inp = get_inputs(pdb_str = pdb_str)
  if(1):
    inp.pdb_hierarchy.adopt_xray_structure(inp.xray_structure)
    inp.pdb_hierarchy.write_pdb_file(file_name="start.pdb")
  # simulate poor starting model: shake coordinates and B-factors
  xrs_poor = inp.xray_structure.deep_copy_scatterers()
  if(refine_xyz):
    xrs_poor.shake_sites_in_place(mean_distance = 0.3)
  if(refine_adp):
    for sc in xrs_poor.scatterers():
      sc.u_iso = adptbx.b_as_u(random.choice([7,8,9,11,12,13]))
  if(1):
    inp.pdb_hierarchy.adopt_xray_structure(xrs_poor)
    inp.pdb_hierarchy.write_pdb_file(file_name="poor.pdb")
  # simulate Fobs
  f_obs = abs(inp.xray_structure.structure_factors(
    d_min = 1.0,
    algorithm="direct").f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  # get fmodel
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = xrs_poor,
    sf_and_grads_accuracy_params = params,
    target_name                  = "ls_wunit_kunit")
  # refinement loop
  print("start r_factor: %6.4f" % fmodel.r_work())
  for macro_cycle in range(100):
    # refine coordinates
    if(refine_xyz):
      minimized = minimizer(fmodel = fmodel, sites = True)
      print("  macro_cycle %3d (sites) r_factor: %6.4f"%(macro_cycle,
        fmodel.r_work()))
    # refine ADPs
    if(refine_adp):
      minimized = minimizer(fmodel = fmodel, u_iso = True)
      print("  macro_cycle %3d (adp)   r_factor: %6.4f"%(macro_cycle, fmodel.r_work()))
  if(1):
    inp.pdb_hierarchy.adopt_xray_structure(fmodel.xray_structure)
    inp.pdb_hierarchy.write_pdb_file(file_name="refined.pdb")

if (__name__ == "__main__"):
  run()
